function [incfac,pdrift_DI,drift_DI,def_based_DI,maxDisplacement]=...
        ElastoPlasticPushoverPlaneFrames(qbary,a,Mp,nbars,nnodes,...
        e,inertia,coordxy,ni,nf,len,co,se,Edof,support,edof,ndof,...
        rxbar,rybar,mbar,seismic_forces,h_floor,dof_forces,dof_disp)

%------------------------------------------------------------------------
% [incfac,pdrift_DI,drift_DI,def_based_DI,maxDisplacement]=...
% ElastoPlasticPushoverPlaneFrames(qbary,a,Mp,nbars,nnodes,...
% e,inertia,coordxy,ni,nf,len,co,se,Edof,support,edof,ndof,...
% rxbar,rybar,mbar,seismic_forces,h_floor,dof_forces,dof_disp)
%
%------------------------------------------------------------------------
% PURPOSE
%  To compute a static non-linear pushover analysis of a plane frame
%  
% 
% INPUT:  a = [area_bar;
%               ...]                 area of all elements
%
%         Mp = [Mpi Mpj;             Plastic Moment for each member 
%               ... ]                (i) initial node, (j) final node
%
%         rxbar = [bar, Fxi, Fxj]    Reacting horizontal forces by 
%                                    external loads
%
%         rybar = [bar, Fyi, Fyj]    Reacting vertical forces by 
%                                    external loads
%
%         mbar = [bar, mi, mj]       Reacting moments by external loads
%
%         nbars, nnodes              number of bars and number of nodes
%
%         e = [e_bar;                Elasticity modulus of each element
%               ...]                    
%
%         inertia = [inertia_bar;    in-plane inertia for all elements'
%                       ...]         cross-section
%                                    
%         coordxy = [coordx coordy;  node coordinates for all nodes
%                       ...];
%
%         ni                         list of initial nodes of all bars,
%         nf                         list of final nodes of all bars:
%                                         size = [nbars,1]
%
%         len, co, se                len: length of each bar:
%                                    co: cos direction of each bar
%                                    se: sen direction of each bar:
%                                        size = [nbar,1]
%
%         Edof                       Degrees of freedom of each bar
%                                    according to its nodes-orientation:
%                                    size = [nbar,7]
% 
%         qbary = [bar, load;
%                   ..    .. ]       uniform distributed load acting 
%                                    downwards (only for beams)
%
%         support = [i, j]           support at each bar's end
%                                    options: "Art" or "Fixed"
%                                    (i) initial node, (j) final node
%
%         edof                       non-restricted dof
%         ndof                       number of dof
%
%         seismic_forces = [f(1);]   lateral forces per floor:
%                          f(n);]    size = [nfloors,1]
%
%         h_floor = [h(1);           Height of each floor from bottom
%                    h(n)]           to top: size = [nfloors,1]
%
%         dof_forces = [dof-f(1),    dof at which the lateral forces are
%                       dof-f(n)]    applied (from bottom to top) - global
%          
%         dof_disp = [dof-f(1),      dof at which lateral forces are
%                     dof-f(n)]      applied in the reduced system of
%                                    dof reference (edof)
%
% OUTPUT: incfac                     incremental load factor:
%                                    Collapse Safety Factor (CSF)
%
%         pdrift_DI                  Plastic inter-story drift Damage 
%                                    Index per floor: size = [nfloors,1]
%
%         drift_DI                   Inter-story drift Damage 
%                                    Index per floor: size = [nfloors,1]
%
%         def_based_DI               Deformation based Damage Index
%                                    per floor size = [nfloors,1]
%
%         maxDisplacement            Max absoloute lateral displacement
%                                    for each floor: size = [nfloors,1]
% 
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2022-01-18
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------
nfloors=length(h_floor);

displacement_history_floor_left=[];
force_floor_left=[];

displacement_history_floor_right=[];
force_floor_right=[];

zero_disp=zeros(nfloors,1);
displacement_history_floor_left=[displacement_history_floor_left,zero_disp];
force_floor_left=[force_floor_left,zero_disp];

displacement_history_floor_right=[displacement_history_floor_right,zero_disp];
force_floor_right=[force_floor_right,zero_disp];

plast_bars=zeros(2,nbars);
bar_plas_node=[];
iter_collection=[];
lambda_collection=[];
looping=0;
incfac=1.0;
iteration=0;
while looping==0
    Kglobal=zeros(3*nnodes);
    fglobal=zeros(3*nnodes,1);

    for i=1:nfloors
        fglobal(dof_forces(i))=seismic_forces(i)*incfac;
    end

    elemental_matrices=zeros(6*nbars,6);

     ep_bars=zeros(nbars,3); 
     eq_bars=zeros(nbars,2);
     for i=1:nbars      

         ex=[coordxy(ni(i),1) coordxy(nf(i),1)];
         ey=[coordxy(ni(i),2) coordxy(nf(i),2)];
         ep=[e(i) a(i) inertia(i)];
         eq=[0 -qbary(i,2)];

         ep_bars(i,:)=ep;
         eq_bars(i,:)=eq;
         [Ke_barra,fe_barra]=beam2e(ex,ey,ep,eq); % This is a CALFEM
                                                  % Download at: 
                                                  % https://www.byggmek.lth.se/english/calfem/

         if support(i,2)=="Fixed" && support(i,3)=="Art"

            ak=(e(i)*a(i)/len(i)*co(i)^2)+(3*e(i)*inertia(i)/len(i)^3*se(i)^2);
            b=e(i)*a(i)/len(i)*se(i)^2+3*e(i)*inertia(i)/len(i)^3*co(i)^2;    
            ck=(e(i)*a(i)/len(i)-3*e(i)*inertia(i)/len(i)^3)*se(i)*co(i);
            d=-3*e(i)*inertia(i)*se(i)/len(i)^2;

            ek=3*e(i)*inertia(i)/len(i)^2*co(i);
            f=3*e(i)*inertia(i)/len(i);

            Ke_barra=[ak ck d -ak -ck 0;
                      ck b ek -ck -b 0;
                      d ek f -d -ek 0;
                      -ak -ck -d ak ck 0;
                      -ck -b -ek ck b 0;
                      0 0 0 0 0 0];
         end
         elemental_matrices((i-1)*6+1:6*i,:)=Ke_barra; %% guardando Ke_barra
         [Kglobal,fglobal]=assem(Edof(i,:),Kglobal,Ke_barra,fglobal,fe_barra);

     end     

    globalKreduced=Kglobal(edof,edof);
    if iteration==0
        det_Kred=det(globalKreduced);
    else
        det_Kred_post=det(globalKreduced);
    end
    
    %%%%%%%%--------- fb (cuando hay cargas especiales entre barras)
    fbreduced=convectiveForceVector(rxbar,rybar,mbar,nbars,nnodes,...
                                    edof,ni,nf);

    %%%%%%%%--------- fd (cargas externas en nodos)---------%%%%%%%%
    
    ext_force_nodes=fglobal(edof);
    %%%%%%%---------- fb + fd (reducidos)-----------%%%%%%%
    forces=ext_force_nodes+fbreduced;

    [displacements,r]=solveq(globalKreduced,forces);

    %%%%%%%%--------displacements---------%%%%%%%%
    
    u=displacements;
    
    %%%%%%%%%%%------Global displacement vector--------%%%%%%%
    Uglobal=zeros(3*nnodes,1);
    Uglobal(edof)=u;
    
    %Ed=extract(Edof,Uglobal);
    ex=coordxy(:,1);
    ey=coordxy(:,2);

    Ex=zeros(nbars,2);
    Ey=zeros(nbars,2);

    for j=1:nbars
        Ex(j,1)=ex(Edof(j,4)/3);
        Ex(j,2)=ex(Edof(j,7)/3);

        Ey(j,1)=ey(Edof(j,4)/3);
        Ey(j,2)=ey(Edof(j,7)/3);

    end 

    reac_collection=zeros(6,nbars);
    % --- computation of mechanic elements at the end of bars --- %
    for i=1:nbars
        ue=Uglobal(Edof(i,2:7));
        ke=elemental_matrices((i-1)*6+1:6*i,:);
        fe=elementalForces(ke,ue,rxbar,rybar,mbar,i);
        
        reac_collection(:,i)=fe;

        rt=[co(i) se(i) 0;
            -se(i) co(i) 0;
            0 0 1];

        fe_local_ni=rt*fe(1:3);
        fe_local_nf=rt*fe(4:6);

        fe_local=[fe_local_ni'; fe_local_nf'];
        
    end


    iteration=iteration+1;
    current_plas=0;

    plastified_bars=[];
    for i=1:nbars

        % Detect if any end of this bar (i) has been plastified___________
        % ________________________________________________________________
        bar_plas_check=0;
        if isempty(bar_plas_node)==0
            for j=1:length(bar_plas_node(:,1))
                if i==bar_plas_node(j,1)
                    bar_plas_check=bar_plas_check+1;

                end
            end

        end

        plastified_bars=[plastified_bars;
                             i,bar_plas_check];

        if bar_plas_check==1
            % Detect if the other end has been also plastified
            if abs(reac_collection(3,i))>=Mp(i,1)

                current_plas=1;
                plast_bars(1,i)=1;
                plast_bars(2,i)=1;

                inertia(i)=1e-10;
                % Change condition Fixed-Art to Art-Art
                support(i,2)="Art";
                support(i,3)="Art";

                mplas=reac_collection(3,i);
                mbar(i,2)=mbar(i,2)+0.5*mplas;

                bar_plas_node=[bar_plas_node;
                                i,1,incfac];
            end

        elseif bar_plas_check==0

            if abs(reac_collection(3,i))>=Mp(i,1) && ...
                    abs(reac_collection(6,i))<Mp(i,2)

                current_plas=1;
                mplas=reac_collection(3,i);
                plast_bars(1,i)=1;

                % Change node condition Fixed-Art....
                node=ni(i);
                ni(i)=nf(i);
                nf(i)=node;

                % Alternate plastification moment
                mp_ext=Mp(i,1);
                Mp(i,1)=Mp(i,2);
                Mp(i,2)=mp_ext;

                % change condition to Fixed-Art
                support(i,2)="Fixed";
                support(i,3)="Art";

                Edof(i,2)=3*ni(i)-2;
                Edof(i,3)=3*ni(i)-1;
                Edof(i,4)=3*ni(i);

                Edof(i,5)=3*nf(i)-2;
                Edof(i,6)=3*nf(i)-1;
                Edof(i,7)=3*nf(i);

                % Equivalent forces
                rxbar(i,2)=rxbar(i,2)+1.5*mplas/len(i)*abs(co(i));
                rxbar(i,3)=rxbar(i,3)-1.5*mplas/len(i)*abs(co(i));

                rybar(i,2)=rybar(i,2)-1.5*mplas/len(i)*abs(se(i));
                rybar(i,3)=rybar(i,3)+1.5*mplas/len(i)*abs(se(i));

                mbar(i,2)=mbar(i,2)+0.5*mplas;
                mbar(i,3)=mbar(i,3)+mplas;

                bar_plas_node=[bar_plas_node;
                                i,1,incfac];
                lambda_collection=[lambda_collection,incfac];

                iter_collection=[iter_collection,iteration];

                co=(coordxy(nf,1)-coordxy(ni,1))./len;       % cos direction vector
                se=(coordxy(nf,2)-coordxy(ni,2))./len;       % sen direction vector
            elseif abs(reac_collection(6,i))>=Mp(i,2) && ...
                    abs(reac_collection(3,i))<Mp(i,1)
                current_plas=1;
                plast_bars(2,i)=1;
                mplas=reac_collection(6,i);

                % change condition to Fixed-Art
                support(i,3)="Art";

                % Equivalent forces
                rxbar(i,2)=rxbar(i,2)+1.5*mplas/len(i)*abs(se(i));
                rxbar(i,3)=rxbar(i,3)-1.5*mplas/len(i)*abs(se(i));

                rybar(i,2)=rybar(i,2)-1.5*mplas/len(i)*abs(co(i));
                rybar(i,3)=rybar(i,3)+1.5*mplas/len(i)*abs(co(i));

                mbar(i,2)=mbar(i,2)+0.5*mplas;
                mbar(i,3)=mbar(i,3)+mplas;

                bar_plas_node=[bar_plas_node;
                                i,2,incfac];

                lambda_collection=[lambda_collection,incfac];

                iter_collection=[iter_collection,iteration];

                co=(coordxy(nf,1)-coordxy(ni,1))./len;       % cos direction vector
                se=(coordxy(nf,2)-coordxy(ni,2))./len;       % sen direction vector
            end
        end
    end

    if current_plas==0
        incfac=incfac+0.01;
    else
        
        if sum(seismic_forces)<0
            disp_iter=[];
            force_iter=[];
            for i=1:nfloors
                if i==1
                    disp_iter=[disp_iter;
                        abs(displacements(dof_disp(i)))];
                else
                    disp_iter=[disp_iter;
                        abs(displacements(dof_disp(i)))-abs(displacements(dof_disp(i-1)))];
                    
                end
                
                force_iter=[force_iter;
                            abs(fglobal(dof_forces(i)))]; 
                        
                figure(1)
                plot(disp_iter(i),...
                    force_iter(i),'b o','MarkerFaceColor','blue')
                xlabel('Lateral displacement (cm) - To the left')
                ylabel('Load (Kg)')
                title('Load-Displacement Historial per floor (Left)')
                hold on
                        
            end
            displacement_history_floor_left=[displacement_history_floor_left,disp_iter];                
            force_floor_left=[force_floor_left,force_iter];
            
        else
            disp_iter=[];
            force_iter=[];
            for i=1:nfloors
                if i==1
                    disp_iter=[disp_iter;
                        abs(displacements(dof_disp(i)))];
                else
                    disp_iter=[disp_iter;
                        abs(displacements(dof_disp(i)))-abs(displacements(dof_disp(i-1)))];
                        
                end
                force_iter=[force_iter;
                            abs(fglobal(dof_forces(i)))]; 
              
                figure(3)
                plot(disp_iter(i),...
                    force_iter(i),'k o','MarkerFaceColor','black')
                xlabel('Lateral displacement (cm) - To the right')
                ylabel('Load (Kg)')
                title('Load-Displacement Historial per floor (Right)')
                hold on
                        
            end
            displacement_history_floor_right=[displacement_history_floor_right,disp_iter];                
            force_floor_right=[force_floor_right,force_iter];
            
        end
        incfac=incfac+0.01;
    end
    
    if iteration>1
        if det_Kred_post/det_Kred<3e-2
            break;

        else

            continue;
        end
    end

end

maxDisplacement=zeros(nfloors,1);
k_direction=zeros(1,nfloors);
pdrift_DI=zeros(1,nfloors);
def_based_DI=zeros(1,nfloors);

drift_DI=zeros(1,nfloors);

if sum(seismic_forces)<0
    nd=length(displacement_history_floor_left(1,:));
    for i=1:nfloors
        k_direction(i)=force_floor_left(i,2)/displacement_history_floor_left(i,2);
        
        % Plastic drift damage
        % index_____________________________________________________________
        
        pdrift_DI(i)=(max(displacement_history_floor_left(i,:))-displacement_history_floor_left(i,2))/...
        (h_floor(i))*100;
    
        % Interstory drift damage index_____________________________

        drift_DI(i)=max(displacement_history_floor_left(i,:))/(h_floor(i))*100;
        
        maxDisplacement(i)=max(displacement_history_floor_left(i,:));
        
        
        delta_u=displacement_history_floor_left(i,nd);
        def_based_DI(i)=(maxDisplacement(i)-displacement_history_floor_left(i,2))/(delta_u-displacement_history_floor_left(i,2));
        
        figure(1)
        plot(displacement_history_floor_left(i,:),force_floor_left(i,:),'b')
        hold on
    end
        
    Ed=extract(Edof,Uglobal);

    %-----Undeformed mesh-----%
    figure(2)
    xlabel('Width [cm]')
    ylabel('Height [cm]')
    title('Deformed-undeformed structure - Forces to the left (Scale x 50)');
    plotpar=[2 1 0];
    eldraw2(Ex,Ey,plotpar);
    %----Deformed mesh---------%
    plotpar=[1 2 1];
    eldisp2(Ex,Ey,Ed,plotpar,50);
    
else
    nd=length(displacement_history_floor_right(1,:));
    for i=1:nfloors
        
        k_direction(i)=force_floor_right(i,2)/displacement_history_floor_right(i,2);
        
        % Plastic drift damage
        % index______________________________________________________________

        pdrift_DI(i)=(max(displacement_history_floor_right(i,:))-displacement_history_floor_right(i,2))/...
            (h_floor(i))*100;
        
        % Interstory drift damage index_____________________________
        drift_DI(i)=max(displacement_history_floor_right(i,:))/(h_floor(i))*100;
        
        maxDisplacement(i)=max(displacement_history_floor_right(i,:));
    
        delta_u=displacement_history_floor_right(i,nd);
        
        def_based_DI(i)=(maxDisplacement(i)-displacement_history_floor_right(i,2))/(delta_u-displacement_history_floor_right(i,2));
        
        figure(3)
        plot(displacement_history_floor_right(i,:),force_floor_right(i,:),'k')
        hold on
    end

    Ed=extract(Edof,Uglobal);

    %-----Undeformed mesh-----%
    figure(4)
    xlabel('Width [cm]')
    ylabel('Height [cm]')
    title('Deformed-undeformed structure - Forces to the right (Scale x 50)');
    plotpar=[2 1 0];
    eldraw2(Ex,Ey,plotpar);
    %----Deformed mesh---------%
    plotpar=[1 2 1];
    eldisp2(Ex,Ey,Ed,plotpar,50);
end

%---------------------------------end----------------------------------