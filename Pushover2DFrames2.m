function [historyIncLoad,pdriftDI,driftDI,defBasedDI,maxDisplacement,...
         barPlasNode]=Pushover2DFrames2(qbary,A,Mp,E,I,coordxy,ni,nf,...
         support,bc,seismicforces,Hfloor,dofForces,dload,kdam)

%------------------------------------------------------------------------
% [historyIncLoad,pdriftDI,driftDI,defBasedDI,maxDisplacement,...
%  barPlasNode]=Pushover2DFrames2(qbary,A,Mp,E,I,coordxy,ni,nf,support,...
%  bc,seismicforces,Hfloor,dofForces,dload,kdam)
%
%------------------------------------------------------------------------
% PURPOSE
%  To compute a static non-linear pushover analysis of a plane frame
%  
% 
% INPUT:  A = [area_bar;
%               ...]                 area of all elements
%
%         Mp = [Mpi Mpj;             Plastic Moment for each member 
%               ... ]                (i) initial node, (j) final node
%
%         E = [e_bar;                Elasticity modulus of each element
%               ...]                    
%
%         I = [inertia_bar;         in-plane inertia for all elements'
%                       ...]         cross-section
%                                    
%         coordxy = [coordx coordy;  node coordinates for all nodes
%                       ...];
%
%         ni                         list of initial nodes of all bars,
%         nf                         list of final nodes of all bars:
%                                         size = [nbars,1]
% 
%         qbary = [bar, load;
%                   ..    .. ]       uniform distributed load acting 
%                                    downwards (only for beams)
%
%         support = [i, j]           support at each bar's end
%                                    options: "Art" or "Fixed"
%                                    (i) initial node, (j) final node
%
%         bc                         restricted dof
%
%         seismicForces = [f(1);]    lateral forces per floor:
%                          f(n);]    size = [nfloors,1]
%
%         Hfloor = [h(1);            Height of each floor from bottom
%                    h(n)]           to top: size = [nfloors,1]
%
%         dofForces = [dof-f(1),     dof at which the lateral forces are
%                       dof-f(n)]    applied (from bottom to top) - global
%
% OUTPUT: historyIncLoad             history of incremental load factors at
%                                    at which plastic moments are reached
%
%         pdriftDI                   Plastic inter-story drift Damage 
%                                    Index per floor: size = [nfloors,1]
%
%         driftDI                    Inter-story drift Damage 
%                                    Index per floor: size = [nfloors,1]
%
%         defBasedDI                 Deformation based Damage Index
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

nbars=length(E);
nnodes=length(coordxy(:,1));

% Topology matrix 
Edof=zeros(nbars,7);
for i=1:nbars
    Edof(i,1)=i;
    Edof(i,2)=ni(i)*3-2;
    Edof(i,3)=ni(i)*3-1;
    Edof(i,4)=ni(i)*3;
    
    Edof(i,5)=nf(i)*3-2;
    Edof(i,6)=nf(i)*3-1;
    Edof(i,7)=nf(i)*3;
    
end

[ndof,edof]=nonRestrcDof(nnodes,bc);

mbar=zeros(nbars,2); % to save the plastic moments at each articulation
                     % of each bar as a plastification occurs
nfloors=length(Hfloor); % will be used to compute the relative floor
                         % displacements and damage indices of each floor

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
barPlasNode=[];
iter_collection=[];
historyIncLoad=[];
looping=0;
incLoad=1.0;
iteration=0;
while looping==0
    Kglobal=zeros(3*nnodes);
    fglobal=zeros(3*nnodes,1);

    for i=1:nfloors
        fglobal(dofForces(i))=seismicforces(i)*incLoad;
    end

    elemental_matrices=zeros(6*nbars,6);

    for i=1:nbars      
        
        ex=[coordxy(ni(i),1) coordxy(nf(i),1)];
        ey=[coordxy(ni(i),2) coordxy(nf(i),2)];
        ep=[E(i) A(i) I(i)];
         
        eq=[0 -qbary(i,2)];

        [Ke_barra,fe_barra]=beam2e(ex,ey,ep,eq); % This is a CALFEM
                                                 % function
                                                 % Download at: 
                                                 % https://www.byggmek.lth.se/english/calfem/

        if support(i,2)=="Fixed" && support(i,3)=="Art"
            Mpl=mbar(i,2);
            eq=[0 -qbary(i,2) Mpl];
            [Ke_barra,fe_barra]=beamArt2e(ex,ey,ep,eq,1);
             
         elseif support(i,2)=="Art" && support(i,3)=="Fixed"

             Mpl=mbar(i,1);
             eq=[0 -qbary(i,2) Mpl];
             [Ke_barra,fe_barra]=beamArt2e(ex,ey,ep,eq,2);
             
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
    
    [Uglobal,Reactions]=solveq(Kglobal,fglobal,bc);
    
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

    reac_bars=zeros(6,nbars);
    
    % --- computation of mechanic elements at the ends of bars --- %
    for i=1:nbars
        ue=Uglobal(Edof(i,2:7));
        ke=elemental_matrices((i-1)*6+1:6*i,:);

        fe=-ke*ue;

        reac_bars(:,i)=fe;
    end

    iteration=iteration+1;
    current_plas=0;

    plastified_bars=zeros(nbars,1);
    for i=1:nbars

        % Detect if any end of this bar (i) has been plastified
        bar_plas_check=0;
        if plast_bars(1,i)~=0 || plast_bars(2,i)~=0
            bar_plas_check=bar_plas_check+1;
        end

        if bar_plas_check==1
            % Detect if the other end has been also plastified
            if abs(reac_bars(3,i))>=Mp(i,1) && ...
               abs(reac_bars(6,i))>=Mp(i,2) % if both ends are plastified
           
                if plast_bars(1,i)==1 && plast_bars(2,i)==0 
                    % The bar is currently Art-Fixed and will be Art-Art
                    current_plas=1;

                    plast_bars(2,i)=1;

                    I(i)=1e-10;

                    % Change condition Fixed-Art to Art-Art
                    support(i,3)="Art";

                    mplas=reac_bars(6,i);
                    mbar(i,2)=mbar(i,2)+mplas;

                    % These are substituiting forces for each plastified bar 
                    fglobal(ni(i)*3-2)=fglobal(ni(i)*3-2)+...
                                            reac_bars(1,i);
                    fglobal(nf(i)*3-2)=fglobal(nf(i)*3-2)+...
                                            reac_bars(4,i);

                    fglobal(ni(i)*3-1)=fglobal(ni(i)*3-1)+...
                                            reac_bars(2,i);
                    fglobal(nf(i)*3-1)=fglobal(nf(i)*3-1)+...
                                            reac_bars(5,i);

                    fglobal(ni(i)*3)=fglobal(ni(i)*3)+...
                                        reac_bars(3,i);
                    fglobal(nf(i)*3)=fglobal(nf(i)*3)+...
                                        reac_bars(6,i);
                    
                    plastified_bars(i,1)=2;
                    
                    historyIncLoad=[historyIncLoad,incLoad];

                    iter_collection=[iter_collection,iteration];
                
                elseif plast_bars(1,i)==0 && plast_bars(2,i)==1
                    % The bar is currently Fixed-Art and will be Art-Art
                    current_plas=1;
                    plast_bars(1,i)=1;

                    I(i)=1e-10;

                    % Change condition Fixed-Art to Art-Art
                    support(i,2)="Art";

                    mplas=reac_bars(3,i);
                    mbar(i,1)=mbar(i,1)+mplas;

                    % These are substituiting forces for each plastified bar 
                    fglobal(ni(i)*3-2)=fglobal(ni(i)*3-2)+...
                                            reac_bars(1,i);
                    fglobal(nf(i)*3-2)=fglobal(nf(i)*3-2)+...
                                            reac_bars(4,i);

                    fglobal(ni(i)*3-1)=fglobal(ni(i)*3-1)+...
                                            reac_bars(2,i);
                    fglobal(nf(i)*3-1)=fglobal(nf(i)*3-1)+...
                                            reac_bars(5,i);

                    fglobal(ni(i)*3)=fglobal(ni(i)*3)+...
                                        reac_bars(3,i);
                    fglobal(nf(i)*3)=fglobal(nf(i)*3)+...
                                        reac_bars(6,i);
                    plastified_bars(i,1)=1;
                                
                    historyIncLoad=[historyIncLoad,incLoad];

                    iter_collection=[iter_collection,iteration];
                end
            end
        elseif bar_plas_check==0

            if abs(reac_bars(3,i))>=Mp(i,1) && ...
                    abs(reac_bars(6,i))<Mp(i,2)

                current_plas=1;
                mplas=reac_bars(3,i);
                plast_bars(1,i)=1;

                % change condition to Art
                support(i,2)="Art";

                % Equivalent forces
                mbar(i,1)=mbar(i,1)+mplas;
                mbar(i,2)=mbar(i,2)+0.5*mplas;

                % These are substituiting forces for each plastified bar 
                fglobal(ni(i)*3-2)=fglobal(ni(i)*3-2)+...
                                        reac_bars(1,i);
                fglobal(nf(i)*3-2)=fglobal(nf(i)*3-2)+...
                                        reac_bars(4,i);

                fglobal(ni(i)*3-1)=fglobal(ni(i)*3-1)+...
                                        reac_bars(2,i);
                fglobal(nf(i)*3-1)=fglobal(nf(i)*3-1)+...
                                        reac_bars(5,i);

                fglobal(ni(i)*3)=fglobal(ni(i)*3)+...
                                    reac_bars(3,i);
                fglobal(nf(i)*3)=fglobal(nf(i)*3)+...
                                    reac_bars(6,i);
                plastified_bars(i,1)=1;
                
                historyIncLoad=[historyIncLoad,incLoad];

                iter_collection=[iter_collection,iteration];

            elseif abs(reac_bars(6,i))>=Mp(i,2) && ...
                    abs(reac_bars(3,i))<Mp(i,1)
                current_plas=1;
                plast_bars(2,i)=1;
                mplas=reac_bars(6,i);

                % change condition to Fixed-Art
                support(i,3)="Art";

                % Equivalent forces
                
                mbar(i,1)=mbar(i,1)+0.5*mplas;
                mbar(i,2)=mbar(i,2)+mplas;

                % These are substituiting forces for each plastified bar 
                fglobal(ni(i)*3-2)=fglobal(ni(i)*3-2)+...
                                        reac_bars(1,i);
                fglobal(nf(i)*3-2)=fglobal(nf(i)*3-2)+...
                                        reac_bars(4,i);

                fglobal(ni(i)*3-1)=fglobal(ni(i)*3-1)+...
                                        reac_bars(2,i);
                fglobal(nf(i)*3-1)=fglobal(nf(i)*3-1)+...
                                        reac_bars(5,i);

                fglobal(ni(i)*3)=fglobal(ni(i)*3)+...
                                    reac_bars(3,i);
                fglobal(nf(i)*3)=fglobal(nf(i)*3)+...
                                    reac_bars(6,i);
                plastified_bars(i,1)=2;

                historyIncLoad=[historyIncLoad,incLoad];

                iter_collection=[iter_collection,iteration];

            end
        end
    end
                
    if current_plas==0
        incLoad=incLoad+dload;
    else
        barPlasNode=[barPlasNode,plastified_bars];
        if sum(seismicforces)<0
            disp_iter=[];
            force_iter=[];
            for i=1:nfloors
                if i==1
                    disp_iter=[disp_iter;
                        abs(Uglobal(dofForces(i)))];
                else
                    % Relative displacement
                    disp_iter=[disp_iter;
                        abs(Uglobal(dofForces(i)))-abs(Uglobal(dofForces(i-1)))];
                    
                end
                
                force_iter=[force_iter;
                            abs(seismicforces(i))*incLoad]; 
                     
            end
            
            displacement_history_floor_left=[displacement_history_floor_left,disp_iter];                
            force_floor_left=[force_floor_left,force_iter];
            
        else
            disp_iter=[];
            force_iter=[];
            for i=1:nfloors
                if i==1
                    disp_iter=[disp_iter;
                        abs(Uglobal(dofForces(i)))];
                else
                    disp_iter=[disp_iter;
                        abs(Uglobal(dofForces(i)))-abs(Uglobal(dofForces(i-1)))];
                        
                end
                force_iter=[force_iter;
                            abs(seismicforces(i))*incLoad]; 
                
            end
            displacement_history_floor_right=[displacement_history_floor_right,disp_iter];                
            force_floor_right=[force_floor_right,force_iter];
            
        end
        incLoad=incLoad+dload;
    end
    
    if iteration>1
        if det_Kred_post/det_Kred<kdam
            break;

        else

            continue;
        end
    end

end
            
maxDisplacement=zeros(nfloors,1);
k_direction=zeros(1,nfloors);
pdriftDI=zeros(1,nfloors);
defBasedDI=zeros(1,nfloors);

driftDI=zeros(1,nfloors);

if sum(seismicforces)<0
    nd=length(displacement_history_floor_left(1,:));
    for i=1:nfloors
        k_direction(i)=force_floor_left(i,2)/displacement_history_floor_left(i,2);
        
        % Plastic drift damage index
        
        pdriftDI(i)=(max(displacement_history_floor_left(i,:))-displacement_history_floor_left(i,2))/...
        (Hfloor(i))*100;
    
        % Interstory drift damage index

        driftDI(i)=max(displacement_history_floor_left(i,:))/(Hfloor(i))*100;
        
        maxDisplacement(i)=max(displacement_history_floor_left(i,:));
        
        delta_u=displacement_history_floor_left(i,nd);
        defBasedDI(i)=(maxDisplacement(i)-displacement_history_floor_left(i,2))/(delta_u-displacement_history_floor_left(i,2));
        
        floorText(i,:)=strcat('Floor ',num2str(i));
        figure(1)
        if i==1
            plot(displacement_history_floor_left(i,:),...
                force_floor_left(i,:),'b -','LineWidth',1.8)
            legend(floorText(i,:))
            hold on
        else
            plot(displacement_history_floor_left(i,:),...
                force_floor_left(i,:),'-',...
                'LineWidth',1.8,'DisplayName',floorText(i,:))
            hold on
        end   
    end
    xlabel('Lateral relative displacement')
    ylabel('Load')
    title('Load-Displacement Historial per floor (Left)')
    hold on
        
    Ed=extract(Edof,Uglobal);

    %-----Undeformed mesh-----%
    figure(2)
    xlabel('Width')
    ylabel('Height')
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
        
        % Plastic drift damage index
        pdriftDI(i)=(max(displacement_history_floor_right(i,:))-displacement_history_floor_right(i,2))/...
            (Hfloor(i))*100;
        
        % Interstory drift damage index
        driftDI(i)=max(displacement_history_floor_right(i,:))/(Hfloor(i))*100;
        
        maxDisplacement(i)=max(displacement_history_floor_right(i,:));
    
        delta_u=displacement_history_floor_right(i,nd);
        
        defBasedDI(i)=(maxDisplacement(i)-displacement_history_floor_right(i,2))/(delta_u-displacement_history_floor_right(i,2));
        
        floorText(i,:)=strcat('Floor ',num2str(i));
        figure(3)
        if i==1
            plot(displacement_history_floor_right(i,:),...
                force_floor_right(i,:),'k -','LineWidth',1.8)
            legend(floorText(i,:))
            hold on
        else
            plot(displacement_history_floor_right(i,:),...
                force_floor_right(i,:),'-','LineWidth',1.8,...
                'DisplayName',floorText(i,:))
            hold on
        end
    end
    xlabel('Lateral relative displacement')
    ylabel('Load')
    title('Load-Displacement Historial per floor (Right)')
    hold on
    
    Ed=extract(Edof,Uglobal);

    %-----Undeformed mesh-----%
    figure(4)
    xlabel('Width')
    ylabel('Height')
    title('Deformed-undeformed structure - Forces to the right (Scale x 50)');
    plotpar=[2 1 0];
    eldraw2(Ex,Ey,plotpar);
    %----Deformed mesh---------%
    plotpar=[1 2 1];
    eldisp2(Ex,Ey,Ed,plotpar,50);
end

%---------------------------------end----------------------------------