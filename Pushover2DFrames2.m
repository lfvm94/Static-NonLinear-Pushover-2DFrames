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

dispHistFloorLeft=[];
force_floor_left=[];

dispHistFloorRight=[];
force_floor_right=[];

zero_disp=zeros(nfloors,1);
dispHistFloorLeft=[dispHistFloorLeft,zero_disp];
force_floor_left=[force_floor_left,zero_disp];

dispHistFloorRight=[dispHistFloorRight,zero_disp];
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
    fglobal(dofForces)=seismicforces*incLoad;
    
    elmMat=zeros(6*nbars,6);

    Ex=zeros(nbars,2);
    Ey=zeros(nbars,2);
    
    for i=1:nbars 
        ex=[coordxy(ni(i),1) coordxy(nf(i),1)];
        ey=[coordxy(ni(i),2) coordxy(nf(i),2)];
        
        Ex(i,:)=ex;
        Ey(i,:)=ey;
        
        ep=[E(i) A(i) I(i)];
         
        eq=[0 qbary(i,2)];

        [Kebar,febar]=beam2e(ex,ey,ep,eq); % This is a CALFEM function
                                           % Download at: 
                                           % https://www.byggmek.lth.se/english/calfem/
        
        if support(i,2)=="Fixed" && support(i,3)=="Art"
            Mpl=mbar(i,2);
            eq=[0 qbary(i,2) Mpl];
            [Kebar,febar]=beamArt2e(ex,ey,ep,eq,1);
             
         elseif support(i,2)=="Art" && support(i,3)=="Fixed"

             Mpl=mbar(i,1);
             eq=[0 qbary(i,2) Mpl];
             [Kebar,febar]=beamArt2e(ex,ey,ep,eq,2);
             
        elseif support(i,2)=="Art" && support(i,3)=="Art"

             Mp1=mbar(i,1);
             Mp2=mbar(i,2);
             eq=[0 qbary(i,2) [Mp1,Mp2]];
             [Kebar,febar]=beamArt2e(ex,ey,ep,eq,3);
        end
        fbe(:,i)=febar; % storing elemental forces for further use
        
        elmMat((i-1)*6+1:6*i,:)=Kebar; % storing Ke of bars for further use
        
        % Assembling global stiffness matrix
        [Kglobal,fglobal]=assem(Edof(i,:),Kglobal,Kebar,fglobal,febar);
     end     

    globalKreduced=Kglobal(edof,edof);
    if iteration==0
        det_Kred=det(globalKreduced);
    else
        det_Kred_post=det(globalKreduced);
    end
    % Solving the system of equations
    [Uglobal,Reactions]=solveq(Kglobal,fglobal,bc);
   
    % --- computation of mechanic elements at the ends of bars --- %
    Ed=extract(Edof,Uglobal);
    for i=1:nbars
        
        es_bar=beam2s(Ex(i,:),Ey(i,:),[E(i) A(i) I(i)],Ed(i,:),...
            [0 qbary(i,2)],2);
        
        ue=Uglobal(Edof(i,2:7));
        ke=elmMat((i-1)*6+1:6*i,:);

        fe=-(ke*ue-fbe(:,i));
        
        reac_bars(:,i)=fe;
    end

    iteration=iteration+1;
    current_plas=0; % to register if there is a plastification in the 
                    % current run

    plastified_bars=zeros(nbars,1); % to register which bars are plastified
                                    % in the current run (if any)
                                    
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

                    % Change condition Art-Fixed to Art-Art
                    support(i,3)="Art";

                    % Equivalent plastic moments
                    mplas=reac_bars(6,i);
                    mbar(i,2)=mplas;
    
                    plastified_bars(i,1)=2;
                
                elseif plast_bars(1,i)==0 && plast_bars(2,i)==1
                    % The bar is currently Fixed-Art and will be Art-Art
                    current_plas=1;
                    plast_bars(1,i)=1;
                    
                    % Change condition Fixed-Art to Art-Art
                    support(i,2)="Art";

                    % Equivalent plastic moments
                    mplas=reac_bars(3,i);
                    mbar(i,1)=mplas;
                    
                    plastified_bars(i,1)=1;
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

                % Equivalent plastic moments
                mbar(i,1)=mplas;
                
                plastified_bars(i,1)=1;

            elseif abs(reac_bars(6,i))>=Mp(i,2) && ...
                    abs(reac_bars(3,i))<Mp(i,1)
                current_plas=1;
                plast_bars(2,i)=1;
                mplas=reac_bars(6,i);

                % change condition to Fixed-Art
                support(i,3)="Art";

                % Equivalent plastic moments
                mbar(i,2)=mplas;
                
                plastified_bars(i,1)=2;
                
            elseif abs(reac_bars(6,i))>=Mp(i,2) && ...
                    abs(reac_bars(3,i))>=Mp(i,1)
                
                current_plas=1;
                plast_bars(2,i)=1; 
                plast_bars(1,i)=1;
                mplas1=reac_bars(3,i); % storing plastic moments
                mplas2=reac_bars(6,i); % for the next iteration
                
                % change condition to Art-Art
                support(i,2)="Art";
                support(i,3)="Art";

                % Equivalent plastic moments
                mbar(i,1)=mplas1;
                mbar(i,2)=mplas2;
                                
                % To have register that both element's ends were
                % plastified in flexure
                plastified_bars(i,1)=3;

            end
        end
    end
    if current_plas==0
        % Updating loads for the next iteration (in case there is one)
        incLoad=incLoad+dload;
    else
        historyIncLoad=[historyIncLoad,incLoad];
        iter_collection=[iter_collection,iteration];
                
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
                               abs(Uglobal(dofForces(i)))-...
                               abs(Uglobal(dofForces(i-1)))];
                    
                end
                
                force_iter=[force_iter;
                            abs(seismicforces(i))*incLoad]; 
                     
            end
            
            dispHistFloorLeft=[dispHistFloorLeft,disp_iter];                
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
                        abs(Uglobal(dofForces(i)))-...
                        abs(Uglobal(dofForces(i-1)))];
                        
                end
                force_iter=[force_iter;
                            abs(seismicforces(i))*incLoad]; 
                
            end
            dispHistFloorRight=[dispHistFloorRight,disp_iter];                
            force_floor_right=[force_floor_right,force_iter];
            
        end
        % Updating loads for the next iteration (in case there is one)
        incLoad=incLoad+dload;
    end
    
    if iteration>1 % to verify is a next run is required
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
    nd=length(dispHistFloorLeft(1,:));
    for i=1:nfloors
        k_direction(i)=force_floor_left(i,2)/dispHistFloorLeft(i,2);
        
        % Plastic drift damage index
        
        pdriftDI(i)=(max(dispHistFloorLeft(i,:))-dispHistFloorLeft(i,2))/...
        (Hfloor(i))*100;
    
        % Interstory drift damage index

        driftDI(i)=max(dispHistFloorLeft(i,:))/(Hfloor(i))*100;
        
        maxDisplacement(i)=max(dispHistFloorLeft(i,:));
        
        % Deformation-Based Damage Index
        du=0.04*Hfloor(i);
        defBasedDI(i)=(maxDisplacement(i)-dispHistFloorLeft(i,2))/...
                      (du-dispHistFloorLeft(i,2));
        
        floorText(i,:)=strcat('Floor ',num2str(i));
        figure(1)
        if i==1
            plot(dispHistFloorLeft(i,:),...
                force_floor_left(i,:),'b -','LineWidth',1.8)
            legend(floorText(i,:))
            hold on
        else
            plot(dispHistFloorLeft(i,:),...
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
    nd=length(dispHistFloorRight(1,:));
    for i=1:nfloors
        
        k_direction(i)=force_floor_right(i,2)/dispHistFloorRight(i,2);
        
        % Plastic drift damage index
        pdriftDI(i)=(max(dispHistFloorRight(i,:))-dispHistFloorRight(i,2))/...
            (Hfloor(i))*100;
        
        % Interstory drift damage index
        driftDI(i)=max(dispHistFloorRight(i,:))/(Hfloor(i))*100;
        
        maxDisplacement(i)=max(dispHistFloorRight(i,:));
        
        % Deformation-Based Damage Index
        du=0.04*Hfloor(i);
        defBasedDI(i)=(maxDisplacement(i)-dispHistFloorRight(i,2))/...
                       (du-dispHistFloorRight(i,2));
        
        floorText(i,:)=strcat('Floor ',num2str(i));
        figure(3)
        if i==1
            plot(dispHistFloorRight(i,:),...
                force_floor_right(i,:),'k -','LineWidth',1.8)
            legend(floorText(i,:))
            hold on
        else
            plot(dispHistFloorRight(i,:),...
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