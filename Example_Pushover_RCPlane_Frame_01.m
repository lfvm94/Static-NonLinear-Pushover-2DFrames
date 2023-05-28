clc
clear all

% Example_Pushover_RCPlane_Frame_01
%----------------------------------------------------------------
% PURPOSE 
%    To compute the non-linear static Pushover analysis for a
%    reinforced concrete plane frame
%
%
%----------------------------------------------------------------

% LAST MODIFIED: L.F.Veduzco    2023-02-23
%                Faculty of Engineering
%                Autonomous University of Queretaro
%----------------------------------------------------------------
clc
clear all

nnodes=8;
nbars=8;

%% Materials
% f'c of each element
fpc=[300;
    300;
    250;
    300;
    300;
    250;
    250;
    300];

% Elasticity modulus of each element in function of f'c
e=zeros(nbars,1);
for i=1:nbars
    e(i)=14000*sqrt(fpc(i));
end

%% Geometry/Topology
dimensions=[40 40;40 40;30 60;40 40;50 50;30 60;30 50;30 30];
        
% cross-section area of each element
a=zeros(nbars,1);
for i=1:nbars
    a(i)=dimensions(i,1)*dimensions(i,2);
end

inertia=zeros(nbars,1);
for i=1:nbars
    inertia(i)=1/12*dimensions(i,1)*dimensions(i,2)^3;
end


%coordinates of each node for each bar
coordxy=[0 -100;0 400;0 800;600 800;600 400;600 -100;1200 400;1200 -100]; 

% final and initial nodes for each element
ni=[1;2;3;4;5;2;5;7];
nf=[2;3;4;5;6;5;7;8];

l=sqrt((coordxy(nf,1)-coordxy(ni,1)).^2+...
      (coordxy(nf,2)-coordxy(ni,2)).^2); % bar-length vector

c=(coordxy(nf,1)-coordxy(ni,1))./l;      % cos direction vector
s=(coordxy(nf,2)-coordxy(ni,2))./l;      % sen direction vector

% prescribed boudnary conditions
bc=[1 0;2 0;3 0;16 0;17 0;18 0;22 0;23 0;24 0];

% habilited nodes
[ndof,edof]=nonRestrcDof(nnodes,bc);

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


supports=[1 "Fixed" "Fixed";
           2 "Fixed" "Fixed";
           3 "Fixed" "Fixed";
           4 "Fixed" "Fixed";
           5 "Fixed" "Fixed";
           6 "Fixed" "Fixed";
           7 "Fixed" "Fixed";
           8 "Fixed" "Fixed"];
   
type_elem=[1 "Col";
           2 "Col";
           3 "Beam";
           4 "Col";
           5 "Col";
           6 "Beam";
           7 "Beam";
           8 "Col"];
       
%% Loads       
beams_LL=[1 100; % Uniformly distributed loads over the beams
          2 100;
          3 100];
      
% To consider the self-weight as uniformly distributed load on the 
% elements
unitWeightElem=zeros(nbars,2);
for i=1:nbars
    unitWeightElem(i,2)=0.0024; % kg/cm3
end


elemcols=[];
elembeams=[];
beams=0;
cols=0;
for j=1:nbars
    if type_elem(j,2)=="Beam"
        beams=beams+1;
        elembeams=[elembeams,j];
    elseif type_elem(j,2)=="Col"
        cols=cols+1;
        elemcols=[elemcols,j];
    end
end

% Uniformly distributed loads considering self weight of the elements
qbary=zeros(nbars,2);
for i=1:beams
    qbary(elembeams(i),2)=1.1*a(elembeams(i))*unitWeightElem(elembeams(i),2)+...
                           1.1*(beams_LL(i,2));
    
end

% Distribution of loads to the end of each element (Do not delete)
rxbar=zeros(nbars,3);
rybar=zeros(nbars,3);
mbar=zeros(nbars,3);
for i=1:nbars
    rybar(i,2)=qbary(i,2)*l(i)*0.5;
    rybar(i,3)=qbary(i,2)*l(i)*0.5;
    
    mbar(i,2)=qbary(i,2)*l(i)^2/12;
    mbar(i,3)=-qbary(i,2)*l(i)^2/12;
end
   
% Plastic moments of each element's ends
Mp=[7680000 7680000;
    6490000 6490000;
    8363000 8976940;
    5490000 5490000;
    8680000 8680000;
    9363000 9976940;
    7363000 7976940;
    5490000 5490000]; %Kg-cm

% Lateral equivalent seismic forces from a modal analysis. The number of
% forces must be equal to the number of floors
seismic_forces=[2000; % upper floor
                1500]; % lower floor
            
% Degrees of freedom over which each seismic force is applied (one for
% each seismic force)
dof_seismic_forces=[4 7];

% Degrees of freedom corresponding to those over which the seismic forces
% are applied, in the same order
dof_disp=[1 4];

% Height of each floor
hfloor=[400; 400];  

%%% IN POSITIVE DIRECTION OF FORCES__________
%____________________________________________________________________


[lambda_der,pdrift_DI_der,drift_DI_der,def_based_di_der,...
max_disp_der]=ElastoPlasticPushoverPlaneFrames(qbary,a,Mp,nbars,...
nnodes,e,inertia,coordxy,ni,nf,l,c,s,Edof,supports,edof,ndof,rxbar,...
rybar,mbar,seismic_forces,hfloor,dof_seismic_forces,dof_disp);

%%% IN NEGATIVE DIRECTION OF FORCES__________
%____________________________________________________________________

seismic_forces=-seismic_forces;
    
[lambda_izq,pdrift_DI_izq,drift_DI_izq,def_based_di_izq,...
max_disp_izq]=ElastoPlasticPushoverPlaneFrames(qbary,a,Mp,nbars,...
nnodes,e,inertia,coordxy,ni,nf,l,c,s,Edof,supports,edof,ndof,rxbar,...
rybar,mbar,seismic_forces,hfloor,dof_seismic_forces,dof_disp);

nfloors=length(hfloor);

% Final results
FS=min([lambda_der, lambda_izq])
pdriftDI=min([sum(pdrift_DI_der)/nfloors,sum(pdrift_DI_izq)/nfloors])
driftDI=min([sum(drift_DI_der)/nfloors,sum(drift_DI_izq)/nfloors])
dbDI=min([sum(def_based_di_der)/nfloors,sum(def_based_di_izq)/nfloors])

Max_Displacement=max(max(max_disp_izq),max(max_disp_der))