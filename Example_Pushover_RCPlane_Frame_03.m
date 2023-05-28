clc 
clear all

% Example_Pushover_RCPlane_Frame_03
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


nnodes=15;
nbars=20;

%% Materials
% f'c of each element
fpc=[300;
     300;
     300;
     300;
     280;
     280;
     280;
     280;
     300;
     300;
     300;
     300;
     280;
     280;
     280;
     280;
     300;
     300;
     300;
     300];

% Elasticity modulus of each element in function of f'c
e=zeros(nbars,1);
for i=1:nbars
    e(i)=14000*(fpc(i))^0.5;
end

%% Geometry/Topology
dimensions=[35 70;
            65 70;
            55 70;
            55 55;
            55 110;
            55 110;
            35 70;
            35 70;
            65 105;
            65 90;
            65 90;
            65 205;
            35 70;
            35 70;
            35 70;
            35 70;
            65 70;
            55 70;
            45 60;
            35 60];


a=zeros(nbars,1); % cross-section area of each element
inertia=zeros(nbars,1);
for i=1:nbars
    a(i)=dimensions(i,1)*dimensions(i,2);
    inertia(i)=1/12*dimensions(i,1)*dimensions(i,2)^3;
end

%coordinates of each node for each bar
coordxy=[0 -150;
         0 400;
         0 800;
         0 1200;
         0 1600;
         600 1600;
         600 1200;
         600 800;
         600 400;
         600 -150;
         1200 -150;
         1200 400;
         1200 800;
         1200 1200;
         1200 1600]; 
     
                  
%%%---- Initial-final node of each bar -----%%%

ni=[1;2;3;4;5;4;3;2;10;9;8;7;6; 7; 8; 9; 11;12;13;14];
nf=[2;3;4;5;6;7;8;9;9; 8;7;6;15;14;13;12;12;13;14;15];

l=sqrt((coordxy(nf,1)-coordxy(ni,1)).^2+...
      (coordxy(nf,2)-coordxy(ni,2)).^2); % bar-length vector

c=(coordxy(nf,1)-coordxy(ni,1))./l;       % cos direction vector
s=(coordxy(nf,2)-coordxy(ni,2))./l;       % sen direction vector

% prescribed boudnary conditions
bc=[1 0;
    2 0;
    3 0;
    28 0;
    29 0;
    30 0;
    31 0;
    32 0;
    33 0];

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

%% Loads  
type_elem=[1 "Col";
           2 "Col";
           3 "Col";
           4 "Col";
           5 "Beam";
           6 "Beam";
           7 "Beam";
           8 "Beam";
           9 "Col";
           10 "Col";
           11 "Col";
           12 "Col";
           13 "Beam";
           14 "Beam";
           15 "Beam";
           16 "Beam";
           17 "Col";
           18 "Col";
           19 "Col";
           20 "Col"];
       
beams_LL=[1 100; % Uniformly distributed loads over the beams
          2 100;
          3 100;
          4 100;
          5 100;
          6 100;
          7 100;
          8 100];

elem_cols=[];
elem_beams=[];
beams=0;
cols=0;
for j=1:nbars
    if type_elem(j,2)=="Beam"
        beams=beams+1;
        elem_beams=[elem_beams,j];
    elseif type_elem(j,2)=="Col"
        cols=cols+1;
        elem_cols=[elem_cols,j];
    end
end


supports=[1 "Empotrado" "Empotrado";
       2 "Empotrado" "Empotrado";
       3 "Empotrado" "Empotrado";
       4 "Empotrado" "Empotrado";
       5 "Empotrado" "Empotrado";
       6 "Empotrado" "Empotrado";
       7 "Empotrado" "Empotrado";
       8 "Empotrado" "Empotrado";
       9 "Empotrado" "Empotrado";
       10 "Empotrado" "Empotrado";
       11 "Empotrado" "Empotrado";
       12 "Empotrado" "Empotrado";
       13 "Empotrado" "Empotrado";
       14 "Empotrado" "Empotrado";
       15 "Empotrado" "Empotrado";
       16 "Empotrado" "Empotrado";
       17 "Empotrado" "Empotrado";
       18 "Empotrado" "Empotrado";
       19 "Empotrado" "Empotrado";
       20 "Empotrado" "Empotrado"];
         
% To consider the self-weight as uniformly distributed load on the 
% elements
unit_weight_elem=zeros(nbars,2);
for i=1:nbars
    unit_weight_elem(i,2)=0.0024; % kg/cm3
end

% Uniformly distributed loads considering self weight of the elements
qbary=zeros(nbars,2);
for i=1:beams
    qbary(elem_beams(i),2)=1.1*a(elem_beams(i))*unit_weight_elem(elem_beams(i),2)+1.1*(beams_LL(i,2));
    
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
Mp=[9680000 9680000;
    8490000 8490000;
    8363000 8976940;
    7490000 7490000;
    5680000 5680000;
    7363000 7976940;
    8363000 8976940;
    9490000 9490000;
    12680000 12680000;
    11490000 11490000;
    10363000 10976940;
    9490000 9490000;
    5680000 5680000;
    7363000 7976940;
    8363000 8976940;
    9490000 9490000;
    9680000 9680000;
    8490000 8490000;
    8363000 8976940;
    7490000 7490000]; %Kg-cm

% Lateral equivalent seismic forces from a modal analysis. The number of
% forces must be equal to the number of floors
seismic_forces=[1500; % upper floor
                2000;
                2500;
                3000]; % lower floor

% Degrees of freedom over which each seismic force is applied (one for
% each seismic force)
dof_seismic_forces=[4 7 10 13];

% Degrees of freedom corresponding to those over which the seismic forces
% are applied, in the same order
dof_disp=[1 4 7 10];

% Height of each floor
hfloor=[400; 400; 400; 400];

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
FS=min([lambda_der, lambda_izq])
pdriftDI=min([sum(pdrift_DI_der)/nfloors,sum(pdrift_DI_izq)/nfloors])
driftDI=min([sum(drift_DI_der)/nfloors,sum(drift_DI_izq)/nfloors])
dbDI=min([sum(def_based_di_der)/nfloors,sum(def_based_di_izq)/nfloors])

Max_Displacement=max(max(max_disp_izq),max(max_disp_der))
