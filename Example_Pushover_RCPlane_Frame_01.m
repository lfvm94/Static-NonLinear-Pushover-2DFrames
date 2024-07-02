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

% LAST MODIFIED: L.F.Veduzco    2023-05-31
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
E=zeros(nbars,1);
for i=1:nbars
    E(i)=14000*sqrt(fpc(i));
end

%% Geometry/Topology
% cross-section dimensions of each element (rectangular geometry)
dimensions=[40 40;40 40;30 60;40 40;50 50;30 60;30 50;30 30];
        
% cross-section area of each element
A=zeros(nbars,1);
for i=1:nbars
    A(i)=dimensions(i,1)*dimensions(i,2);
end

% Cross-section inertia
I=zeros(nbars,1);
for i=1:nbars
    I(i)=1/12*dimensions(i,1)*dimensions(i,2)^3;
end

% coordinates of each node
coordxy=[0 -100;0 400;0 800;600 800;600 400;600 -100;1200 400;1200 -100]; 

%% Topology (connectivity)
ni=[1;2;3;4;5;2;5;7];
nf=[2;3;4;5;6;5;7;8];

%% Prescribed boudnary conditions [dof, displacement]
bc=[1 0;2 0;3 0;16 0;17 0;18 0;22 0;23 0;24 0];

supports=[1 "Fixed" "Fixed";
           2 "Fixed" "Fixed";
           3 "Fixed" "Fixed";
           4 "Fixed" "Fixed";
           5 "Fixed" "Fixed";
           6 "Fixed" "Fixed";
           7 "Fixed" "Fixed";
           8 "Fixed" "Fixed"];
       
%% Additional data (optional)
type_elem=[1 "Col";
           2 "Col";
           3 "Beam";
           4 "Col";
           5 "Col";
           6 "Beam";
           7 "Beam";
           8 "Col"];
       
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

%% Loads       
beams_LL=[1 -50; % Uniformly distributed loads over the beams
          2 -50;
          3 -50];

% Assignation of distributed loads on beams
qbary=zeros(nbars,2);
for i=1:beams
    qbary(elembeams(i),2)=beams_LL(i,2);
end
   
% Lateral equivalent seismic forces from a modal analysis. The number of
% forces must be equal to the number of floors
baseShearForces=[1000; % lower floor
                1500]; % upper floor
            
% Degrees of freedom over which each seismic force is applied (one for
% each seismic force)
dofSeismicForces=[4 7];

%% Plastic moments of each element's ends
Mp=[9680000 9680000;
    8490000 8490000;
    3363000 3276940;
    9490000 9490000;
    9680000 9680000;
    3363000 3976940;
    3363000 3976940;
    8490000 8490000]; %Kg-cm

% Height of each floor
hfloor=[400; 400];  
nfloors=length(hfloor);

%% PUSHOVER IN POSITIVE DIRECTION OF FORCES

[lambdaRight,pdriftDIRight,driftDIRight,defBasedDIRight,maxDispRight,...
 barPlasNodeRight]=Pushover2DFrames2(qbary,A,Mp,E,I,coordxy,ni,nf,...
supports,bc,baseShearForces,hfloor,dofSeismicForces,0.01,0.009);

%% PUSHOVER IN NEGATIVE DIRECTION OF FORCES

baseShearForces=-baseShearForces;
    
[lambdaLeft,pdriftDILeft,driftDILeft,defBasedDILeft,maxDispLeft,...
 barPlasNodeLeft]=Pushover2DFrames2(qbary,A,Mp,E,I,coordxy,ni,nf,...
supports,bc,baseShearForces,hfloor,dofSeismicForces,0.01,0.009);

%% Final results
SafetyFac=min([max(lambdaRight), max(lambdaLeft)])
pdriftDI=min([sum(pdriftDIRight)/nfloors,sum(pdriftDILeft)/nfloors])
driftDI=min([sum(driftDIRight)/nfloors,sum(driftDILeft)/nfloors])
dbDI=min([sum(defBasedDIRight)/nfloors,sum(defBasedDILeft)/nfloors])

Max_Displacement=max(max(maxDispLeft),max(maxDispRight))
