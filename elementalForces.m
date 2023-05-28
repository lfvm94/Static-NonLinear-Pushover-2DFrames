function fe=elementalForces(ke,ue,rxbar,rybar,mbar,i)
%------------------------------------------------------------------------
% Syntax:
% fe=elementalForces(ke,ue,rxbar,rybar,mbar,i)
%
%------------------------------------------------------------------------
% PURPOSE
%  To compute the the equivalent forces at each end of each element i
% according to the distributed forces rxbar, rybar, mbar due to
% the application of the distributed loads
% 
% INPUT:  ke:                   stiffness matrix of the element 6x6
%
%         ue:                   Displacement vector of the element 6x1
%
%         rxbar:                Distributed loads at end of the element
%                               in the local x direction
%                               due to the distributed loads. Vector 6x2
%
%         rybar:                Distributed loads at end of the element
%                               in the local y direction
%                               due to the distributed loads. Vector 6x2
%
%         mbar:                 Ditributed bending moments at the ends of
%                               the element due to the
%                               distributed loads. Vector 6x2
%
%         i:                    is the index of each element according to
%                               the topology established initially
%
% OUTPUT: fe:                   equivalent elemental forces. Vector 6x1
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-23
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

fbe=[rxbar(i,2);
     rybar(i,2);
     mbar(i,2);
     rxbar(i,3);
     rybar(i,3);
     mbar(i,3)];
 
fe=ke*ue+fbe;
