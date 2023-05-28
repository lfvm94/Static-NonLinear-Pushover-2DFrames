function fbreduced=convectiveForceVector(rxbar,rybar,mbar,...
                                         nbars,nnodes,edof,ni,nf)

%------------------------------------------------------------------------
% Syntax:
% fbreduced=convectiveForceVector(rxbar,rybar,mbar,...
%                                 nbars,nnodes,edof,ni,nf)
%
%------------------------------------------------------------------------
% PURPOSE
%  To compute the equivalent forces at each restricted degree of freedom
%  of the structure in the global system of reference according to the
%  distributed forces rxbar,rybar,mbar at each end of each element due to
%  the application of the distributed loads
% 
% INPUT:  ni,nf:                are the initial and final nodes of each 
%                               element. Vectors of size: nbars x 1
%
%         edof:                 is vector containing the non-restricted
%                               DOF. Size: [Non-restric-DOF,1]
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
% OUTPUT: fbreduced:            equivalent forces at each degree of fredom
%                               of the global structure. Vector of size:
%                               [nedof x 1]
%
%------------------------------------------------------------------------
% LAST MODIFIED: L.F.Veduzco    2023-02-23
%                Faculty of Engineering
%                Autonomous University of Queretaro
%------------------------------------------------------------------------

fb=zeros(3*nnodes,1);
for i=1:nbars

    fb(3*ni(i)-2)=fb(3*ni(i)-2)-rxbar(i,2);
    fb(3*ni(i)-1)=fb(3*ni(i)-1)-rybar(i,2);
    fb(3*ni(i))=fb(3*ni(i))-mbar(i,2);
    fb(3*nf(i)-2)=fb(3*nf(i)-2)-rxbar(i,3);
    fb(3*nf(i)-1)=fb(3*nf(i)-1)-rybar(i,3);
    fb(3*nf(i))=fb(3*nf(i))-mbar(i,3);
    
end

ndof=length(edof);
fbreduced=zeros(ndof,1);

for j=1:ndof
    fbreduced(j)=fb(edof(j));
end
