function Ce=plantml(ex,ey,c)
% Ce=plantml(ex,ey,c)
%-------------------------------------------------------------
% PURPOSE
%  Compute the quantity: Ce=c*int(N^T*N)dA
%
% INPUT:  ex,ey;       Element coordinates
%	
%	  c
%
% OUTPUT: Ce :      Matix 3 x 3
%-------------------------------------------------------------


Area=1/2*det([ones(3,1) ex' ey']);



L1=[0.5 0 0.5];
L2=[0.5 0.5 0];
L3=[0 0.5 0.5];

NtN=zeros(3);


for i=1:3
	NtN=NtN+1/3*[L1(i)^2     L1(i)*L2(i) L1(i)*L2(i)
		  		 L2(i)*L1(i) L2(i)^2     L2(i)*L3(i)
		 		 L3(i)*L1(i) L3(i)*L2(i) L3(i)^2];
end

Ce=NtN*Area*c;