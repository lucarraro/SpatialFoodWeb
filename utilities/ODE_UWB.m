function [t,y] = ODE_local(parameters_local,tspan,y0)

[Q,V,kN,As,r,a_mat]=v2struct(parameters_local);
nSpecies=size(a_mat,1); 
V = sum(V);
Q = max(Q);
A=sum(As);

% known term
PHI=zeros(nSpecies,1);
PHI(nSpecies)=kN*A/V;

% loss term
L=zeros(nSpecies,1);
L(1:nSpecies-1)=r;
L(nSpecies)=-Q/V;

% interaction matrix
Amat=a_mat;

   
[t,y]=ode23(@(t,y)eqs(t,y),tspan,y0,odeset('NonNegative',1:nSpecies));  %'NonNegative',1:(N_reach*N_var),  ) ,,odeset('RelTol',1e-6,'AbsTol',1e-12) odeset('NonNegative',1:(N_reach*6) 

function dy=eqs(t,y)
    % y is the CONCENTRATION in the reach!! hence the division by V
    dy=zeros(nSpecies,1);
    
    dy = PHI + L.*y + (Amat*y).*y;  
    

end

end