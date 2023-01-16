function [t,y] = ODE_full(parameters,tspan,y0)

[Q,V,Wt,kN,As,r,a_mat,vUptake,depth,p_matrix,dispersalRate]=v2struct(parameters);
nSpecies=size(a_mat,1); nReach=length(Q);

% known term
phi_mat=zeros(nSpecies,nReach);
phi_mat(nSpecies,:)=(kN.*As./V)';
PHI=phi_mat(:);

% loss term
L_mat=zeros(nSpecies,nReach);
L_mat(1:nSpecies-1,:)=repmat(r-dispersalRate,1,nReach);
L_mat(nSpecies,:)=-vUptake./depth'-(Q./V)';
L=L_mat(:);

% interaction matrix
Amat=sparse(nSpecies*nReach,nSpecies*nReach);
for i=1:nReach
    ind=(i-1)*nSpecies+1:i*nSpecies;
    Amat(ind,ind)=a_mat;
end

% dispersal matrix
Dmat=sparse(nSpecies*nReach,nSpecies*nReach);
for sp=1:nSpecies-1
    ind_sp=sp:nSpecies:(nSpecies*nReach);
    Dmat(ind_sp,ind_sp)=p_matrix(:,:,sp)*dispersalRate(sp);
end
ind_N=nSpecies:nSpecies:nReach*nSpecies;
TransportMatrix=Wt.*repmat(Q',nReach,1); % only upstream transport (don't multiply by 86400 as a's, b's are already in s^-1)
TransportMatrix=sparse(TransportMatrix);
Dmat(ind_N,ind_N)  =  1./V.*TransportMatrix;
        
   
[t,y]=ode23(@(t,y)eqs(t,y),tspan,y0);  %'NonNegative',1:(N_reach*N_var),  ) ,,odeset('RelTol',1e-6,'AbsTol',1e-12) odeset('NonNegative',1:(N_reach*6) 

function dy=eqs(t,y)
    % y is the CONCENTRATION in the reach!! hence the division by V
    dy=zeros(nSpecies*nReach,1);
    
    dy = PHI + L.*y + Dmat*y + (Amat*y).*y;  
    

end

end

