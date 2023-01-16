clear all; close all; clc

addpath('utilities')
load('100FW_nSp100.mat')
load('OCN.mat')

nJob=input('Enter nJob: ');
nSpecies=length(r_list.w1)+1;
nFW=length(fieldnames(r_list));

L(L==0)=100; % fix L=0 at the outlet
nReach=length(downNode);
Q=width.*depth.*velocity;
V=width.*depth.*L;
Wt=W';

kN_vec=[1e-2 1e-1 1e-0]; % m-2/day
dispRate_vec=[1e-8 1e-7 1e-6];
downBias_vec=[1e4 1e1 1e-2];
vUptake_vec=[1e-6 1e-5 1e-4]; % upatake velocity [m/s]
LUscore_mat=ones(nReach,3); % flat / downstream / random
LUscore_mat(:,2)=1-0.99*distToOutlet/max(distToOutlet);

nPar=11;
simParam=2*ones(nPar,5); % 11 simulations, 5 parameters; default value is 2
simParam(2,1)=1; simParam(3,1)=3; % alternative nutrient load
simParam(4,2)=1; simParam(5,2)=3; % alternative nutrient distribution
simParam(6,3)=1; simParam(7,3)=3; % alternative dispersal rate
simParam(8,4)=1; simParam(9,4)=3; % alternative downstream bias
simParam(10,5)=1; simParam(11,5)=3; % alternative uptake velocity
simParam(12,3)=1; simParam(12,4)=1; % simulation with low dispersal, no dB

y_thr=5e-7;

for ind_FW=1:nFW
    a_mat = A_list.(['w',num2str(ind_FW)]);
    bodymass = bodymass_list.(['w',num2str(ind_FW)]);
    r = r_list.(['w',num2str(ind_FW)]);
    rng(ind_FW); LUscore_mat(:,3)=rand(nReach,1); % randomize nutrient distribution 
    for indPar=1:nPar
        
        fnam=['r_',num2str(ind_FW),'_',num2str(indPar),'.mat'];
        if not(isfile(['results/',fnam]))
            tmp=[];
            save(['results/',fnam],'tmp')
            
            kN_mean=kN_vec(simParam(indPar,1));
            LUscore=LUscore_mat(:,simParam(indPar,2));
            meanDispRate=dispRate_vec(simParam(indPar,3));
            downBias=downBias_vec(simParam(indPar,4));
            vUptake=vUptake_vec(simParam(indPar,5));
            
            kN=kN_mean*LUscore*max(A)/sum(As.*LUscore)/86400;
            dispersalRate=nSpecies*meanDispRate*(bodymass.^0.36)./sum(bodymass.^0.36);
            pD_vec = 0.5*(1 + exp(-downBias*bodymass));
            p_matrix=eval_p_matrix(nReach,nSpecies,pD_vec,downNode,V,depth,width);
            
            % simulation
            %if isDone(indSim)==0
            
            parameters=v2struct(Q,V,Wt,kN,As,r,a_mat,vUptake,depth,p_matrix,dispersalRate);
            tic;
            y_mat=zeros(nSpecies*nReach,1); y_mat(:,1)=1e-4*ones(nSpecies*nReach,1);
            reltol=1; ind_time=1;
            while reltol>5e-3 % change in value after 10 days
                y0=y_mat(:,ind_time);
                tspan=[1+10*(ind_time-1)*86400 86400*10*ind_time];
                [t,y] = ODE_full(parameters,tspan,y0);
                y0=y(end,:)';
                y0(y0<y_thr)=0; % variables with density <y_thr are set to 0
                y_mat(:,ind_time+1)=y0;
                reltol=abs(y_mat(:,ind_time+1)-y_mat(:,ind_time))./y_mat(:,ind_time+1);
                reltol(reltol==Inf)=-Inf; reltol(isnan(reltol))=-Inf;
                ind_max=find(reltol==max(reltol)); ind_max=ind_max(1);
                node_max=floor(ind_max/nSpecies)+1; sp_max=mod(ind_max,nSpecies); sp_max(sp_max==0)=nSpecies;
                reltol=max(reltol);
                disp(sprintf('Job: %d  -  FW: %d  -  sim: %d  -  Elapsed time: %.2f s -  Sim time: %d d  - reltol: %.2e  -  species %d  -  node %d  - value = %.2e',...
                    nJob,ind_FW,indPar,toc,10*ind_time,reltol,sp_max,node_max,y_mat(ind_max,ind_time+1)))
                ind_time = ind_time+1;
            end
            y=y_mat(:,end);
            timeSim=10*ind_time; timeElapsed=toc;
            param=simParam(indPar,:);
            save(['results/',fnam],'y','timeSim','timeElapsed','param','nJob',...
                'kN_vec','dispRate_vec','downBias_vec','vUptake_vec','LUscore_mat','y_thr','bodymass')
        else
            disp(sprintf('Job: %d  -  FW: %d  -  sim: %d  -  done!',nJob,ind_FW,indPar))
        end
        
    end
    disp(' ')
end

%% prepare results for R

FW_vec=zeros(nPar*nFW,1); dR_vec=FW_vec; dB_vec=FW_vec; nL_vec=FW_vec;
nT_vec=FW_vec; pos_vec=FW_vec; vUp_vec=FW_vec;
alphaDiv_vec=FW_vec;  alphaDivNorm_vec=FW_vec; 
nNutFeed_vec=FW_vec;

y_raw=zeros(nSpecies*nReach,nFW*nPar); 
FW_raw=zeros(1,nFW*nPar); nL_raw=zeros(1,nFW*nPar);
nT_raw=zeros(1,nFW*nPar); dR_raw=zeros(1,nFW*nPar);
dB_raw=zeros(1,nFW*nPar); vUp_raw=zeros(1,nFW*nPar);

downHW = A < median(A) & distToOutlet < median(distToOutlet);
upHW = A < median(A) & distToOutlet > median(distToOutlet);
downMain = A > median(A) & distToOutlet < median(distToOutlet);
midReach = A > median(A) & distToOutlet > median(distToOutlet);
position=upHW+ 2*midReach + 3*downHW + 4*downMain;

pos_raw=repmat(position',nSpecies,1); pos_raw=pos_raw(:);

ind=1;
for indFW=1:nFW
    a_mat = A_list.(['w',num2str(indFW)]);
    r = r_list.(['w',num2str(indFW)]);
    for indSim=1:nPar
        fnam=['r_',num2str(indFW),'_',num2str(indSim),'.mat'];
        load(['results/',fnam])
                
        y_raw(:,ind)=y; FW_raw(ind)=indFW;
        nL_raw(ind)=param(1); nT_raw(ind)=param(2); dR_raw(ind)=param(3);
        dB_raw(ind)=param(4); vUp_raw(ind)=param(5);
        
        % y=5e-7  threshold density 
        y(y<5e-7)=0;
        PA=reshape(y>0,nSpecies,nReach); 
        dietMat = dietMatrix_list.(['w',num2str(indFW)]);
        nutFeed = nutrientFeeders_list.(['w',num2str(indFW)]); 
        vUptake=vUptake_vec(param(5)); meanDispRate=dispRate_vec(param(3));
        dispersalRate=nSpecies*meanDispRate*(bodymass.^0.36)./sum(bodymass.^0.36);
        alphaDiv_vec((ind-1)*nReach+1:ind*nReach)=sum(PA)'; % remove nutrients
        nNutFeed_vec((ind-1)*nReach+1:ind*nReach)= sum(PA(nutFeed,:))';
        if indSim==1; alphaDiv_med=median(sum(PA)); end
        alphaDivNorm_vec((ind-1)*nReach+1:ind*nReach) =  alphaDiv_vec((ind-1)*nReach+1:ind*nReach)/alphaDiv_med;
        FW_vec((ind-1)*nReach+1:ind*nReach)=indFW;
        pos_vec((ind-1)*nReach+1:ind*nReach)=position;
        nL_vec((ind-1)*nReach+1:ind*nReach)=param(1);
        nT_vec((ind-1)*nReach+1:ind*nReach)=param(2);
        dR_vec((ind-1)*nReach+1:ind*nReach)=param(3);
        dB_vec((ind-1)*nReach+1:ind*nReach)=param(4);
        vUp_vec((ind-1)*nReach+1:ind*nReach)=param(5);
        ind=ind+1;
    end
end

save('utilities/results_for_R.mat','alphaDiv_vec','alphaDivNorm_vec','nNutFeed_vec',...
    'FW_vec','dR_vec','dB_vec','nL_vec','nT_vec','vUp_vec','pos_vec',...
    'FW_raw','y_raw','nL_raw','nT_raw','dR_raw','dB_raw','vUp_raw','pos_raw')

