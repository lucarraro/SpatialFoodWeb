clear all; close all; clc

addpath('utilities');
load('100FW_nSp100.mat')
load('OCN.mat')

nSpecies=length(r_list.w1)+1;
nFW=length(fieldnames(r_list));

L(L==0)=100; % fix L=0 at the outlet
Q=width.*depth.*velocity;
V=width.*depth.*L;

kN_vec=[1e-5 1e-4 1e-3]; % equivalent to [10-2 10-1 10] kg_N/km2/day

nPar=3;
simParam=2*ones(nPar,5); 
simParam(2,1)=1; simParam(3,1)=3; % alternative nutrient load

y_thr=5e-7;
for ind_FW=1:nFW
    a_mat = A_list.(['w',num2str(ind_FW)]);
    bodymass = bodymass_list.(['w',num2str(ind_FW)]);
    r = r_list.(['w',num2str(ind_FW)]);
    
    for indPar=1:nPar
        
        fnam=['r_',num2str(ind_FW),'_',num2str(indPar),'.mat'];
        if not(isfile(['results_UWB/',fnam]))
            tmp=[];
            save(['results_UWB/',fnam],'tmp')
            
            kN_mean=kN_vec(simParam(indPar,1));
         
            parameters_local=v2struct(Q,V,kN_mean,As,r,a_mat);
            tic;
            y_mat=zeros(nSpecies,1); y_mat(:,1)=1e-4*ones(nSpecies,1);
            reltol=1; ind_time=1;
            while reltol>5e-3 % change in value after 10 days
                y0=y_mat(:,ind_time);
                tspan=[1+10*(ind_time-1)*86400 86400*10*ind_time];
                [t,y] = ODE_UWB(parameters_local,tspan,y0);
                y0=y(end,:)';
                y0(y0<y_thr)=0; % variables with density <y_thr are set to 0
                y_mat(:,ind_time+1)=y0;
                reltol=abs(y_mat(:,ind_time+1)-y_mat(:,ind_time))./y_mat(:,ind_time+1);
                reltol(reltol==Inf)=-Inf; reltol(isnan(reltol))=-Inf;
                ind_max=find(reltol==max(reltol)); ind_max=ind_max(1);
                node_max=floor(ind_max/nSpecies)+1; sp_max=mod(ind_max,nSpecies); sp_max(sp_max==0)=nSpecies;
                reltol=max(reltol);
                disp(sprintf('FW: %d  -  sim: %d  -  Elapsed time: %.2f s -  Sim time: %d d  - reltol: %.2e  -  species %d  -  node %d  - value = %.2e',...
                    ind_FW,indPar,toc,10*ind_time,reltol,sp_max,node_max,y_mat(ind_max,ind_time+1)))
                ind_time = ind_time+1;
            end
            y=y_mat(:,end);
            timeSim=10*ind_time; timeElapsed=toc;
            param=simParam(indPar,:);
            save(['results_UWB/',fnam],'y','timeSim','timeElapsed','param',...
                'kN_vec','y_thr','bodymass')  
        else
            disp(sprintf('FW: %d  -  sim: %d  -  done!',ind_FW,indPar))
        end
        
    end
    disp(' ')
end

%% prepare results for R
FW_vec=zeros(nFW,1);  nL_vec=FW_vec;
alphaDiv_vec=FW_vec;  alphaDivNorm_vec=FW_vec; 
nNutFeed_vec=FW_vec;

y_raw=zeros(nSpecies,nFW); 
FW_raw=zeros(1,nFW); nL_raw=zeros(1,nFW);

ind=1;
for indFW=1:nFW
    a_mat = A_list.(['w',num2str(indFW)]);
    r = r_list.(['w',num2str(indFW)]);
    for indSim=1
        fnam=['r_',num2str(indFW),'_',num2str(indSim),'.mat'];
        load(['results_UWB/',fnam])
                
        y_raw(:,ind)=y; FW_raw(ind)=indFW;
        nL_raw(ind)=param(1); 
        
        % y=5e-7  threshold density 
        y(y<5e-7)=0;
        PA=y>0; 
        dietMat = dietMatrix_list.(['w',num2str(indFW)]);
        nutFeed = nutrientFeeders_list.(['w',num2str(indFW)]);
       
        alphaDiv_vec(ind)=sum(PA); % remove nutrients
        nNutFeed_vec(ind)= sum(PA(nutFeed));
        if indSim==1; alphaDiv_med=median(sum(PA)); end
        alphaDivNorm_vec(ind) =  alphaDiv_vec(ind)/alphaDiv_med;

        FW_vec(ind)=indFW;

        nL_vec(ind)=param(1);
        ind=ind+1;
    end
end

save('utilities/resultsUWB_for_R.mat','alphaDiv_vec','alphaDivNorm_vec','nNutFeed_vec',...
    'FW_vec','nL_vec','FW_raw','y_raw','nL_raw')
