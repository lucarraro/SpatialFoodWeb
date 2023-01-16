clear all; close all; clc;

addpath('utilities')
load('utilities/delta_SFW_UWB.mat')

colmap.c1=[linspace(59/255,0,100)' linspace(153/255,0,100)' linspace(177/255,0,100)'];
colmap.c2=[linspace(159/255,0,100)' linspace(192/255,0,100)' linspace(149/255,0,100)'];
colmap.c3=[linspace(232/255,0,100)' linspace(164/255,0,100)' linspace(25/255,0,100)'];
colmap.c4=[linspace(245/255,0,100)' linspace(25/255,0,100)' linspace(28/255,0,100)'];

figure('units','centimeters','position',[0 0 30 20]); 
for pos=1:4
ax.(['a',num2str(pos)])=subplot(6,4,pos); 
binscatter(dfDelta.alphaDiv(dfDelta.pos==pos),dfDelta.Connectance(dfDelta.pos==pos),...
    [-80:5:0],[-0.1:0.025:0.25],1000,ax.(['a',num2str(pos)]),colmap.(['c',num2str(pos)])); 
hold on; plot([-80 0],[0 0],'--m'); 
set(gca,'xtick',[-80:20:0],'ytick',[-0.1:0.1:0.3],'ylim',[-0.1 0.3])
xlabel('\Delta Sp. Richness'); ylabel('\Delta Connectance')
ax.(['a',num2str(pos)])=subplot(6,4,4+pos); 
binscatter(dfDelta.alphaDiv(dfDelta.pos==pos),dfDelta.LinkDensity(dfDelta.pos==pos),...
    [-80:5:0],[-16:1:0],1000,ax.(['a',num2str(pos)]),colmap.(['c',num2str(pos)]));  
hold on; plot([-80 0],[0 0],'--m'); set(gca,'xtick',[-80:20:0],'ytick',[-16:4:0])
xlabel('\Delta Sp. Richness'); ylabel('\Delta Link Density')
ax.(['a',num2str(pos)])=subplot(6,4,8+pos); 
binscatter(dfDelta.alphaDiv(dfDelta.pos==pos),dfDelta.Modularity(dfDelta.pos==pos),...
    [-80:5:0],[-0.35:0.035:0.15],1000,ax.(['a',num2str(pos)]),colmap.(['c',num2str(pos)])); 
hold on; plot([-80 0],[0 0],'--m'); set(gca,'xtick',[-80:20:0],'ytick',[-0.4:0.2:0.2],'ylim',[-0.4 0.2])
xlabel('\Delta Sp. Richness'); ylabel('\Delta Modularity')
ax.(['a',num2str(pos)])=subplot(6,4,12+pos); 
binscatter(dfDelta.alphaDiv(dfDelta.pos==pos),dfDelta.Nestedness(dfDelta.pos==pos),...
    [-80:5:0],[-0.3:0.05:0.5],1000,ax.(['a',num2str(pos)]),colmap.(['c',num2str(pos)]));  
hold on; plot([-80 0],[0 0],'--m'); set(gca,'xtick',[-80:20:0],'ytick',[-0.4:0.2:0.6],'ylim',[-0.4 0.6])
xlabel('\Delta Sp. Richness'); ylabel('\Delta Nestedness')
ax.(['a',num2str(pos)])=subplot(6,4,16+pos); 
binscatter(dfDelta.alphaDiv(dfDelta.pos==pos),dfDelta.NicheOverlap(dfDelta.pos==pos),...
    [-80:5:0],[-0.2:0.05:0.6],1000,ax.(['a',num2str(pos)]),colmap.(['c',num2str(pos)])); 
hold on; plot([-80 0],[0 0],'--m'); set(gca,'xtick',[-80:20:0],'ytick',[-0.2:0.2:0.6],'ylim',[-0.2 0.6])
xlabel('\Delta Sp. Richness'); ylabel('\Delta Niche Overlap')
ax.(['a',num2str(pos)])=subplot(6,4,20+pos);
binscatter(dfDelta.alphaDiv(dfDelta.pos==pos),dfDelta.Omnivory(dfDelta.pos==pos),...
    [-80:5:0],[-0.7:0.05:0.3],1000,ax.(['a',num2str(pos)]),colmap.(['c',num2str(pos)]))
hold on; plot([-80 0],[0 0],'--m'); set(gca,'xtick',[-80:20:0],'ytick',[-0.8:0.4:0.4],'ylim',[-0.8 0.4])
xlabel('\Delta Sp. Richness'); ylabel('\Delta Omnivory')
end

% add histograms of fw metrics
figure('units','centimeters','position',[0 0 10 20]); 
for pos=1:4
ax.(['a',num2str(pos)])=subplot(6,4,pos); 
histogram(dfDelta.Connectance(dfDelta.pos==pos),[-0.1:0.025:0.25],'FaceColor',colmap.(['c',num2str(pos)])(1,:),...
    'Normalization','probability')
set(gca,'tickdir','out','view',[90 -90],'xlim',[-0.1 0.3],'ylim',[0 0.7])
ax.(['a',num2str(pos)])=subplot(6,4,4+pos); 
histogram(dfDelta.LinkDensity(dfDelta.pos==pos),[-16:1:0],'FaceColor',colmap.(['c',num2str(pos)])(1,:),...
    'Normalization','probability')
set(gca,'tickdir','out','view',[90 -90],'xlim',[-16 0],'ylim',[0 0.7])
ax.(['a',num2str(pos)])=subplot(6,4,8+pos); 
histogram(dfDelta.Modularity(dfDelta.pos==pos),[-0.35:0.035:0.15],'FaceColor',colmap.(['c',num2str(pos)])(1,:),...
    'Normalization','probability')
set(gca,'tickdir','out','view',[90 -90],'xlim',[-0.4 0.2],'ylim',[0 0.7])
ax.(['a',num2str(pos)])=subplot(6,4,12+pos); 
histogram(dfDelta.Nestedness(dfDelta.pos==pos),[-0.3:0.05:0.5],'FaceColor',colmap.(['c',num2str(pos)])(1,:),...
    'Normalization','probability')
set(gca,'tickdir','out','view',[90 -90],'xlim',[-0.4 0.6],'ylim',[0 0.7])
ax.(['a',num2str(pos)])=subplot(6,4,16+pos); 
histogram(dfDelta.NicheOverlap(dfDelta.pos==pos),[-0.2:0.05:0.6],'FaceColor',colmap.(['c',num2str(pos)])(1,:),...
    'Normalization','probability')
set(gca,'tickdir','out','view',[90 -90],'xlim',[-0.2 0.6],'ylim',[0 0.7])
ax.(['a',num2str(pos)])=subplot(6,4,20+pos);
histogram(dfDelta.Omnivory(dfDelta.pos==pos),[-0.7:0.05:0.3],'FaceColor',colmap.(['c',num2str(pos)])(1,:),...
'Normalization','probability')
set(gca,'tickdir','out','view',[90 -90],'xlim',[-0.8 0.4],'ylim',[0 0.7])

end
% add histograms of sp richness
figure('units','centimeters','position',[0 0 30 5]); 
for pos=1:4
    subplot(2,4,pos)
histogram(dfDelta.alphaDiv(dfDelta.pos==pos),[-80:5:0],'FaceColor',colmap.(['c',num2str(pos)])(1,:),...
    'Normalization','probability')
set(gca,'tickdir','out','xlim',[-80 0],'ylim',[0 0.4])
end



