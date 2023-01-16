function binscatter(x,y,binX,binY,maxCount,ax,colmap)

if nargin==6
    colmap=flipud(copper);
end

for i=2:length(binX)
   midX(i-1)=0.5*(binX(i)+binX(i-1)); 
end
for i=2:length(binY)
   midY(i-1)=0.5*(binY(i)+binY(i-1)); 
end
mat=zeros(length(binY),length(binX));
for i=2:length(binX)
    for j=2:length(binY)
        mat(j-1,i-1)=length(find(x<=binX(i) & x>binX(i-1) & y<=binY(j) & y>binY(j-1)));
    end
end
mat(mat==0)=NaN;
pp=pcolor(binX,binY,mat); pp.EdgeColor='none'; colormap(ax,colmap); caxis([0 maxCount]); colorbar;
set(gca,'xlim',[min(binX) max(binX)],'ylim',[min(binY) max(binY)],'tickdir','out'); box off; 
