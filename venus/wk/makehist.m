pos=load('positions.out');
th=linspace(0,2*pi,nbinspol+2);
es=linspace(0,1,nbinsrad+1);
[a,b] = find(pos(:,6)==0);
lst=pos(a,:);
sz = size(pos);
vec_to_keep = 1:sz(1);
vec_to_keep(a) = [];
ps=pos(vec_to_keep,:);

mHist=hist2d([mod(ps(:,2),2*pi), ps(:,1)],th,es);

nYBins = length(th);
nXBins = length(es);
vYLabel = 0.5*(th(1:(nYBins-1))+th(2:nYBins))/pi;
vXLabel = 0.5*(es(1:(nXBins-1))+es(2:nXBins));

figure;
pcolor(vXLabel,vYLabel,mHist);
colorbar
xlabel('s');
ylabel('theta [pi]');

figure;
vXEdge = linspace(1/2/(length(es)-1),1-1/2/(length(es)-1),length(es)-1);
vYEdge = linspace(-1/(length(th)-1),2-1/(length(th)-1),length(th)-1)*pi;
x=vXEdge'*cos(vYEdge);
y=vXEdge'*sin(vYEdge);
surf(x',y',mHist/ntot*nbinspol*nbinsrad*100);view(0,90)
colorbar
axis equal tight
title('Ptcles end of simulation');

mHistL=hist2d([mod(lst(:,2),2*pi), lst(:,1)],th,es);
figure;
vXEdge = linspace(1/2/(length(es)-1),1-1/2/(length(es)-1),length(es)-1);
vYEdge = linspace(-1/(length(th)-1),2-1/(length(th)-1),length(th)-1)*pi;
x=vXEdge'*cos(vYEdge);
y=vXEdge'*sin(vYEdge);
surf(x',y',mHistL/ntot*nbinspol*nbinsrad*100);view(0,90)
colorbar
axis equal tight
title('Lost ptcles');
colormap hsv

