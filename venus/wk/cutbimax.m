function [RB,ZB,BM]=cutbimax(plotter);
%PLOTS CROSS SECTIONAL CUTS GENERATED FROM VMEC OR CUR3DPPPL;
fid8=fopen('fort.45','rt');
[im,ncount]=fscanf(fid8,'%3i',3);
nper=1;
ntheta=im(3);
nphi=im(1);
ns=im(2);
ntotal=ntheta.*nphi.*ns;
[b,ncount]=fscanf(fid8,'%g %g %g %g',[4,inf]);
ncount;
numcoil=ncount./4;
fclose(fid8);
 bmint=min(b');
 bmaxt=max(b');
 bsc(1)=bmint(4);
 bsc(2)=bmaxt(4);
%clf;
kij=0;
for k=1:nphi;
if(plotter~=0),figure;end
 for i=1:ns;
  for j=1:ntheta;
   kij=j+(i-1)*ntheta+(k-1)*ns*ntheta;
   XB(i,j)=b(1,kij);
   YB(i,j)=b(2,kij);
   RB(i,j)=sqrt(b(1,kij).*b(1,kij)+b(2,kij).*b(2,kij));
   ZB(i,j)=b(3,kij);
   BM(i,j)=b(4,kij);
   CL(i,j)=0.0;
  end;
% end;
end;
if(plotter~=0)
colormap('copper');
colormap('default');
colormap('gray');
colormap('jet');
%surf(XB,YB,ZB,BM);
%mesh(RB,ZB,CL);
%shading interp;
%axis([2.5 4.5 -1.6 1.6]);
%axis([0.75 2.05 -0.8 0.8]);
%axis([4.0 11.5 -4.5 4.5]);
%axis([17.5 22.5 -3.0 3.0]);
%axis([5.5 11.5 -3.0 3.0]);
%xlabel('X');
%ylabel('Y');
%zlabel('Z');
set(gcf,'renderer','zbuffer')
surf(RB,ZB,BM);
xlabel('R');
ylabel('Z');
view(0,90);
axis('equal');
%title('Mod-B^2 distribution in Sphellamak');
%title('Field line bending distribution in 3-period QAS device');
%title('Normal curvature distribution in 3-period QAS device');
%title('j\cdot B/B^2 distribution in 3-period QAS device');
%title('Geodesic curvature distribution in 3-period QAS device');
%title('Local shear distribution in 3-period QAS device');
%title('Mod-B^2 distribution in 3-period QAS device');
%title('Mod-B^2 distribution in 6-period {\it J}-optimised QPS device');
%title('Total pressure distribution in 2-period QAS system');
title('Mod-B distribution','FontSize',16);
%title('Total p_{\_|\_} distribution in a Tokamak');
%title('Total p_{||} distribution in a Tokamak');
%title('Normal curvature distribution in ST/Sphellamak');
%if (k==1)text(4.2,-4.50,'\phi=0');
%elseif (k==2)text(7.5,-4.50,'\phi=\pi/6');
%elseif (k==2)text(18.5,-2.40,'\phi=\pi/12');
%elseif (k==2)text(4.2,-4.50,'\phi=\pi/4');
%elseif (k==3)text(7.5,-4.50,'\phi=\pi/3');
%elseif (k==3)text(18.5,-2.40,'\phi=\pi/6');
%elseif (k==3)text(4.2,-4.50,'\phi=\pi/2');
%end;
%text(4.2,4.6,'2/3 p_{\_|\_}+1/3 p_{||}')
%title('Mod-B distribution in Configuration with Inclined Coils');
h=axes('Position',[0.865 +0.55 0.03 0.4]);
%bsc=[min(min(BM)) max(max(BM))];
c=bsc(1):(bsc(2)-bsc(1))./32:bsc(2);
C=[c' c'];
pcolor([0 1],c,C);
shading interp;
set(h,'Xtick',[]);
end;
end;