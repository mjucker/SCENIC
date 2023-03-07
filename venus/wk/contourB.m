[RB,ZB,BM]=cutbimax(0);
nse=100;
figure;
contour(RB,ZB,BM,nse);colorbar;axis equal;
hold on
[tmp,ind]=min(abs(BM(:)-Bc));
k=0;
RBc=0;ZBc=0;BMc=0;
[im,jm]=size(BM);
for i=1:im
    for j=1:jm
        if(abs(BM(i,j)-BM(ind))<0.01)
            k=k+1;
            RBc(k)=RB(i,j);
            ZBc(k)=ZB(i,j);
            BMc(k)=BM(i,j);
        end
    end
end
plot(RBc,ZBc,'*m');
title('B(R,Z) contours');

[RP,ZP,PM]=cutbimaxPperp(0);
figure;hold on;
PM=2*PM/(B0^2)*100;
contour(RP,ZP,PM,nse);colorbar;%axis equal;
plot(RBc,ZBc,'*m');
axis equal;
title('\beta_{\perp}(R,Z) contours [%]');

figure;hold on;
colormap('jet');
set(gcf,'renderer','zbuffer')
surf(RP,ZP,PM);
shading interp;
axis('equal');
a=axis;a(1)=min(min(RP));a(2)=max(max(RP));axis(a);
%axis equal
xlabel('R');
ylabel('Z');
view(0,90);
title('\beta_{\perp} [%]','FontSize',16);

Hp=griddata(RP,ZP,PM,RBc,ZBc);
[ZBc,ix]=sort(ZBc);
plot3(RBc(ix),ZBc(ix),Hp(ix)+0.1,'w','MarkerSize',3);
h=axes('Position',[0.865 +0.55 0.03 0.4]);
bsc(2)=max(max(PM));bsc(1)=min(min(PM));
c=bsc(1):(bsc(2)-bsc(1))./32:bsc(2);
C=[c' c'];
pcolor([0 1],c,C);
shading interp;
set(h,'Xtick',[]);