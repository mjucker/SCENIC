clear manual_min manual_max intervals manual_minE manual_maxE intervalsE
clear ctr_high ctrl_highE freq mi1_SI Qi1_SI Bres1
clear fif v niT njT nkT nperiod Njk v2 tb2 crx crz Btot coorRes1 coorRes1s
clear fid s_hp cR_v cZ_v NbDistr k_lkp AP_nb Pres_Sp Ppla_rv Ppla_iv Ppla_fi E3D_rv E3D_iv
clear B3D_rv B3D_iv A3D_rv A3D_iv Ppla3D_rv Ppla3D_iv PplaSp_rv PplaSp_iv
clear eps_nn_rv eps_nn_iv eps_bn_rv eps_bn_iv eps_par_rv eps_par_iv
clear kpar_v
clear A3_s B3_s E3D E3_s1 E3_s2 Emavg Epavg MAX_E MIN_E NelR OutSelect P3_s PSp_s Ppla3D PplaSp cR cZ
clear ctr_highE eR electrons eps_bn_s eps_nn_s eps_par_s fast_ions first_species hx hy
clear kpar kpar_iv kpar_rv kpar_s n_ant second_species temp_PSp temp_kpar x3 x3_s xb y3 y3_s yb z3 z3_s zb

clear BM BMc a b indt it rb rp C Hp PM RB RP ZB ZP bm bsc c cl delta fid8 h im ix jm pm ppm rpp tmp zp zpp

manual_min =  0.0e0;
manual_max =  4.0e2;
intervals  = 200;
manual_minE = 0.0e0;
manual_maxE = 1.0e0;
intervalsE  = 200;

for i = 1:intervals+1
  ctr_high(i) = ((manual_max-manual_min)/intervals)*(i-1)+manual_min;
end

for i = 1:intervalsE+1
  ctr_highE(i) = ((manual_maxE-manual_minE)/intervalsE)*(i-1)+manual_minE;
end

%freq=2.04e8
%freq=1.74e8
%freq=2.03e8
freq=3.1e8
%freq=2.642e8;
%freq=2.0651e8;
%B0=2.89;
%B0=2.8871;
%mi1_SI=1.67e-27;
%mi1_SI=1.617*1.67e-27;
%mi1_SI=3*1.67e-27;
mi1_SI=1*1.67e-27;
Qi1_SI=2*1.6e-19;
Qi1_SI=1.6e-19;

Bres1 = mi1_SI*freq/(Qi1_SI*B0)
%Bres1 = 1.2367*mi1_SI*freq/(1.6e-19*B0)


% --------------------------------------------------------------------------
% Read the equilibrium from fort.37.DAT

fif=fopen('fort.37.DAT.N','r','l');

fread(fif,1,'uint');
v=fread(fif,30,'char');
fread(fif,1,'uint');

fread(fif,1,'uint');
v=fread(fif,4,'int');
fread(fif,1,'uint');

niT     = v(1);
njT     = v(2);
nkT     = v(3);
nperiod = v(4);
clear v;
Njk = njT * nkT;

disp(sprintf('TERPSICHORE grid: (ni=%d,  nj=%d,  nk=%d)',niT,njT,nkT))


fread(fif,1,'uint');
v=fread(fif,niT*7,'double');
fread(fif,1,'uint');

v2 = reshape(v,niT,7);

fread(fif,1,'uint');
v=fread(fif,niT*Njk*11,'double');
fread(fif,1,'uint');

v2 = reshape(v,Njk,niT,11);

tb2   = v2(:,:,11);
crx   = v2(:,:,9);
crz   = v2(:,:,10);

Btot=0;
for i = 1:Njk
  Btot = Btot+tb2(i,1);
end
Btot = Btot/Njk;
tb2 = sqrt(tb2/Btot);

indt=0;
it=1;
for i = 1:niT
for ip = 1:njT
ind = (it-1)*njT + ip;
if (ip==1 && (tb2(ind,i)-Bres1)*(tb2(it*njT,i)-Bres1)<0)
  indt=indt+1;
  coorRes1(indt,1)=(crx(ind,i)*(tb2(it*njT,i)-Bres1)-crx(it*njT,i)*(tb2(ind,i)-Bres1))/(tb2(it*njT,i)-tb2(ind,i));
  coorRes1(indt,2)=1;
  coorRes1(indt,3)=(crz(ind,i)*(tb2(it*njT,i)-Bres1)-crz(it*njT,i)*(tb2(ind,i)-Bres1))/(tb2(it*njT,i)-tb2(ind,i));
  indx(indt,1)=ind;
  indx(indt,2)=i;
end
if (ip~=1 && (tb2(ind,i)-Bres1)*(tb2(ind-1,i)-Bres1)<0)
  indt=indt+1;
  coorRes1(indt,1)=(crx(ind,i)*(tb2(ind-1,i)-Bres1)-crx(ind-1,i)*(tb2(ind,i)-Bres1))/(tb2(ind-1,i)-tb2(ind,i));
  coorRes1(indt,2)=1;
  coorRes1(indt,3)=(crz(ind,i)*(tb2(ind-1,i)-Bres1)-crz(ind-1,i)*(tb2(ind,i)-Bres1))/(tb2(ind-1,i)-tb2(ind,i));
end
end
end


%coorRes1(:,4)=coorRes1(:,1)+coorRes1(:,3);

  coorRes1s=sortrows(coorRes1,3);
%id = 0;
%iu = 0;
%for i = 1:indt
%  if (coorRes1(i,3)<0)
%    id=id+1;
%    coorRes1d(id,:)=coorRes1(i,:);
%  else
%    iu=iu+1;
%    coorRes1u(iu,:)=coorRes1(i,:);
%  end
%end
%coorRes2d=sortrows(coorRes1d,1);
%coorRes2u=sortrows(coorRes1u,1);


clear v v2;

fclose(fif);

%----------------------------------------------------

fid=fopen('Ppla3D.dat','r',endian);

fread(fid,1,'uint');
v=fread(fid,3,'int');
fread(fid,1,'uint');
NelR = v(1);
njT  = v(2);
nkT  = v(3);
clear v;
Njk = njT*nkT;

fread(fid,1,'uint');
s_hp=fread(fid,NelR,'double');
fread(fid,1,'uint');

fread(fid,1,'uint');
cR_v=fread(fid,NelR*Njk,'double');
fread(fid,1,'uint');
fread(fid,1,'uint');
cZ_v=fread(fid,NelR*Njk,'double');
fread(fid,1,'uint');

fread(fid,1,'uint');
OutSelect=fread(fid,8,'int');
fread(fid,1,'uint');

fread(fid,1,'uint');
NbDistr=fread(fid,1,'uint');
fread(fid,1,'uint');

k_lkp(1)=1;
for i = 2:NbDistr-1
fread(fid,1,'uint');
k_lkp(i)=fread(fid,1,'int');
fread(fid,1,'uint');
end
k_lkp(NbDistr)=7;
k_lkp(NbDistr+1)=8;

fread(fid,1,'uint');
AP_nb=fread(fid,1,'uint');
fread(fid,1,'uint');

for i = 1:AP_nb
fread(fid,1,'uint');
Pres_Sp(i)=fread(fid,1,'int');
fread(fid,1,'uint');
end

fread(fid,1,'uint');
Ppla_rv=fread(fid,1,'double');
fread(fid,1,'uint');

fread(fid,1,'uint');
Ppla_iv=fread(fid,1,'double');
fread(fid,1,'uint');

if(Pres_Sp(AP_nb-1)==3)
  fread(fid,1,'uint');
  Ppla_fi=fread(fid,1,'double');
  fread(fid,1,'uint');
end

if(Pres_Sp(AP_nb)==4)
  fread(fid,1,'uint');
  Ppla_fi=fread(fid,1,'double');
  fread(fid,1,'uint');
end

if(OutSelect(1)==1)
  fread(fid,1,'uint');
  E3D_rv=fread(fid,3*NelR*Njk,'double');
  fread(fid,1,'uint');
  fread(fid,1,'uint');
  E3D_iv=fread(fid,3*NelR*Njk,'double');
  fread(fid,1,'uint');
end

if(OutSelect(2)==1)
  fread(fid,1,'uint');
  B3D_rv=fread(fid,3*NelR*Njk,'double');
  fread(fid,1,'uint');
  fread(fid,1,'uint');
  B3D_iv=fread(fid,3*NelR*Njk,'double');
  fread(fid,1,'uint');
end

if(OutSelect(3)==1)
  fread(fid,1,'uint');
  A3D_rv=fread(fid,4*NelR*Njk,'double');
  fread(fid,1,'uint');
  fread(fid,1,'uint');
  A3D_iv=fread(fid,4*NelR*Njk,'double');
  fread(fid,1,'uint');
end

if(OutSelect(4)==1)
  fread(fid,1,'uint');
  Ppla3D_rv=fread(fid,NelR*Njk,'double');
  fread(fid,1,'uint');
  fread(fid,1,'uint');
  Ppla3D_iv=fread(fid,NelR*Njk,'double');
  fread(fid,1,'uint');

  PplaSp_rv=zeros(AP_nb,NelR*Njk);
  PplaSp_iv=zeros(AP_nb,NelR*Njk);
  for i = 1:AP_nb
    fread(fid,1,'uint');
    PplaSp_rv(i,:,:)=fread(fid,NelR*Njk,'double');
    fread(fid,1,'uint');
    fread(fid,1,'uint');
    PplaSp_iv(i,:,:)=fread(fid,NelR*Njk,'double');
    fread(fid,1,'uint');
  end
end

if(OutSelect(5)==1)
  fread(fid,1,'uint');
  eps_nn_rv=fread(fid,NelR*Njk,'double');
  fread(fid,1,'uint');
  fread(fid,1,'uint');
  eps_nn_iv=fread(fid,NelR*Njk,'double');
  fread(fid,1,'uint');
end

if(OutSelect(6)==1)
  fread(fid,1,'uint');
  eps_bn_rv=fread(fid,NelR*Njk,'double');
  fread(fid,1,'uint');
  fread(fid,1,'uint');
  eps_bn_iv=fread(fid,NelR*Njk,'double');
  fread(fid,1,'uint');
end

if(OutSelect(7)==1)
  fread(fid,1,'uint');
  eps_par_rv=fread(fid,NelR*Njk,'double');
  fread(fid,1,'uint');
  fread(fid,1,'uint');
  eps_par_iv=fread(fid,NelR*Njk,'double');
  fread(fid,1,'uint');
end

if(OutSelect(8)==1)
  fread(fid,1,'uint');
  n_ant=fread(fid,1,'int');
  fread(fid,1,'uint');
  kpar_rv=zeros(NbDistr,NelR*Njk);
  kpar_iv=zeros(NbDistr,NelR*Njk);
  for i = 1:NbDistr+1
    fread(fid,1,'uint');
    kpar_v(i,:,:)=fread(fid,NelR*Njk,'double');
    fread(fid,1,'uint');
  end
end

fclose(fid);

clear E3D_v B3D_v A3D_v Ppla3D_v eps_nn_v eps_bn_v eps_par_v

if(OutSelect(1)==1)
  E3D_v = E3D_rv + j*E3D_iv;  clear E3D_rv E3D_iv;
end
if(OutSelect(2)==1)
  B3D_v = B3D_rv + j*B3D_iv;  clear B3D_rv B3D_iv;
end
if(OutSelect(3)==1)
  A3D_v = A3D_rv + j*A3D_iv;  clear A3D_rv A3D_iv;
end
if(OutSelect(4)==1)
  Ppla3D_v = Ppla3D_rv + j*Ppla3D_iv;  clear Ppla3D_rv Ppla3D_iv;
  PplaSp_v = PplaSp_rv + j*PplaSp_iv;  clear PplaSp_rv PplaSp_iv;
end
if(OutSelect(5)==1)
  eps_nn_v = eps_nn_rv + j*eps_nn_iv; clear eps_nn_rv eps_nn_iv
end
if(OutSelect(6)==1)
  eps_bn_v = eps_bn_rv + j*eps_bn_iv; clear eps_bn_rv eps_bn_iv
end
if(OutSelect(7)==1)
  eps_par_v = eps_par_rv + j*eps_par_iv; clear eps_par_rv eps_par_iv
end

clear cR cZ E3D B3D A3D Ppla3D PplaSp eps_nn eps_bn eps_par kpar


cR  = reshape(cR_v, NelR,njT,nkT);
cZ  = reshape(cZ_v, NelR,njT,nkT);
if(OutSelect(1)==1)
  E3D = reshape(E3D_v,3,NelR,njT,nkT);
end
if(OutSelect(2)==1)
  B3D = reshape(B3D_v,3,NelR,njT,nkT);
end
if(OutSelect(3)==1)
  A3D = reshape(A3D_v,4,NelR,njT,nkT);
end
if(OutSelect(4)==1)
  Ppla3D = reshape(Ppla3D_v,NelR,njT,nkT);
  PplaSp = reshape(PplaSp_v,AP_nb,NelR,njT,nkT);
end
if(OutSelect(5)==1)
  eps_nn = reshape(eps_nn_v,NelR,njT,nkT);
end
if(OutSelect(6)==1)
  eps_bn = reshape(eps_bn_v,NelR,njT,nkT);
end
if(OutSelect(7)==1)
  eps_par = reshape(eps_par_v,NelR,njT,nkT);
end
if(OutSelect(8)==1)
  kpar = reshape(kpar_v,NbDistr+1,NelR,njT,nkT);
end

clear Ppla3D_v E3D_v B3D_v A3D_v cR_v cZ_v PplaSp_v eps_nn_v eps_bn_v eps_par_v kpar_v;



% Close poloidally

njT = njT+1;
for eR = 1:NelR
  for ip = 1:nkT
    cR(eR,njT,ip)  =  cR(eR,1,ip);
    cZ(eR,njT,ip)  =  cZ(eR,1,ip);
    if(OutSelect(1)==1)
      E3D(:,eR,njT,ip) = E3D(:,eR,1,ip);
    end
    if(OutSelect(2)==1)
      B3D(:,eR,njT,ip) = B3D(:,eR,1,ip);
    end
    if(OutSelect(3)==1)
      A3D(:,eR,njT,ip) = A3D(:,eR,1,ip);
    end
    if(OutSelect(4)==1)
      Ppla3D(eR,njT,ip) = Ppla3D(eR,1,ip);
      PplaSp(:,eR,njT,ip) = PplaSp(:,eR,1,ip);
    end
    if(OutSelect(5)==1)
      eps_nn(eR,njT,ip) = eps_nn(eR,1,ip);
    end
    if(OutSelect(6)==1)
      eps_bn(eR,njT,ip) = eps_bn(eR,1,ip);
    end
    if(OutSelect(7)==1)
      eps_par(eR,njT,ip) = eps_par(eR,1,ip);
    end
    if(OutSelect(8)==1)
      kpar(:,eR,njT,ip) = kpar(:,eR,1,ip);
    end
  end 
end

clear x3 y3 z3

for eR = 1:NelR
  for it = 1:njT
  for ip = 1:nkT
%    x3(eR,it,ip) = cR(eR,it,ip) .* cos((ip-1)/nkT*(2*pi));
%    y3(eR,it,ip) = cR(eR,it,ip) .* sin((ip-1)/nkT*(2*pi));
    x3(eR,it,ip) = cR(eR,it,ip);
    y3(eR,it,ip) = cR(eR,it,ip);
    z3(eR,it,ip) = cZ(eR,it,ip);
  end 
  end
end


x3_s  = zeros(NelR,njT);
y3_s  = zeros(NelR,njT);
z3_s  = zeros(NelR,njT);
E3_s1  = zeros(NelR,njT);
E3_s2  = zeros(NelR,njT);
B3_s  = zeros(NelR,njT);
A3_s  = zeros(NelR,njT);
P3_s  = zeros(NelR,njT);
PSp_s = zeros(AP_nb,NelR,njT);
eps_nn_s = zeros(NelR,njT);
eps_bn_s = zeros(NelR,njT);
eps_par_s = zeros(NelR,njT);
kpar_s = zeros(NbDistr+1,NelR,njT);



% Definition of the components under investigation

for eR = 1:NelR

x3_s(eR,:)  = x3 (eR,:,1);
y3_s(eR,:)  = y3 (eR,:,1);
z3_s(eR,:)  = z3 (eR,:,1);

if(OutSelect(1)==1)
  E3_s1(eR,:) = E3D(1,eR,:,1) + complex(0,1)*E3D(2,eR,:,1);
  E3_s2(eR,:) = E3D(1,eR,:,1) - complex(0,1)*E3D(2,eR,:,1);
  %E3_s(eR,:) = E3D(2,eR,:,1);
end
if(OutSelect(2)==1)
  B3_s(eR,:) = B3D(2,eR,:,1);
end
if(OutSelect(3)==1)
  A3_s(eR,:) = A3D(2,eR,:,1);
end
if(OutSelect(4)==1)
  P3_s(eR,:) = Ppla3D(eR,:,1);
  PSp_s(:,eR,:) = PplaSp(:,eR,:,1);
end
if(OutSelect(5)==1)
  eps_nn_s(eR,:) = eps_nn(eR,:,1);
end
if(OutSelect(6)==1)
  eps_bn_s(eR,:) = eps_bn(eR,:,1);
end
if(OutSelect(7)==1)
  eps_par_s(eR,:) = eps_par(eR,:,1);
end
if(OutSelect(8)==1)
  kpar_s(:,eR,:) = kpar(:,eR,:,1);
end

end

%P3_line = zeros(2*NelR);
%half=round(njT/2+1)

%for eR = 1:NelR
%P3_line(eR) = P3_s(NelR+1-eR,half);
%P3_line(NelR+eR) = P3_s(eR,1);
%end

%figure
%plot(imag(P3_line))

% Fill the hole on the axis

%tmp1D = zeros(1,njT);
%tmp1D(:) = sum(x3_s(1,:))/njT;
%x3_s = '[tmp1D' x3_s']';
%tmp1D(:) = sum(y3_s(1,:))/njT;
%y3_s = '[tmp1D' y3_s']';
%tmp1D(:) = sum(z3_s(1,:))/njT;
%z3_s = '[tmp1D' z3_s']';
%tmp1D(:) = sum(E3_s(1,:))/njT;
%E3_s = '[tmp1D' E3_s']';
%tmp1D(:) = sum(P3_s(1,:))/njT;
%P3_s = '[tmp1D' P3_s']';


%E3_s(NelR-10:NelR+1,:) = 0;

% Border
xb = zeros(1,njT);
yb = zeros(1,njT);
zb = zeros(1,njT);
xb(:) = x3_s(NelR,:);
zb(:) = z3_s(NelR,:);

MAX_E=max(max(real(E3_s1(:))),max(real(E3_s2(:))))
MIN_E=min(min(real(E3_s1(:))),min(real(E3_s2(:))))

% Visualisation
% Available data: E3_s, B3_s, A3_s, P3_s, eps_nn_s, eps_bn_s, eps_par_s

figure
set(gcf,'renderer','zbuffer');
%subplot(2,3,1)
%mesh(x3_s,z3_s,real(E3_s)+imag(E3_s2));
%contour(x3_s,z3_s,abs(E3_s1),200);
contour(x3_s(1:NelR-wind,:),z3_s(1:NelR-wind,:),abs(E3_s1(1:NelR-wind,:)),200);
%contour(x3_s,z3_s,abs(E3_s),ctr_highE);
%colormap autumn
axis equal;
%axis([2.75 4.75 -0.8 0.8]) 
%axis([5 11 -3 3 min(min(min(real(E3_s(:,:,:))))) max(max(max(real(E3_s(:,:,:)))))]) 
%axis([5 11 -3 3 -0.1 0.1]) 
%shading interp;
view(0,90);
title('E_+');
hold on;
plot3(xb,zb,yb,'k',coorRes1s(:,1),coorRes1s(:,3),coorRes1s(:,2),'k-');
%plot3(xb,zb,yb,'k',coorRes2u(:,1),coorRes2u(:,3),coorRes2u(:,2),'k-',coorRes2d(:,1),coorRes2d(:,3),coorRes2d(:,2),'k-');
% 
 figure
 set(gcf,'renderer','zbuffer');
% %subplot(2,3,1)
% %mesh(x3_s,z3_s,real(E3_s)+imag(E3_s2));
 contour(x3_s(1:NelR-wind,:),z3_s(1:NelR-wind,:),real(E3_s2(1:NelR-wind,:)),200);
% %contour(x3_s,z3_s,abs(E3_s),ctr_highE);
 axis equal;
% %axis([2.75 4.75 -0.8 0.8]) 
% %axis([5 11 -3 3 min(min(min(real(E3_s(:,:,:))))) max(max(max(real(E3_s(:,:,:)))))]) 
% %axis([5 11 -3 3 -0.1 0.1]) 
% %shading interp;
 view(0,90);
 title('E_-');
 hold on;
 plot3(xb,zb,yb,'k',coorRes1s(:,1),coorRes1s(:,3),coorRes1s(:,2),'k-');
% %plot3(xb,zb,yb,'k',coorRes2u(:,1),coorRes2u(:,3),coorRes2u(:,2),'k-',coorRes2d(:,1),coorRes2d(:,3),coorRes2d(:,2),'k-');
% figure;
%  for i=1:NelR,Epavg(i)=mean(E3_s1(i,:));Emavg(i)=mean(E3_s2(i,:));end;
%  if(Bres1<=1)
%     plot(x3_s(:,1),abs(Epavg./Emavg))
%  else
%    plot(x3_s(:,floor(njT/2)),abs(Epavg./Emavg))  
%  end
%  title('|E_+/E_-|');


figure
set(gcf,'renderer','zbuffer');
%subplot(2,3,4)
%mesh(x3_s,z3_s,real(P3_s));
%[c,ch,cf]=contourf(x3_s,z3_s,real(P3_s),100);
%[c,ch,cf]=contourf(x3_s,z3_s,imag(P3_s),ctr_high);
%contour(x3_s,z3_s,real(P3_s),100);
contour(x3_s(1:NelR-wind,:),z3_s(1:NelR-wind,:),imag(P3_s(1:NelR-wind,:)),100);
%for i = 1:size(ch,1)
%set(ch(i),'LineStyle','none')
%end
axis equal;
%axis([2.75 4.75 -0.8 0.8]) 
%axis([1.5 4.5 -2.1 2.1 min(min(min(real(P3_s(:,:,:))))) max(max(max(real(P3_s(:,:,:)))))]) 
%shading interp;
view(0,90);
hold on;
plot3(xb,zb,yb,'k',coorRes1s(:,1),coorRes1s(:,3),coorRes1s(:,2),'k-');
%plot3(xb,zb,yb,'k',coorRes2u(:,1),coorRes2u(:,3),coorRes2u(:,2),'k-',coorRes2d(:,1),coorRes2d(:,3),coorRes2d(:,2),'k-');
hx=xlabel('r [m]');
hy=ylabel('z [m]');
set(gca,'Fontsize',[18]);
set(hx,'Fontsize',[24]);
set(hy,'Fontsize',[24]);


%figure
%mesh(x3_s,z3_s,imag(P3_s))

%break;


% Visualisation of power per species and evaluation of k||

for i = 1:AP_nb
temp_PSp(:,:)=imag(PSp_s(i,:,:));
figure
set(gcf,'renderer','zbuffer');
if(mean(temp_PSp(:))>1e-20),contour(x3_s(1:NelR-wind,:),z3_s(1:NelR-wind,:),temp_PSp(1:NelR-wind,:),200);end
axis equal;
view(0,90);
hold on;
plot3(xb,zb,yb,'k')

if (Pres_Sp(i)==1)
 title('Electrons')
 electrons=sum(temp_PSp(:))/sum(imag(PSp_s(:)))*100
end
if (Pres_Sp(i)==2)
 title('First Species')
 first_species=sum(temp_PSp(:))/sum(imag(PSp_s(:)))*100
end
if (Pres_Sp(i)==3)
 title('Second Species')
 second_species=sum(temp_PSp(:))/sum(imag(PSp_s(:)))*100
end
if (Pres_Sp(i)==4)
 title('Fast ions')
 fast_ions=sum(temp_PSp(:))/sum(imag(PSp_s(:)))*100
end

end

if(OutSelect(8)==1)
    
for i = 1:NbDistr+1
temp_kpar(:,:)=kpar_s(i,:,:);
[a,b]=find(temp_kpar>100);
temp_kpar(a,b)=0;
figure
set(gcf,'renderer','zbuffer');
%contour(x3_s,z3_s,real(temp_kpar),200);
surf(x3_s,z3_s,real(temp_kpar),'LineStyle','none');
axis equal;
view(0,90);
hold on;
plot3(xb,zb,yb,'k')

if (k_lkp(i)==1)
    if(a)
        title(['k_{||} (Initial, n_{ant}=',num2str(n_ant),', max=100)'])
    else
        title(['k_{||} (Initial, n_{ant}=',num2str(n_ant),')']) 
    end
end
if (k_lkp(i)==2)
    if(a)
        title(['k_{||} (Electrons parallel, n_{ant}=',num2str(n_ant),', max=100)'])
    else
        title(['k_{||} (Electrons parallel, n_{ant}=',num2str(n_ant),')'])
    end
end
if (k_lkp(i)==3)
    if(a)
        title(['k_{||} (Ions 1 parallel, n_{ant}=',num2str(n_ant),', max=100)'])
    else
        title(['k_{||} (Ions 1 parallel, n_{ant}=',num2str(n_ant),')'])
    end
end
if (k_lkp(i)==4)
    if(a)
        title(['k_{||} (Ions 2 parallel, n_{ant}=',num2str(n_ant),', max=100)'])
    else
       title(['k_{||} (Ions 2 parallel, n_{ant}=',num2str(n_ant),')']) 
    end
end
if (k_lkp(i)==5)
    if(a)
        title(['k_{||} (Ions 1 resonance, n_{ant}=',num2str(n_ant),', max=100)'])
    else
        title(['k_{||} (Ions 1 resonance, n_{ant}=',num2str(n_ant),')'])
    end
end
if (k_lkp(i)==6)
    if(a)
        title(['k_{||} (Ions 2 resonance, n_{ant}=',num2str(n_ant),', max=100)'])
    else
        title(['k_{||} (Ions 2 resonance, n_{ant}=',num2str(n_ant),')'])
    end
end
if (k_lkp(i)==7)
    if(a)
        title(['k_{||} (based on phi, n_{ant}=',num2str(n_ant),', max=100)'])
    else
        title(['k_{||} (based on phi, n_{ant}=',num2str(n_ant),')'])
    end
end
if (k_lkp(i)==8)
    if(a)
        title(['k_{perp} (based on phi, n_{ant}=',num2str(n_ant),', max=100)'])
    else
        title(['k_{perp} (based on phi, n_{ant}=',num2str(n_ant),')'])
    end
end

end

end


for i = 5:4+sum(OutSelect(5:7))
    if(i*OutSelect(i)==5)
        temp_eps(:,:)=eps_nn_s(:,:);
    end
    if(i*OutSelect(i)==6)
        temp_eps(:,:)=eps_bn_s(:,:);
    end
    if(i*OutSelect(i)==7)
        temp_eps(:,:)=eps_par_s(:,:);
    end
    figure
    set(gcf,'renderer','zbuffer');
    %contour(x3_s,z3_s,abs(temp_eps),200);
    surf(x3_s,z3_s,abs(temp_eps),'LineStyle','none')
    colormap hsv
    axis equal;
    view(0,90);
    hold on;
    plot3(xb,zb,yb,'k')

    if (i==5)
        title('\epsilon_{nn}')
    end
    if (i==6)
        title('\epsilon_{bn}')
    end
    if (i==7)
        title('\epsilon_{||}')
    end


end

