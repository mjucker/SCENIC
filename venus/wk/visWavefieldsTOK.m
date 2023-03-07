clear all;
endian='b'

manual_min = -1.0e0;
manual_max =  1.0e0;
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

freq=2.03e8
%freq=2.642e8;
%freq=2.0575e8;
B0=2.9;
%mi1_SI=1.2367*1.67e-27;
%mi1_SI=1.617*1.67e-27;
mi1_SI=3*1.67e-27;
%mi1_SI=1.67e-27;
charge=2*1.6e-19;

Bres1 = mi1_SI*freq/(charge*B0)
%Bres1 = 1.2367*mi1_SI*freq/(1.6e-19*B0)


% --------------------------------------------------------------------------
% Read the equilibrium from fort.37.DAT

fif=fopen('fort.37.DAT.N','r',endian);

fread(fif,1,'uint');
v=fread(fif,30,'char');
fread(fif,1,'uint');

fread(fif,1,'uint');
v=fread(fif,4,'int');
fread(fif,1,'uint');

niT     = v(1)
njT     = v(2)
nkT     = v(3)
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
  coorRes1(indt,2)=0;
  coorRes1(indt,3)=(crz(ind,i)*(tb2(it*njT,i)-Bres1)-crz(it*njT,i)*(tb2(ind,i)-Bres1))/(tb2(it*njT,i)-tb2(ind,i));
end
if (ip~=1 && (tb2(ind,i)-Bres1)*(tb2(ind-1,i)-Bres1)<0)
  indt=indt+1;
  coorRes1(indt,1)=(crx(ind,i)*(tb2(ind-1,i)-Bres1)-crx(ind-1,i)*(tb2(ind,i)-Bres1))/(tb2(ind-1,i)-tb2(ind,i));
  coorRes1(indt,2)=0;
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
Ppla3D_rv=fread(fid,NelR*Njk,'double');
fread(fid,1,'uint');
fread(fid,1,'uint');
Ppla3D_iv=fread(fid,NelR*Njk,'double');
fread(fid,1,'uint');

fread(fid,1,'uint');
E3D_rv=fread(fid,3*NelR*Njk,'double');
fread(fid,1,'uint');
fread(fid,1,'uint');
E3D_iv=fread(fid,3*NelR*Njk,'double');
fread(fid,1,'uint');

fread(fid,1,'uint');
B3D_rv=fread(fid,3*NelR*Njk,'double');
fread(fid,1,'uint');
fread(fid,1,'uint');
B3D_iv=fread(fid,3*NelR*Njk,'double');
fread(fid,1,'uint');

fread(fid,1,'uint');
A3D_rv=fread(fid,4*NelR*Njk,'double');
fread(fid,1,'uint');
fread(fid,1,'uint');
A3D_iv=fread(fid,4*NelR*Njk,'double');
fread(fid,1,'uint');

fread(fid,1,'uint');
Ppla_rv=fread(fid,1,'double');
fread(fid,1,'uint');

fread(fid,1,'uint');
Ppla_iv=fread(fid,1,'double');
fread(fid,1,'uint');

fread(fid,1,'uint');
AP_nb=fread(fid,1,'uint');
fread(fid,1,'uint');

Sp_name = zeros(15,AP_nb);
for i = 1:AP_nb
strlen=fread(fid,1,'uint');
text_tmp=fread(fid,strlen,'char');
Sp_name(1:strlen,i)=text_tmp;
fread(fid,1,'uint');
end
lgd=transpose(Sp_name);

PplaSp_rv=zeros(AP_nb,NelR*Njk);
PplaSp_iv=zeros(AP_nb,NelR*Njk);
for i = 1:AP_nb
fread(fid,1,'uint');
PplaSp_rv(i,:)=fread(fid,NelR*Njk,'double');
fread(fid,1,'uint');
fread(fid,1,'uint');
PplaSp_iv(i,:)=fread(fid,NelR*Njk,'double');
fread(fid,1,'uint');
end


fclose(fid);

Ppla3D_v = Ppla3D_rv + j*Ppla3D_iv;  clear Ppla3D_rv Ppla3D_iv;
E3D_v = E3D_rv + j*E3D_iv;  clear E3D_rv E3D_iv;
B3D_v = B3D_rv + j*B3D_iv;  clear B3D_rv B3D_iv;
A3D_v = A3D_rv + j*A3D_iv;  clear A3D_rv A3D_iv;
PplaSp_v = PplaSp_rv + j*PplaSp_iv;  clear PplaSp_rv PplaSp_iv;

cR  = reshape(cR_v, NelR,njT,nkT);
cZ  = reshape(cZ_v, NelR,njT,nkT);
Ppla3D = reshape(Ppla3D_v,NelR,njT,nkT);
E3D = reshape(E3D_v,3,NelR,njT,nkT);
B3D = reshape(B3D_v,3,NelR,njT,nkT);
A3D = reshape(A3D_v,4,NelR,njT,nkT);
PplaSp = reshape(PplaSp_v,AP_nb,NelR,njT,nkT);
clear Ppla3D_v E3D_v B3D_v A3D_v cR_v cZ_v PplaSp_v;



% Close poloidally

njT = njT+1;
for eR = 1:NelR
  for ip = 1:nkT
    cR(eR,njT,ip)  =  cR(eR,1,ip);
    cZ(eR,njT,ip)  =  cZ(eR,1,ip);
    E3D(:,eR,njT,ip) = E3D(:,eR,1,ip);
    Ppla3D(eR,njT,ip) = Ppla3D(eR,1,ip);
    PplaSp(:,eR,njT,ip) = PplaSp(:,eR,1,ip);
  end 
end


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
E3_s  = zeros(NelR,njT);
P3_s  = zeros(NelR,njT);
PSp_s = zeros(AP_nb,NelR,njT);

for eR = 1:NelR

x3_s(eR,:)  = x3 (eR,:,1);
y3_s(eR,:)  = y3 (eR,:,1);
z3_s(eR,:)  = z3 (eR,:,1);
E3_s(eR,:) = E3D(1,eR,:,1) - complex(0,1)*E3D(2,eR,:,1);
%E3_s(eR,:) = E3D(2,eR,:,1);
P3_s(eR,:) = Ppla3D(eR,:,1);
PSp_s(:,eR,:) = PplaSp(:,eR,:,1);

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
%x3_s = [tmp1D' x3_s']';
%tmp1D(:) = sum(y3_s(1,:))/njT;
%y3_s = [tmp1D' y3_s']';
%tmp1D(:) = sum(z3_s(1,:))/njT;
%z3_s = [tmp1D' z3_s']';
%tmp1D(:) = sum(E3_s(1,:))/njT;
%E3_s = [tmp1D' E3_s']';
%tmp1D(:) = sum(P3_s(1,:))/njT;
%P3_s = [tmp1D' P3_s']';


%E3_s(NelR-10:NelR+1,:) = 0;

% Border
xb = zeros(1,njT);
yb = zeros(1,njT);
zb = zeros(1,njT);
xb(:) = x3_s(NelR,:);
zb(:) = z3_s(NelR,:);

MAX_E=max(max(max(real(E3_s(:,:,:)))))
MIN_E=min(min(min(real(E3_s(:,:,:)))))

figure
set(gcf,'renderer','zbuffer');
%subplot(2,3,1)
%mesh(x3_s,z3_s,real(E3_s)+imag(E3_s2));
contour(x3_s,z3_s,real(E3_s),200);
%contour(x3_s,z3_s,abs(E3_s),ctr_highE);
%axis([2.75 4.75 -0.8 0.8]) 
%axis([5 11 -3 3 min(min(min(real(E3_s(:,:,:))))) max(max(max(real(E3_s(:,:,:)))))]) 
%axis([5 11 -3 3 -0.1 0.1]) 
%shading interp;
axis equal;
view(0,90);
hold on;
plot3(xb,zb,yb,'k',coorRes1s(:,1),coorRes1s(:,3),coorRes1s(:,2),'k-');
%plot3(xb,zb,yb,'k',coorRes2u(:,1),coorRes2u(:,3),coorRes2u(:,2),'k-',coorRes2d(:,1),coorRes2d(:,3),coorRes2d(:,2),'k-');


figure
set(gcf,'renderer','zbuffer');
%subplot(2,3,4)
%contourf(x3_s,z3_s,abs(P3_s),200.,'LineStyle','none');
contourf(x3_s,z3_s,-imag(P3_s),200,'LineStyle','none');
%contourf(x3_s,z3_s,real(P3_s),ctr_high);
axis equal;
%axis([2.75 4.75 -0.8 0.8]) 
%axis([1.5 4.5 -2.1 2.1 min(min(min(real(P3_s(:,:,:))))) max(max(max(real(P3_s(:,:,:)))))]) 
%shading interp;
view(0,90);
colorbar;
hold on;
plot3(xb,zb,yb,'w',coorRes1s(:,1),coorRes1s(:,3),coorRes1s(:,2),'w-');
%plot3(xb,zb,yb,'k',coorRes2u(:,1),coorRes2u(:,3),coorRes2u(:,2),'k-',coorRes2d(:,1),coorRes2d(:,3),coorRes2d(:,2),'k-');

%figure
%mesh(x3_s,z3_s,imag(P3_s))

for i = 1:AP_nb
temp_PSp(:,:)=PSp_s(i,:,:);
figure
set(gcf,'renderer','zbuffer');
contourf(x3_s,z3_s,-imag(temp_PSp),200,'LineStyle','none');
axis equal;
view(0,90);
colorbar;
hold on;
plot3(xb,zb,yb,'w',coorRes1s(:,1),coorRes1s(:,3),coorRes1s(:,2),'w-');
title(setstr(lgd(i,:)))
end


