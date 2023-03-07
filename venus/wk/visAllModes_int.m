function [x_3d,y_3d,temp_3d]=visAllModes_int(es,iota,re_im)
%clear all;

% For Hermite cubics visualization:
PPC = 10;  % Internal points per radial cell


% --------------------------------
% calculation results


fid=fopen('RawA.dat','r','b');



fread(fid,1,'uint');
v=fread(fid,8,'double');
fread(fid,1,'uint');
Om_SI = v(4);
clear v;

fread(fid,1,'uint');
v=fread(fid,3,'int');
fread(fid,1,'uint');
m_ant = v(1);
n_ant = v(2);
nperiod = v(3);
clear v;

fread(fid,1,'uint');
v=fread(fid,2,'int');
fread(fid,1,'uint');
Nmesh = v(1);
mnAFTot = v(2);
clear v;

fread(fid,1,'uint');
mLkp=fread(fid,mnAFTot,'int');
fread(fid,1,'uint');
clear v;

fread(fid,1,'uint');
nLkp=fread(fid,mnAFTot,'int');
fread(fid,1,'uint');
clear v;

fread(fid,1,'uint');
s_A=fread(fid,Nmesh,'double');
fread(fid,1,'uint');


fread(fid,1,'uint');
XXc_rv=fread(fid,Nmesh*2*mnAFTot*4,'double');
fread(fid,1,'uint');

fread(fid,1,'uint');
XXc_iv=fread(fid,Nmesh*2*mnAFTot*4,'double');
fread(fid,1,'uint');


fclose(fid);

XXc_v = XXc_rv + j*XXc_iv;  clear XXc_rv XXc_iv;
XXc = reshape(XXc_v,Nmesh*2,mnAFTot,4);
clear XXc_v;


NRfine = (Nmesh-1) * (PPC+1);
A = zeros(4,mnAFTot,NRfine);
pf = -1;

for eR = 1:Nmesh-1
for ie = 0:PPC

  i_fine = (eR-1)*(PPC+1)+ie;

  if (i_fine > 0)

  s_fine(i_fine) = s_A(eR) + (s_A(eR+1)-s_A(eR))*ie/(PPC+1);
  ksi = 2.0*ie/(PPC+1) - 1;

  Psi_al1 = 0.25  * ( ksi*ksi*ksi - 3*ksi + 2);
  Psi_al2 = 0.125 * ( ksi*ksi*ksi - ksi*ksi - ksi + 1)*(s_A(eR+1)-s_A(eR));
  Psi_al3 = 0.25  * (-ksi*ksi*ksi + 3*ksi + 2);
  Psi_al4 = 0.125 * ( ksi*ksi*ksi + ksi*ksi - ksi - 1)*(s_A(eR+1)-s_A(eR));

  sfr(1) = s_fine(i_fine)^pf;
  sfr(2) = s_fine(i_fine)^pf;
  sfr(3) = 1;
  sfr(4) = 1;
  for k = 1:4
  for mn = 1:mnAFTot
    A(k,mn,i_fine) = (XXc(2*(eR-1)+1,mn,k)*Psi_al1 + XXc(2*(eR-1)+2,mn,k)*Psi_al2 + ...
                      XXc(2*(eR-1)+3,mn,k)*Psi_al3 + XXc(2*(eR-1)+4,mn,k)*Psi_al4) * sfr(k);
  end;
  end;

  end;

end;
end;

% Last point:
i_fine = NRfine;
s_fine(i_fine) = 1;
for k = 1:4
for mn = 1:mnAFTot
A(k,mn,i_fine) = XXc(2*(eR-1)+3,mn,k)*Psi_al3;
end;
end;



temp = zeros(1,NRfine); 


if (re_im == 're')
  Ac = real(A);
else
%  Ac = imag(A);
Ac = abs(A);
end

clear A;
A = Ac;


% Normalise
%for k = 1:4
%  A(k,:,:,:) = A(k,:,:,:) / max(max(max(abs(A(k,:,:,:)))));
%end
GlobMaxA3 = max(max(abs(A(3,:,:))));


figure
title('Potentials');

temp = zeros(1,NRfine);
plot_sA = s_fine;
%plot_sA = 18*sqrt(s_fine);

for mn = 1:mnAFTot

  cc = [rand rand rand];

  sm = strcat(int2str(mLkp(mn)));
  sn = strcat(int2str(nLkp(mn)));
  slabel = strcat(sm,',',sn);



  temp(:) = A(1,mn,:);  subplot(2,2,1); hold on;
  h=plot(plot_sA(:),(real(temp(:,:)))); grid on;
  set(h,'LineWidth',[2]);  ylabel('A_n','Rotation',[0]);  xlabel('s'); 
  set(h,'Color',cc);
  [my,itext] = max(abs(real(temp)));  xtext = plot_sA(itext);
  ytext = real(temp(itext));
  subplot(2,2,1); h=text(xtext,ytext,slabel); set(h,'Color',cc);

  
  temp(:) = A(2,mn,:);  subplot(2,2,2); hold on;
  h=plot(plot_sA(:),(real(temp(:,:)))); grid on;
  set(h,'LineWidth',[2]);  ylabel('A_b','Rotation',[0]);  xlabel('s'); 
  set(h,'Color',cc);
  [my,itext] = max(abs(real(temp)));  xtext = plot_sA(itext);
  ytext = real(temp(itext));
  subplot(2,2,2); h=text(xtext,ytext,slabel); set(h,'Color',cc);


  temp(:) = A(3,mn,:);  subplot(2,2,3); hold on;
  h=plot(plot_sA(:),(real(temp(:,:)))); grid on;
  set(h,'LineWidth',[2]);  ylabel('A_{||}','Rotation',[0]);  xlabel('s'); 
  set(h,'Color',cc);
  [my,itext] = max(abs(real(temp)));  xtext = plot_sA(itext);
  ytext = real(temp(itext));
  subplot(2,2,3); h=text(xtext,ytext,slabel); set(h,'Color',cc);


  temp(:) = A(4,mn,:);  subplot(2,2,4); hold on;
  h=plot(plot_sA(:),(real(temp(:,:)))); grid on;
  set(h,'LineWidth',[2]);  ylabel('\phi','Rotation',[0]);  xlabel('s'); 
  set(h,'Color',cc);
  [my,itext] = max(abs(real(temp)));  xtext = plot_sA(itext);
  ytext = real(temp(itext));
  subplot(2,2,4); h=text(xtext,ytext,slabel); set(h,'Color',cc);

end;



%%subplot(3,2,5);
%%plot(s_T,iota_T,'--',s_T,1./iota_T); grid on;
%
%%ylabel('i,q','Rotation',[0]); xlabel('s');
%%legend('i','q');
%
%h=subplot(3,2,6);
%axis off;
%% Parameters:
%str_param = strvcat('Equilibrium parameters:', ...
%	     sprintf('n0           = %10.4e ,m^{-3}', n0_SI), ...
%	     sprintf('B0           = %10.4e ,T',      B0), ...
%	     sprintf('Mi/me        = %10.4e',         Mie), ...
%	     sprintf('C/Ca         = %10.4e',         3.e8/cA_SI), ...
%	     sprintf('Om(ci)*a/Ca  = %10.4e',         OmCi_SI/cA_SI), ...
%	     sprintf('Omega        = %10.4e ,1/rad',  Om_SI), ...
%	     sprintf('Omega.Imag.  = %10.4e',         om_im), ...
%	     sprintf('Omega/OmCi   = %10.4e',         om), ...
%	     sprintf('Major radius = %10.4e ,m',      MajRadius), ...
%	     sprintf('Ant. tor. m. = %10.4e,',        n_ant), ...
%	     sprintf('Ant. por. m. = %10.4e,',        m_ant) ...
%	     );
%
%scale_axis = axis;
%xtext = (scale_axis(1)+scale_axis(2))/2;
%ytext = (scale_axis(4)+scale_axis(3))/2;
%
%t_handle=text(0,ytext,str_param);
%set(t_handle,'FontName','Courier');



%Ppla = Ppla + PdivA;

%PdivA

for mn = 1:mnAFTot
A_mn(mn) = mLkp(1)+mn-1;
end

temp_surf = zeros(mnAFTot);
temp_surf = abs(A(3,:,200));

for mn = 1:mnAFTot
for s_index = 1:Nmesh-1
iotasp=spline(es,iota,plot_sA((PPC+1)*(s_index-1)+PPC/2));
x_3d(mn,s_index) = (mLkp(1)+mn-1)*iotasp;
y_3d(mn,s_index) = sqrt(plot_sA((PPC+1)*(s_index-1)+PPC/2));

temp_3d(mn,s_index) = abs(A(1,mn,(PPC+1)*(s_index-1)+PPC/2));
end
end

size(x_3d)
size(y_3d)
size(temp_3d)

figure
bar(A_mn,temp_surf);
figure
%mesh(x_3d,y_3d,temp_3d);
pcolor(x_3d,y_3d,temp_3d);
figure
contour(x_3d,y_3d,temp_3d,200);

% A parallel
figure
title(sprintf('%s,  Omega=%9.4e',re_im,Om_SI));
temp = zeros(1,NRfine);
hold on;


for mn = 1:mnAFTot

  cc = [rand rand rand];

  sm = strcat(int2str(mLkp(mn)));
  sn = strcat(int2str(nLkp(mn)));
  slabel = strcat(sm,',',sn);


  temp(:) = A(3,mn,:); 
  h=plot(plot_sA(:),(real(temp(:)))); grid on;
  set(h,'LineWidth',[2]);
  set(h,'Color',cc);
  [my,itext] = max(abs(real(temp)));  xtext = plot_sA(itext);

  maxModeAmp = max(max(abs(A(3,mn,:))));
%  if ((mLkp(m) >= -7) & (mLkp(m) <= -1))
  if (maxModeAmp >= 0.10*GlobMaxA3)
  ytext = real(temp(itext));
  h=text(xtext,ytext,slabel); set(h,'Color',cc);
  set(h,'Fontsize',[18]);  
  end;

  ha=gca;
  set(ha,'Fontsize',[14]);  
%  hy=ylabel('A_{||}','Rotation',[0]);
%  set(hy,'Fontsize',[14]);  
%  hx=xlabel('s');
%  set(hx,'Fontsize',[14]);  
end;
