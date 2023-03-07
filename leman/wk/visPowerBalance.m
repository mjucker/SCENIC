clear all;


% --------------------------------
% calculation results

fid=fopen('PowersDiff.txt');
fseek(fid,0,-1);

Nmesh = fscanf(fid,'%d \n',1);
mdat  = fscanf(fid,'%f \n',[9,Nmesh]);
s_fe  = mdat(1,:);
Ppla  = (mdat(2,:) + j*mdat(3,:));
Pant  = -(mdat(4,:) + j*mdat(5,:));
Pdiff = mdat(6,:) + j*mdat(7,:);
Pdiv = mdat(8,:) + j*mdat(9,:); clear mdat;


NelR = fscanf(fid,'%d \n',1);
mdat  = fscanf(fid,'%f \n',[3,NelR]);
s_hp  = mdat(1,:);
Ppoynt  = -(mdat(2,:) + j*mdat(3,:));; clear mdat;

NelR = fscanf(fid,'%d \n',1);
mdat  = fscanf(fid,'%f \n',[3,NelR]);
s_hp  = mdat(1,:);
PAdivA = (mdat(2,:) + j*mdat(3,:));; clear mdat;

Ntot = fscanf(fid,'%d \n',1);
mdat  = fscanf(fid,'%f \n',[2,Ntot]);
s_tot  = mdat(1,:);
divAint  = mdat(2,:); clear mdat;

%Read power deposition of all species

nb_sp = fscanf(fid,'%d \n',1);

lgd = zeros(nb_sp,15);
lgd(1,1:5) = 'Total';
for i = 1:nb_sp
num_tmp=fscanf(fid,'%d \n',1);

if(num_tmp==1)
lgd(i+1,1:9) = 'Electrons';
end
if(num_tmp==2)
lgd(i+1,1:13) = 'First Species';
end
if(num_tmp==3)
lgd(i+1,1:14) = 'Second Species';
end
if(num_tmp==4)
lgd(i+1,1:9) = 'Fast Ions';
end

end

for i = 1:nb_sp
PplaSp(i,:,:) = fscanf(fid,'%f \n',[2,Nmesh]);
end

fclose(fid);


% Comparing two functions: Ppla+Pant (integer mesh points) and Ppoynt (half-mesh points)
% Should be:  Ppla_hp + Pant_hp + Ppoynt = 0
% Compare on the s_int grid

Nint = 100000;
s_int = eps:(1.0-2*eps)/(Nint-1):1.0-eps;

% Functions to compare
F1_hp = -Ppoynt;
F2_fe = Ppla + Pant;


% Interpolate F1,F2 to s_int grid
% Boundary conditions: F1_fe(0)=F1_fe(1)=0;  F2_fe(0)=F2_fe(1)=0

F1_hp00  = [0,F1_hp,0];
s_hp00   = [0,s_hp,1];
F1_int   = spline(s_hp00,F1_hp00,s_int);
F2_int   = spline(s_fe,F2_fe,s_int);

Error_func_re = abs(real((F1_int - F2_int)));
Error_func_im = abs(imag((F1_int - F2_int)));
Error_func_abs = abs((F1_int - F2_int));
Error_integral_re = sum(Error_func_re)/Nint;
Error_integral_im = sum(Error_func_im)/Nint;
Error_integral_abs = sum(Error_func_abs)/Nint;

MaxError_func_re = max(abs(real(F1_int - F2_int)));
MaxError_func_im = max(abs(imag(F1_int - F2_int)));
MaxError_func_abs = max(abs(abs(F1_int - F2_int)));




disp(sprintf('Radial grid: %5d',Nmesh));
disp(sprintf('Global energy balance, rel.'));
disp(sprintf('  Real part:      %0.5e',abs(real(Ppla(Nmesh)+Pant(Nmesh)))/abs(Pant(Nmesh))));
disp(sprintf('  Imaginary part: %0.5e',abs(imag(Ppla(Nmesh)+Pant(Nmesh)))/abs(Pant(Nmesh))));
disp(sprintf('  Abs:            %0.5e',abs(Ppla(Nmesh)+Pant(Nmesh))/abs(Pant(Nmesh))));
disp(sprintf('Relative error in the local energy balance, integrated'));
disp(sprintf('  Real part:      %0.5e',Error_integral_re/abs(Pant(Nmesh))));
disp(sprintf('  Imaginary part: %0.5e',Error_integral_im/abs(Pant(Nmesh))));
disp(sprintf('  Abs.:           %0.5e',Error_integral_abs/abs(Pant(Nmesh))));
disp(sprintf('Relative error in the local energy balance, maximum'));
disp(sprintf('  Real part:      %0.5e',MaxError_func_re/abs(Pant(Nmesh))));
disp(sprintf('  Imaginary part: %0.5e',MaxError_func_im/abs(Pant(Nmesh))));
disp(sprintf('  Abs.:           %0.5e',MaxError_func_abs/abs(Pant(Nmesh))));
disp(sprintf('Avrg. Pant value used for the error estimate'));
disp(sprintf('  Real part:      %0.5e',real(Pant(Nmesh))));
disp(sprintf('  Imaginary part: %0.5e',imag(Pant(Nmesh))));



fid = fopen('convergence.txt','a');
%fprintf(fid,'%5d \n',NelR);
%fprintf(fid,'%0.12e \n',real(Ppla(Nmesh)+Pant(Nmesh))); % Global power balance, real
%fprintf(fid,'%0.12e \n',imag(Ppla(Nmesh)+Pant(Nmesh))); % Global power balance, imag

fprintf(fid,'%s \n','Local error, integral'); % Local power balance
fprintf(fid,'%22.12e %22.12e \n',Error_integral_re, Error_integral_im);
fprintf(fid,'%s \n','Local error, maximum'); % Local power balance
fprintf(fid,'%22.12e %22.12e \n',MaxError_func_re, MaxError_func_im);

%fprintf(fid,'%0.12e \n',real(Ppla(Nmesh))); % Absorbed power, real
%fprintf(fid,'%0.12e \n',imag(Ppla(Nmesh))); % Absorbed power, imag
fclose(fid);



figure;
subplot(4,1,1);
plot(s_fe,real(Ppla),'-o', s_fe,real(-Pant),'-x'); grid on;
ylabel('Re[Ppla,Pant]','Rotation',[90]); xlabel('s');

subplot(4,1,2);
plot(s_int,real(F2_int),'-', s_hp,real(F1_hp),'-'); grid on;
ylabel('Re[Ppla+Pant,-Ppoynt]','Rotation',[90]); xlabel('s');

subplot(4,1,3);
plot(s_int,real(F1_int-F2_int),'-o'); grid on;
ylabel('Re[Ppla+Pant-Ppoynt]','Rotation',[90]); xlabel('s');

subplot(4,1,4);
plot(s_fe,real(Pdiv),'-o',s_hp,real(PAdivA),'-x'); grid on;
ylabel('Re[Pdiv, PAdivA]','Rotation',[90]); xlabel('s');


figure;
subplot(4,1,1);
plot(s_fe,imag(Ppla),'-o', s_fe,imag(-Pant),'-x'); grid on;
ylabel('Im[Ppla,Pant]','Rotation',[90]); xlabel('s');

subplot(4,1,2);
plot(s_int,imag(F2_int),'-', s_hp,imag(F1_hp),'-'); grid on;
ylabel('Im[Ppla+Pant,-Ppoynt]','Rotation',[90]); xlabel('s');

subplot(4,1,3);
plot(s_int,imag(F1_int-F2_int),'-o'); grid on;
ylabel('Im[Ppla+Pant-Ppoynt]','Rotation',[90]); xlabel('s');

subplot(4,1,4);
plot(s_fe,imag(Pdiv),'-o',s_hp,imag(PAdivA),'-x'); grid on;
ylabel('Im[Pdiv, PAdivA]','Rotation',[90]); xlabel('s');


figure;
plot(s_fe,imag(Ppla),'k-'); hold on; grid on;

cc(1,:) = [0 0 1];
cc(3,:) = [1 0 0];
cc(2,:) = [0.7 0 1];
cc(4,:) = [0 0.7 0.7];


for i = 1:nb_sp
Ppla_plot(1,:) = PplaSp(i,2,:);
h=plot(s_fe,Ppla_plot,'-');
ylabel('Im[Ppla]','Rotation',[90]); xlabel('s');
set(h,'Color',cc(i,:));
end
legend(setstr(lgd(:,:)));



break;




clear tmp;  tmp = zeros(1,Nmesh);













Ppla_hp = (interp1(s_fe,Ppla, s_hp));
Pant_hp = (interp1(s_fe,Pant, s_hp));
P_sum_hp = Ppla_hp + Pant_hp + Ppoynt;


Pdiv_hp = (interp1(s_fe,Pdiv, s_hp));




% ---  Compare ---

figure;
subplot(3,1,1);
plot(s_fe,real(Pdiv),'-o'); grid on;
ylabel('Re[Pdiv]','Rotation',[90]); xlabel('s');

subplot(3,1,2);
plot(s_hp,real(PAdivA),'-o'); grid on;
ylabel('Re[AdivA]','Rotation',[90]); xlabel('s');

subplot(3,1,3);
plot(s_hp,real(PAdivA + Pdiv_hp),'-o'); grid on;
ylabel('Re[AdivA - Pdiv]','Rotation',[90]); xlabel('s');


figure;
subplot(3,1,1);
plot(s_fe,real(Ppla),'-o',s_fe,real(-Pant),'-o'); grid on;
ylabel('Re[Ppla,Pant]','Rotation',[90]); xlabel('s');

subplot(3,1,2);
title('Approximate balance: sum=0');
plot(s_fe,real(Ppla+Pant),'-o', s_hp,real(-Ppoynt),'-o'); grid on;
ylabel('Re[Ppla+Pant]','Rotation',[90]); xlabel('s');

subplot(3,1,3);
title('Approximate balance: sum=0');
plot(s_fe,real(Ppla+Pant+Pdiv),'-o', s_hp,real(-Ppoynt-PAdivA),'-o'); grid on;
ylabel('Re[Ppla+Pant]','Rotation',[90]); xlabel('s');




break;






figure;
subplot(3,1,1);
plot(s_fe,real(Pdiv),'-o'); grid on;
ylabel('Re[Pdiv]','Rotation',[90]); xlabel('s');

subplot(3,1,2);
plot(s_hp,real(PAdivA),'-o'); grid on;
ylabel('Re[AdivA]','Rotation',[90]); xlabel('s');

subplot(3,1,3);
plot(s_hp,real(PAdivA + Pdiv_hp),'-o'); grid on;
ylabel('Re[AdivA - Pdiv]','Rotation',[90]); xlabel('s');


figure;
subplot(4,1,1);
plot(s_fe,real(Ppla),'-o',s_fe,real(-Pant),'-o'); grid on;
ylabel('Re[Ppla,Pant]','Rotation',[90]); xlabel('s');

subplot(4,1,2);
title('Approximate balance: sum=0');
plot(s_fe,real(Ppla+Pant),'-o'); grid on;
ylabel('Re[Ppla+Pant]','Rotation',[90]); xlabel('s');

subplot(4,1,3);
title('Approximate balance: sum=0');
plot(s_hp,real(Ppoynt),'-o'); grid on;
ylabel('Re[Ppoynt]','Rotation',[90]); xlabel('s');

subplot(4,1,4);
title('Approximate balance: sum=0');
plot(s_hp,real(P_sum_hp),'-o'); grid on;
ylabel('Re[Ppla+Pant+Ppoynt]','Rotation',[90]); xlabel('s');


break;


figure;
subplot(4,1,1);
plot(s_fe,imag(Ppla),'-o',s_fe,imag(-Pant),'-o'); grid on;
ylabel('Im[Ppla,Pant]','Rotation',[90]); xlabel('s');

subplot(4,1,2);
title('Approximate balance: sum=0');
plot(s_fe,imag(Ppla+Pant),'-o'); grid on;
ylabel('Im[Ppla+Pant]','Rotation',[90]); xlabel('s');

subplot(4,1,3);
title('Approximate balance: sum=0');
plot(s_hp,imag(Ppoynt),'-o'); grid on;
ylabel('Im[Ppoynt]','Rotation',[90]); xlabel('s');

subplot(4,1,4);
title('Approximate balance: sum=0');
plot(s_hp,imag(P_sum_hp),'-o'); grid on;
ylabel('Im[Ppla+Pant+Ppoynt]','Rotation',[90]); xlabel('s');





figure;
plot(s_tot,divAint,'-o'); grid on;
ylabel('divA, integral on \theta,\phi','Rotation',[90]); xlabel('s');







break

figure;
subplot(2,1,1);
plot(s_fe,imag(Ppla),'-o',s_fe,imag(-Pant),'-o'); grid on;
ylabel('Im[Ppla,Pant]','Rotation',[90]); xlabel('s');

subplot(2,1,2);
title('Approximate balance: sum=0');
plot(s_hp,imag(P_sum_hp),'-o'); grid on;
ylabel('Im[sum]','Rotation',[90]); xlabel('s');




break


figure;
subplot(2,1,1);
plot(s_fe,real(Ppla),'-o',s_fe,real(Pant),'-o'); grid on;
ylabel('Re[Ppla,Pant]','Rotation',[90]); xlabel('s');

subplot(2,1,2);
title('Approximate balance: sum=0');
plot(s_fe,real(Pdiff),'-o'); grid on;
ylabel('Re[sum]','Rotation',[90]); xlabel('s');


figure;
subplot(2,1,1);
plot(s_fe,imag(Ppla),s_fe,imag(Pant)); grid on;
ylabel('Im[Ppla,Pant]','Rotation',[90]); xlabel('s');

subplot(2,1,2);
plot(s_fe,imag(Pdiff)); grid on;
ylabel('Im[sum]','Rotation',[90]); xlabel('s');


