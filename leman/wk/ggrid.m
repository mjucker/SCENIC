%  Radial grid generation


clear all;

NelR = 150;

Ni = 100000;

x_f = 1.e-3:1./Ni:1.;
%x_f = 1.e-6:1./Ni:1.;


weight = zeros(1,length(x_f));
weight(:) = 30.;

%for ii = 1:length(x_f)
%  if (x_f(ii)>0.65) 
%    weight(ii) = weight(ii) + 220;
%  end
%end

%weight(:) = weight(:) + 200*exp(-(x_f(:)-0.40).^2/0.012);

%weight(:) = weight(:) + 100*exp(-(x_f(:)-0.5625).^2/0.01);
%weight(:) = weight(:) + 400*exp(-(x_f(:)-0.30).^2/0.002);
%weight(:) = weight(:) + 200*exp(-(x_f(:)-0.0484).^2/0.01);

%for ii = 1:length(x_f)
%  if ((sqrt(x_f(ii))>0.69)&(sqrt(x_f(ii))<0.81)) 
%    weight(ii) = weight(ii) + 30;
%  end
%end

%weight(:) = weight(:) + (10.^(10-9*x_f(:)));



%weight(:) = weight(:) + 3./(x_f(:));



%weight(:) = weight(:) + 3./(exp(x_f(:)));

%weight(:) = weight(:) + (3.^(1.0e-1./(x_f(:))));

accw = cumtrapz(x_f,weight);
accw(:) = accw(:)/accw(length(accw));

y_n = 0:1./NelR:1;
x_n = interp1(accw,x_f,y_n);

% rescaling (x(1) -> 0)
x_n(:) = (x_n(:) - x_n(1))/(x_n(length(x_n)) - x_n(1));

x_n(:) = x_n(:).^2

fid=fopen('grid.dat', 'W');
fprintf(fid,'%d \n',NelR);
fprintf(fid,'%f \n',x_n);
fclose(fid);



%figure;
%plot(x_f, weight); grid on;% axis([0 1 0 4]);


figure;
plot(x_n, x_n, 'o'); grid on;% axis([0 1 0 1]);
disp(  sprintf('%d / %d',length(find(x_n>0.65)),NelR) )

length(find(x_n>0.75))/96 * NelR
