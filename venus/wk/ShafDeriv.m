function [Shafp,Shafpp]=ShafDeriv(r,Shaf,plotter);

x=linspace(r(1),max(r),1000);
ShafSp=spline(r,Shaf,x);
pShafSp=diff(ShafSp)./diff(x); %one value at the end is missing!
Shafp=spline(x(1:length(x)-1),pShafSp,r);

%second derivative
ShafpSp=spline(r,Shafp,x);
pShafpSp=diff(ShafpSp)./diff(x);
Shafpp=spline(x(1:length(x)-1),pShafpSp,r);


if plotter>0
Shafplen=length(Shafp);
Shafpplen=length(Shafpp);
figure;plot(r(1:Shafplen),Shafp);
title('\Delta \prime (spline interp)');grid;
figure;plot(r(1:Shafpplen),Shafpp);
title('\Delta \prime \prime (spline interp)');grid;
end