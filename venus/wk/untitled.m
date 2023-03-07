f=load('fort.49');t=f(:,1);s=f(:,2);V=f(:,3);dens=f(:,4);ppar=f(:,5);pperp=f(:,6);
for i=1:1:ndiag
    sloc=s((i-1)*96+1:i*(96));
    p=polyfit(sloc,dens((i-1)*96+1:i*(96)).*V((i-1)*96+1:i*(96)),4);
    pdens=polyval(p,s((i-1)*96+1:i*(96)));
    plot(a00*sqrt(sloc),pdens);
    pause(.1);
end
