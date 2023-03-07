function [r,R0,Shaf,Sm,kappa,delta]=FourMath(plotter,ns)
%PLOTS SHAFRANOV SHIFT, ELONGATION AND TRIANGULARITY
% GENERATED FROM VMEC AS DEFINED IN THE ANALYTICAL WORK
fid8=fopen('fort.49','rt');
[b,ncount]=fscanf(fid8,'%g %g %g %g',[4,inf]);
fclose(fid8);
%ns=129;
nmodes=ncount/4/ns;
ntotal=ncount/4;
%b
r=zeros(1,ns);
R0=b(3,1);

Shaf=zeros(1,ns);
Sm=zeros(nmodes,ns);
s=0;re=0;
m=zeros(1,nmodes);
for i=1:ntotal;
 if b(1,i)==0.0
  s=s+1;
  Shaf(s)=R0-b(3,i);
 elseif b(1,i)==1.0
  re=re+1;
  r(re)=0.5*(b(3,i)+b(4,i));
  Sm(2,re)=0.5*(b(3,i)-b(4,i));
 elseif b(1,i)>nmodes-1
  'some really bad error with the number of modes!'
  return
 elseif b(1,i)>2
  mode=b(1,i);
  m(mode)=m(mode)+1;
  for bi=1:nmodes
      if b(1,i-bi)==b(1,i)-1
          b1=i-bi;
          break
      end
  end
  Sm(mode,m(mode))=b(3,b1);
 end
end

r(1)=1e-14;
[SmNew]=TruncExp(Shaf,Sm,r);
Sm=SmNew;
kappa=(r-Sm(2,:))./(r+Sm(2,:));
delta=4.*Sm(3,:)./r;

r(1)=0;
Sm(2,1)=0;
%r=[0,r];
%Shaf=[0,Shaf];
%sz=size(Sm);
%Sm=[zeros(sz(1),1),Sm];

if plotter==1

figure;
%plot(r,Shaf./r);
plot(r,Shaf);
axis([0,max(r),min(Shaf),max(Shaf)]);
xlabel('r');
%ylabel('Delta(r)/r')
ylabel('\Delta(r)')
title('Shafranov shift divided by r');
title('Shafranov shift')
grid;

figure;
%plot(r,Sm(2,:)./r);
plot(r,Sm(2,:));
axis([0,max(r),min(Sm(2,:)),max(Sm(2,:))]);
%axis([1,2,3,4]);
xlabel('r');
ylabel('S2(r)/r')
ylabel('S2(r)')
title('Elongation divided by r');
title('Elongation');
grid;

figure;
%plot(r,Sm(3,:)./r);
plot(r,Sm(3,:));
axis([0,max(r),min(Sm(3,:)),max(Sm(3,:))]);
%axis([1,2,3,4]);
xlabel('r');
%ylabel('S3(r)/r')
ylabel('S3(r)')
title('Triangularity divided by r');
title('Triangularity')
grid;

end


