clear fif fif2 fif3 in vperpm vparm vparbox vperpbox B A v w As Ap distmatrix
clear f1 f2 sum_f b

fif=fopen('fort.39','r',endian);
fif3=fopen('fort.69','r',endian);

fread(fif,1,'uint');
in=fread(fif,7,'int');
fread(fif,1,'double');
nbinsvpar=in(1),nbinsvperp=in(2),nbinspol=in(3),nbinsrad=in(4),ndiag=in(5);tfin=in(6);
nlines=in(1)*in(2)*in(3)*in(4);
% tst=exist('intplot');
% if(tst)
%     intplot=max(intplot,ndiag/100); %at most 100 time step plots
% else
%     intplot=max(1,ndiag/100);
% end
if(indpol==0)indpol=(nbinspol+1)/2;end
%start=nbinsvperp*nbinsvpar*((indrad-1)*nbinspol+indpol-1)+1;

fread(fif3,1,'uint');
in3=fread(fif3,2,'int');
fread(fif3,1,'uint');
fread(fif3,1,'uint');
in4=fread(fif3,3,'double');
fread(fif3,1,'uint');
nbinse=in3(1);
nall=in4(2);
ntot=in4(3);

fif2=fopen('fort.59','r',endian);

if(plotsep<1)
    f1=rand;
    f2=rand;
    figure('Tag',num2str(f1));
    set(gcf,'Position',[5 541 1000 500]);
    figure('Tag',num2str(f2));
    set(gcf,'Position',[361    36   916   319]);
end


A=zeros(nbinsvperp,nbinsvpar);
B=zeros(nbinsvperp,nbinsvpar);

if(~exist('ndiagnostics'))
    ndiagl=ndiag;
else
    ndiagl=ndiagnostics+2;
end
if(plotsum==0),ndiagl=ndiag-1;end;
for i=1:ndiagl

%%%%%%%%% read fort.59 %%%%%%%%%%%
    
    fread(fif2,1,'uint');
    t=fread(fif2,1,'double');
    fread(fif2,1,'uint');

    fread(fif2,1,'uint');   
    ss=fread(fif2,nbinsrad+1,'double');
    fread(fif2,1,'uint');

    fread(fif2,1,'uint');
    dens=fread(fif2,nbinsrad+1,'double');
    fread(fif2,1,'uint');

    fread(fif2,1,'uint');
    denserr=fread(fif2,nbinsrad,'double');
    fread(fif2,1,'uint');
    denserr=[0;denserr];        
    
    fread(fif2,1,'uint');
    densth=fread(fif2,nbinsrad+1,'double');
    fread(fif2,1,'uint');
    
    fread(fif2,1,'uint');
    ppar=fread(fif2,nbinsrad+1,'double');
    fread(fif2,1,'uint');

    fread(fif2,1,'uint');
    pparerr=fread(fif2,nbinsrad,'double');
    fread(fif2,1,'uint');

    fread(fif2,1,'uint');
    pperp=fread(fif2,nbinsrad+1,'double');
    fread(fif2,1,'uint');

    fread(fif2,1,'uint');
    pperperr=fread(fif2,nbinsrad,'double');
    fread(fif2,1,'uint');
    
    fread(fif2,1,'uint');
    currt=fread(fif2,nbinsrad+1,'double');
    fread(fif2,1,'uint');
    
    fread(fif2,1,'uint');
    currp=fread(fif2,nbinsrad+1,'double');
    fread(fif2,1,'uint');
    
%%%%%% read fort.39 %%%%%%%%%%%%%%%%%    
    
    fread(fif,1,'uint');
    vparm=fread(fif,nbinsvpar,'double');
    fread(fif,1,'uint');
    
    fread(fif,1,'uint');
    vperpm=fread(fif,nbinsvperp,'double');
    fread(fif,1,'uint');   
    
    fread(fif3,1,'uint');
    Ebin=fread(fif3,nbinse,'double');
    fread(fif3,1,'uint');
    
%    for rx=1:nbinsrad
%        for p=1:nbinspol
            fread(fif,1,'uint');    
            distmatrix=fread(fif,[nbinsvpar,nbinsvperp],'double');
            fread(fif,1,'uint');   
 %           if(rx==indrad && p==indpol)
 %               distmatrix=temp;
 %           end
 %       end
 %   end
    fread(fif3,1,'uint');    
    distE=fread(fif3,nbinse,'double');
    fread(fif3,1,'uint');   
    
     fread(fif3,1,'uint');    
     distEth=fread(fif3,nbinse,'double');
     fread(fif3,1,'uint');   
%    distEth=zeros(size(distE));

    for n=1:nbinsvpar
        for m=1:nbinsvperp
            distmatrix(n,m)=distmatrix(n,m)/vperpm(m);
        end
    end
    
    if(i==1)
%        r=round(100*sqrt((indrad-1)/nbinsrad+1/(2*nbinsrad))*a00)/100;
%        r=round(100*a00*sqrt(ss(indrad)))/100;
        r=round(100*sqrt(ss(indrad+1)))/100;
        indpl=indpol;
        if(indpol==0),indpl=ceil(nbinspol/2);end
        theta=round(10*((indpl-1)*360/nbinspol+180/nbinspol))/10-180;
       
%        theta=round(10*(indpol*360/nbinspol))/10;
        theta=mod(theta,360);
    end
%    vperpm=f2(2,start:nbinsvpar:start+nbinsvpar*nbinsvperp-1);
%    vparm=f2(1,start:start+nbinsvpar-1);
    vparmin=min(vparm)-(vparm(1)-vparm(2))/2;
    %vparmin=min(vparm)+(vparm(1)-vparm(2))/2;
    vparmax=max(vparm)-(vparm(end)-vparm(end-1))/2;
    %vparmax=max(vparm)+(vparm(end)-vparm(end-1))/2;
    vperpmin=min(vperpm)+(vperpm(1)-vperpm(2))/2;
    vperpmax=max(vperpm)+(vperpm(end)-vperpm(end-1))/2;
    ax=[vparmin vparmax vperpmin vperpmax];
    axmax=max(vparmax,vperpmax);
    ax=[-axmax axmax vperpmin axmax];
    vparbox(1)=vparmin;
    vparbox(2:nbinsvpar)=-vparm(1:end-1)-diff(vparm)/2;
    %vparbox(2:nbinsvpar)=+vparm(1:end-1)+diff(vparm)/2;
    vparbox(nbinsvpar+1)=vparmax;
    l=ceil(log10(vparbox(end)));
   clear tickst
     tickst(1,:)=[num2str(round(vparbox(1)/(10^(l-3)))*10^(-2)),'e0',num2str(l-1)];
     tickst(ceil(length(vparbox)/2),1)='0';
     tickst(length(vparbox),:)=['+',num2str(round(vparbox(end)/(10^(l-3)))*10^(-2)),'e0',num2str(l-1)];
    vperpbox(1)=vperpmin;
    vperpbox(2:nbinsvperp)=vperpm(1:end-1)+diff(vperpm)/2;
    vperpbox(nbinsvperp+1)=vperpmax;
    Bdet=zeros(nbinsvperp,nbinsvpar);
    k=-1;v=0;w=0;
%    start=nbinsvperp*nbinsvpar*(indpol-1+nbinspol*(indrad-1+nbinsrad*(i-1)))+1;
%    for n=1:nbinsvperp
%        for m=1:nbinsvpar
%            k=k+1;
%            Bdet(n,m)=Bdet(n,m)+distmatrix(m,n)/ndiag;
    Bdet=Bdet+distmatrix'/ndiag;
%        end
%    end
    A=A+Bdet;
    if (mod(i,intplot)==0 || i==ndiagl)
        if(plotsep<1)
         figure(findobj('Tag',num2str(f1)));
                
         suptitle(['t=',num2str(t),'s']);
         if(i==ndiag)
             suptitle(['<f>_t']);
         end
        end
         
         B=A-B;
         if(plotsep<1);subplot(3,4,[1 2 5 6]);end
         v=unique(B);
         step=max(1,floor(length(v)/10));
 %        step=1;
         w=[v(2:step:end)];
         if(plotsep<1)
         if(length(v)>1),contour3(vparm,vperpm,B,w);hold on;end
         contour3(vparm,vperpm,B,5);
         trapcone;hold off;
         view(0,90);
         axis(ax);
         set(gca,'Xtick',vparbox,'Ytick',vperpbox);
         set(gca,'XTickLabel',tickst);
         %axis equal tight
         title(['f(\rho=',num2str(r),',\theta=',num2str(theta),',vpar,vperp)']);
        
        subplot(3,4,[9 10]);
        if(~exist('sum_f')),sum_f=sum(distE+distEth);end
 %       distE=distE/sum_f;
 %       distEth=distEth/sum_f;
        %loglog(Ebin,distE,'-x')
        semilogy(Ebin,distE+distEth,'-x')
        title('f(E)');
        xlabel('E [eV]')
        
        subplot(3,4,[3 4]);
%        sloc=s((i-1)*nbinsrad+1:i*(nbinsrad));
%        p=polyfit(s,dens,4);%((i-1)*nbinsrad+1:i*(nbinsrad)),4);
%        pdens=polyval(p,s);%((i-1)*nbinsrad+1:i*(nbinsrad)));
        errorbar(sqrt(ss),dens,denserr,'-x');hold on;
        plot(sqrt(ss),densth,'k-x',sqrt(ss),dens+densth,'r-x');hold off
        title(['density']);
        legend('hot','thermal','tot');
%        xlabel('r');
        set(gca,'XTickLabel',[]);
        %ylabel('n(r)');
        a=axis;a(1)=0;a(2)=1;a(3)=0;axis(a);%a(4)=densmax;axis(a);
    
        subplot(3,4,[7 8]);
 %       p=polyfit(s,ppar,4);%((i-1)*nbinsrad+1:i*(nbinsrad)),4);
 %       pppar=polyval(p,s);%((i-1)*nbinsrad+1:i*(nbinsrad)));
 %       p=polyfit(s,pperp,4);%((i-1)*nbinsrad+1:i*(nbinsrad)),4);
 %       ppperp=polyval(p,s);%((i-1)*nbinsrad+1:i*(nbinsrad)));
 %       errorbar(a00*sqrt(ss),ppar,pparerr);hold on;
 %       errorbar(a00*sqrt(ss),pperp,pperperr,'k');hold off;
        plot(sqrt(ss),ppar,'-x',sqrt(ss),pperp,'-x');
        title(['Ppar, Pperp']);
 %       xlabel('r');
        set(gca,'XTickLabel',[]);
        %ylabel('Ppar(r), Pperp(r)');
        a=axis;a(1)=0;a(2)=1;a(3)=0;axis(a);%a(4)=pmax;axis(a);
        legend('Ppar','Pperp');
        
%         subplot(2,2,4);
%         plot(a00*sqrt(sloc(1:end-1)),ppperp(1:end-1)./pppar(1:end-1));
%         title(['Pperp/Ppar']);
%         xlabel('r');
%         ylabel('Pperp(r)/Ppar(r)');
%         a=axis;a(3)=0;a(4)=2;axis(a);

         subplot(3,4,[11 12]);
         plot(sqrt(ss),currt,sqrt(ss),currp,sqrt(ss),currt+currp,'-x');
         title('Parallel current');
         xlabel('\rho');
         %ylabel('j(r)');
         legend('trapped','passing','total')
         if(i>1);if(~exist('b')),b=axis;end
         a=axis;a(1)=0;amax=max(-min(a(3),0),a(4));a(3)=-amax;a(4)=amax;
         b(1)=0;b(2)=1;b(3)=-max(-b(3),-a(3));b(4)=max(b(4),a(4));b(4)=max(b(4),-b(3));
         axis(b);
         else a=axis;a(1)=0;a(2)=a00;axis(a);end
        
        drawnow;
    
    figure(findobj('Tag',num2str(f2)))
    subplot(1,2,1);
    %semilogy(Ebin,distE,'-x')
    semilogx(Ebin,sqrt(Ebin).*(distE+distEth).*gradient(Ebin),'-x')
    hold on;
    semilogx(Ebin,sqrt(Ebin).*distEth.*gradient(Ebin),'k-x')
    semilogx(Ebin,sqrt(Ebin).*distE.*gradient(Ebin),'r-x')
    %plot(Ebin(2:end),sqrt(Ebin(2:end)).*distE(2:end).*diff(Ebin),'-x')
    hold off
    legend('total','thermal','hot')
    title('E density');
    xlabel('E [eV]')
    
    subplot(1,2,2);
    %loglog(Ebin(2:end),-(diff(log(distE))./diff(Ebin)).^(-1),'-x'); %Teff
    %title('T_{eff}')
    semilogx(Ebin,cumsum(sqrt(Ebin).*(distE+distEth).*gradient(Ebin)),'-x')
    title('N(<E)');
    
    drawnow;
    if(wait>0),pause(wait),end;
    end
    end 
end
 
if(plotsum>1)
    figure; 
%    v=unique(A);
%    step=max(1,floor(length(v)/10));
%    w=[v(2:step:end)];
%    ivperpp=spline(vperpm,vperpm,linspace(min(vperpm),max(vperpm),100));f(length(v)>1),contour3(vparm,vperpm,A,w);hold on;end
%    contour3(vparm,vperpm,A,5);hold off;
    vparp=spline(vparm,vparm,linspace(min(vparm),max(vparm),49));
    vperpp=spline(vperpm,vperpm,linspace(min(vperpm),max(vperpm),7));
    [vp,vb]=meshgrid(vparp,vperpp);
    Ap=griddata(vparm,vperpm,A,vp,vb);
    for i=1:length(Ap(1,:)),As(:,i)=smooth(Ap(:,i),10,'sgolay',2);end
%    surfc(vp,vb,As);
%    set(gco,'LineStyle','none');
%    view(0,90);
%    contour(vp,vb,As,50);
%    set(gca,'Xtick',vparbox,'Ytick',vperpbox);grid on;axis(ax);
%    figure;
    contour(vparm,vperpm,A-Bdet,w);hold on;
    contour(vparm,vperpm,A-Bdet,5);hold off;
    %axis(ax);
    set(gca,'Xtick',vparbox,'Ytick',vperpbox);
%    set(gca,'XTickLabel',tickst);
    grid on;
%    axis equal
    axis(ax)
    title(['sum f(vpar,vperp)']);
end

if(plotsep==1)
         figure;
         if(length(v)>1),contour3(vparm,vperpm,B,w);hold on;end
        contour3(vparm,vperpm,B,5);
        trapcone;hold off;
        view(0,90);
        axis(ax);
        set(gca,'Xtick',vparbox,'Ytick',vperpbox);
        set(gca,'XTickLabel',tickst);
        %axis equal tight
        title(['f(\rho=',num2str(r),',\theta=',num2str(theta),',vpar,vperp)']);
        
%         Etot=cumsum(sqrt(Ebin).*(distE+distEth).*gradient(Ebin));
%         Eth=cumsum(sqrt(Ebin).*distEth.*gradient(Ebin));
%         Eh=cumsum(sqrt(Ebin).*distE.*gradient(Ebin));
        Etot=cumsum(sqrt(Ebin(2:end)).*(distE(2:end)+distEth(2:end)).*diff(Ebin));
        Eth=cumsum(sqrt(Ebin(2:end)).*distEth(2:end).*diff(Ebin));
        Eh=cumsum(sqrt(Ebin(2:end)).*distE(2:end).*diff(Ebin));
        hot_particles=Eh(end)/Etot(end)
        
        figure;
        set(gcf,'Position',[ 5   158   692   623]);
        subplot(3,1,[1]);
%        sloc=s((i-1)*nbinsrad+1:i*(nbinsrad));
%        p=polyfit(s,dens,4);%((i-1)*nbinsrad+1:i*(nbinsrad)),4);
%        pdens=polyval(p,s);%((i-1)*nbinsrad+1:i*(nbinsrad)));
        errorbar(sqrt(ss),dens,denserr,'-x');hold on;
        plot(sqrt(ss),densth,'k-x',sqrt(ss),dens+densth,'r-x');hold off
        title(['density']);
        legend('hot','thermal','tot');
%        xlabel('r');
         xlabel('\rho');
%        set(gca,'XTickLabel',[]);
        %ylabel('n(r)');
        a=axis;a(1)=0;a(2)=1;a(3)=0;axis(a);%a(4)=densmax;axis(a);
    
        subplot(3,1,[2]);
 %       p=polyfit(s,ppar,4);%((i-1)*nbinsrad+1:i*(nbinsrad)),4);
 %       pppar=polyval(p,s);%((i-1)*nbinsrad+1:i*(nbinsrad)));
 %       p=polyfit(s,pperp,4);%((i-1)*nbinsrad+1:i*(nbinsrad)),4);
 %       ppperp=polyval(p,s);%((i-1)*nbinsrad+1:i*(nbinsrad)));
 %       errorbar(a00*sqrt(ss),ppar,pparerr);hold on;
 %       errorbar(a00*sqrt(ss),pperp,pperperr,'k');hold off;
        plot(sqrt(ss),ppar,'-x',sqrt(ss),pperp,'-x');
        title(['Ppar, Pperp']);
 %       xlabel('r');
         xlabel('\rho');
  %      set(gca,'XTickLabel',[]);
        %ylabel('Ppar(r), Pperp(r)');
        a=axis;a(1)=0;a(2)=1;a(3)=0;axis(a);%a(4)=pmax;axis(a);
        legend('Ppar','Pperp');
        
%         subplot(2,2,4);
%         plot(a00*sqrt(sloc(1:end-1)),ppperp(1:end-1)./pppar(1:end-1));
%         title(['Pperp/Ppar']);
%         xlabel('r');
%         ylabel('Pperp(r)/Ppar(r)');
%         a=axis;a(3)=0;a(4)=2;axis(a);

         subplot(3,1,[3]);
         plot(sqrt(ss),currt,sqrt(ss),currp,sqrt(ss),currt+currp,'-x');
         title('Toroidal current');
         xlabel('\rho');
         %ylabel('j(r)');
         legend('trapped','passing','total')
         if(i>1);if(~exist('b')),b=axis;end
         a=axis;a(1)=0;amax=max(-min(a(3),0),a(4));a(3)=-amax;a(4)=amax;
         b(1)=0;b(2)=1;b(3)=-max(-b(3),-a(3));b(4)=max(b(4),a(4));b(4)=max(b(4),-b(3));
         axis(b);
         else a=axis;a(1)=0;a(2)=a00;axis(a);end
         
         figure;
         semilogx(Ebin,sqrt(Ebin).*(distE+distEth).*gradient(Ebin),'-x')
%          semilogx(Ebin,sqrt(Ebin).*(distE+distEth),'-x')
         hold on;
         semilogx(Ebin,sqrt(Ebin).*distEth.*gradient(Ebin),'k-x')
         semilogx(Ebin,sqrt(Ebin).*distE.*gradient(Ebin),'r-x')
%          semilogx(Ebin,sqrt(Ebin).*distEth,'k-x')
%          semilogx(Ebin,sqrt(Ebin).*distE,'r-x')
         %plot(Ebin(2:end),sqrt(Ebin(2:end)).*distE(2:end).*diff(Ebin),'-x')
         hold off
         legend('total','thermal','hot')
         title('E density');
         xlabel('E [eV]');
         
         figure;
         subplot(3,1,1);
         subint=length(Ebin)/nbinsvpar;
         j=0;clear Ebin_sub distE_sub
         for i=1:subint:length(Ebin);j=j+1;Ebin_sub(j)=sum(Ebin(i:i+subint-1))/subint;distE_sub(j)=sum(distE(i:i+subint-1)+distE(i:i+subint-1))/subint;end;
         loglog(Ebin_sub,-(gradient(log(distE_sub))./gradient(Ebin_sub)).^(-1),'-x'); %Teff
         title('T_{eff}')
         subplot(3,1,2);
         semilogx(Ebin,cumsum(sqrt(Ebin).*(distE+distEth).*gradient(Ebin)),'-x');
         aa=axis;aa(4)=1;axis(aa);
         title('N(<E)');
         subplot(3,1,3);
         semilogy(Ebin,distE+distEth,'-x');
         title('f(E)');
         
end
        