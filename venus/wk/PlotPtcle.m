function [R,Z,phi,x,y,z]=PlotPtcle(nstep,choice)
%PLOTS PARTICLE ORBIT AS GIVEN IN FORT.98 FROM VENUS/tpr3d.f
%NEEDS ALSO FORT.37 FROM TERPSICHORE for flux surfaces

%read in data from fort.98
f98=load('fort.98');
R=f98(:,1);
Z=f98(:,2);
phi=f98(:,3);

if choice>0 %plot 2D orbit if choice >=1
    figure;
    plot(R,Z);
    axis equal;%grid
end

if choice>1 %plot flux surfaces in 2D plot if choice>=3
    fid8=fopen('fort.37','rt');
    [firstline,ncount]=fscanf(fid8,'%i %i',2);
    ni=firstline(1);           %number of radial grid points in TERPSICHORE
    nmodes=firstline(2);       %number of poloidal modes

    [fcoords,ncount]=fscanf(fid8,'%i %i %g %g %g %g %g %g %g %g',[10,nmodes*ni]);
    fclose(fid8);
    m=fcoords(1,:);             %vector with poloidal mode numbers
    n=fcoords(2,:);             %vector with toroidal mode numbers
    j=1;
    nmlb=nmodes*ni;
    for i=1:nmlb
        if (i>1) && (m(i)<m(i-1))
            j=j+1;              %j=radial grid point index
        end
        FR(m(i)+1,n(i)+1,j)=fcoords(3,i); %Fourier amplitudes R_{mn}(x)
        FZ(m(i)+1,n(i)+1,j)=fcoords(4,i); %Fourier amplitudes Z_{mn}(x)
        PHV(m(i)+1,n(i)+1,j)=fcoords(5,i);%Fourier amplitudes PHV_{mn}(x)
    end
    
    %Plotting the flux surfaces around the particle orbit
    nj=500 %number of poloidal grid points
    twopi=2*pi;
    dtheta=twopi/(nj-1);
    RH=zeros(ni,nj);
    ZH=zeros(ni,nj);
    for li=2:nmlb
        if (m(li)<m(li-1))
            mmodes=li-1;
            break;
        end
    end
    for line=1:mmodes
        for j=1:nj
            RH(1,j)=RH(1,j)+FR(m(line)+1,n(line)+1,1)*cos(m(line)*dtheta*(j-1));
            ZH(1,j)=ZH(1,j)+FZ(m(line)+1,n(line)+1,1)*sin(m(line)*dtheta*(j-1));
        end
    end
    for i=2:ni
        k=1;l=0;
        for line=mmodes+1:nmlb
            if (m(line)<m(line-1))
                if (l==1)
                    break;
                end
                k=k+1;
                if (i==k)
                    l=1;
                end
            end
            if (l==1)
                for j=1:nj
                    RH(i,j)=RH(i,j)+FR(m(line)+1,n(line)+1,i)*cos(m(line)*dtheta*(j-1));
                    ZH(i,j)=ZH(i,j)+FZ(m(line)+1,n(line)+1,i)*sin(m(line)*dtheta*(j-1));
                end 
            end
        end
    end


    for i=1:nstep:ni
        plot(RH(i,:),ZH(i,:),'k','Color',[0.8,0.8,0.8]);hold on;
    end
    plot(R,Z,'LineWidth',1.5);
    axis equal;
    hold off;
end

if choice>2
    [x,y,z]=Tok2Cart(R,Z,phi);
    figure;
    plot3(x,y,z);
    %grid;
    axis equal
end