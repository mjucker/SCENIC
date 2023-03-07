function [SmNew]=TruncExp(Shaf,Sm,r)

si=size(Sm);
SmNew=zeros(3,si(2));
theta=linspace(0,2*pi,100);
for k=1:si(2)
    R=zeros(1,100);
    for i=2:si(1)
        R=R+Sm(i,k)*cos((i-1)*theta);
    end
    R=R-Shaf(k)*ones(1,100)+r(k)*cos(theta);

    % --- Create fit
    ok_ = ~(isnan(theta) | isnan(R));
    ft_ = fittype({'-1', 'cos(x)', 'cos(2*x)'},...
         'dependent',{'y'},'independent',{'x'},...
         'coefficients',{'a0', 'a1', 'a2'});

    % Fit this model using new data
    cf = fit(theta(ok_)',R(ok_)',ft_);

    % Get coefficients
    SmNew(1,k)=cf.a0;
    SmNew(2,k)=cf.a1-r(k);
    SmNew(3,k)=cf.a2;
end