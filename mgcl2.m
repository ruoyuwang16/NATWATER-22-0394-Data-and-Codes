function y = mgcl2(m)
    alpha = 2.0;
    beta0 = 0.4698*3/4;
    beta1 = 2.242*3/4;
    Ct = 0.00979*3/2^2.5;
    
    zc = 2;
    za = -1;
    vc = 1;
    va = 2;
    v = vc + va;
    mc = vc*m;
    ma = va*m;
    Im = 0.5*(ma*za^2+mc*zc^2);
    At = 0.3903; % (kg/mol)^0.5
    b = 1.2;
    
    ft = -At*Im.^0.5./(1+b*Im.^0.5);
    fg = -At*(Im.^0.5./(1+b*Im.^0.5)+2/b*log(1+b*Im.^0.5));
    Bt = beta0 + beta1*exp(-alpha*Im.^0.5);
    Bg = 2*beta0 + 2*beta1*(1-exp(-alpha*Im.^0.5).*(1+alpha*Im.^0.5-0.5*alpha^2*Im))./(alpha^2*Im);
    Cg = 1.5*Ct;
    
    y1 = abs(zc*za)*fg + m*(2*vc*va/v).*Bg + m.^2*(2*(vc*va)^1.5/v)*Cg; % ln(gamma)
    theta = 1 + abs(zc*za)*ft + m*(2*vc*va/v).*Bt + m.^2*(2*(vc*va)^1.5/v)*Ct; % theta
    y2 = -v*m*0.018.*theta; % ln(aH2O)

    y = theta;
end