function R = SDEM(P,X)
    
    Jv = X(1);
    c1 = X(2);
    c2 = X(3);
    c3 = c1 + 2*c2;

    p1 = P(1);
    p2 = P(2);
    p3 = p1;

    z1 = 1;
    z2 = 2;
    z3 = -1;

    c0 = c1+c2+c3;
    u0 = (z1^2*c1+z2^2*c2+z3^2*c3)/c0;

    pi = z1*z2*z3;
    sigma = z1+z2+z3;

    b02 = (z2^2-z3^2)/p1-(z1^2-z3^2)/p2+(z1^2-z2^2)/p3;
    b12 = z1*(z2^2-z3^2)/p1-z2*(z1^2-z3^2)/p2+z3*(z1^2-z2^2)/p3;
    b01 = (z2-z3)/p1-(z1-z3)/p2+(z1-z2)/p3;
    b11 = z1*(z2-z3)/p1-z2*(z1-z3)/p2+z3*(z1-z2)/p3;
    b_11 = 1/z1*(z2-z3)/p1-1/z2*(z1-z3)/p2+1/z3*(z1-z2)/p3;
    
    fun = @(x)upsolver(x,pi,b02,b12,b01,b11,b_11,u0,c0,z1,z2,z3,Jv);
    x0 = -pi*(u0^2*b_11+u0*b02+pi*b01)/(u0^2*b01+u0*b12+pi*b11)+1e-6;
    up = fsolve(fun,x0);

    m1 = (pi*b02+b12*up)^2-4*pi*(pi*b01+b11*up)*(pi*b_11+b01*up);
    m2 = pi*b02+b12*up;
    m3 = 2*(pi*b_11+b01*up);
    fp = (m1^0.5-m2)/m3;
    fm = (m1^0.5+m2)/m3;

    n1 = (fm+u0)/(fm+up);
    n2 = fm/(fm+fp);
    n3 = (fp-u0)/(fp-up);
    n4 = fp/(fm+fp);
    cp = c0*n1^n2*n3^n4;

    cp1 = cp*(up+z2*z3)/(z1-z2)/(z1-z3);
    cp2 = cp*(up+z1*z3)/(z2-z1)/(z2-z3);
    cp3 = cp*(up+z2*z1)/(z3-z2)/(z3-z1);

    R1 = 1-cp1/c1;
    R2 = 1-cp2/c2;
    R3 = 1-cp3/c3;
    R = [R1,R2];

end