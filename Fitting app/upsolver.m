function y=upsolver(x,pi,b02,b12,b01,b11,b_11,u0,c0,z1,z2,z3,Jv)
    m1 = (pi*b02+b12*x)^2-4*pi*(pi*b01+b11*x)*(pi*b_11+b01*x);
    m2 = pi*b02+b12*x;
    m3 = 2*(pi*b_11+b01*x);
    fp = (m1^0.5-m2)/m3;
    fm = (m1^0.5+m2)/m3;

    n1 = (fm+u0)/(fm+x);
    n2 = fm/(fm+fp);
    n3 = (fp-u0)/(fp-x);
    n4 = fp/(fm+fp);
    cp = c0*n1^n2*n3^n4;

    jv = (c0/cp-1)*(z1-z2)*(z2-z3)*(z1-z3)/(pi*b_11+b01*x);

    y = jv-Jv;
end