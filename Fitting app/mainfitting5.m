% fitting Li-Mg rejctions with SDEM model
function y = mainfitting5(x, cfli, cfmg, jv, Rli, Rmg, k)
    % units: mM, mM, um/s, 0-1, 0-1, LMH

    c1 = cfli*(1-Rli+Rli*exp(jv*3.6/k)); % li conc at feed surface
    c2 = cfmg*(1-Rmg+Rmg*exp(jv*3.6/k)); % mg conc at feed surface
    xf = [jv,c1,c2];
    R = SDEM(x/1e4,xf);
    r1 = 1/((1/Rli-1)*exp(-jv*3.6/k)+1); % li intrinsic rejection
    r2 = 1/((1/Rmg-1)*exp(-jv*3.6/k)+1); % mg intrinsic rejection
    y = [1*(R(1)-r1),50*(R(2)-r2)];
end


