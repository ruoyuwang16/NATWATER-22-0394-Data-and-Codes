function [plilist,pmglist] = main(cfli, cfmg, jv, Rli, Rmg, k)

    plilist = zeros(1,length(cfli));
    pmglist = zeros(1,length(cfli));
    
    for i = 1:length(cfli)
        x0 = [100,5]*1e4; % initial guess
        options  = optimoptions(@lsqnonlin,'Display','iter');
        f = @(x)mainfitting5(x, cfli(i), cfmg(i), jv(i), Rli(i), Rmg(i), k);
        p = lsqnonlin(f,x0,zeros(size(x0)),[],options); % solve nonlinear regression
        plilist(i) = p(1)/1e4;
        pmglist(i) = p(2)/1e4;
    end
end
