function pi_u = NiN_upstreamprofit(zu,qum,pm,par)
    w = par.w;
    pi_u = (pm-w/zu)*qum;
end