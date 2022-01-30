function [qum_sol,qd_sol] = NiN_bargaining(zd,zu,q_um,q_m,q_d,par)


% % set qum,zd,q_um,q_m,q_d
% qum     = .11;
% zd      = 1;
% zu      = 20;
% q_um    = [0.22; 0.43];
% q_m     = [0.8;0.3];
% q_d     = [0.3;0.9; 0.43];
% [pi_d,pi_u] = NiN_profits(zd,zu,qum,q_um,q_m,q_d,par);
% [pi_d_second,~] = NiN_profits(zd,zu,0,q_um,q_m,q_d,par);

%%
qum0 = 0.1;
qum_sol = fminsearch(@(qum) BargainingObjective(zd,zu,qum,q_um,q_m,q_d,par),...
    qum0);
qum_sol = exp(qum_sol);
[~,~,qd_sol] = NiN_profits(zd,zu,qum_sol,q_um,q_m,q_d,par);
end
function val = BargainingObjective(zd,zu,qum,q_um,q_m,q_d,par)
    [pi_d,pi_u,~] = NiN_profits(zd,zu,qum,q_um,q_m,q_d,par);
    [pi_d_second,~] = NiN_profits(zd,zu,0,q_um,q_m,q_d,par);
    gamma = par.gamma;
    val = -(pi_d-pi_d_second)^gamma*(pi_u)^(1-gamma);
end

