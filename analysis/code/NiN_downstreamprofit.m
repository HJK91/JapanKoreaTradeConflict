function [pi_d,pm] = NiN_downstreamprofit(zd,qum,q_um,q_m,q_d,par)
% find the firm's profit given sourcing quantity choice and competitor's quantity
% supply.
% qum: quantity supplied by upstream firm u on input m.
% q_um: quantity supplied on input m by unstream firms other thna suplier
% u.
% zd: productivity of the downstream firm.
% q_m: quantity of inputs other than input m supplied to the downstream
% firm.
% q_d: quantity of competing downstream firms supplied to the market 
% par: various paramters 
    [pm,p_m] = price_upstream(qum,q_um,q_m,par);
    cd = unit_cost(pm,p_m,par);
    [epsilond,qd] = cournot_markup(cd,q_d,par);
    pi_d = (epsilond-1)\(cd/zd)*qd;
end

function [epsilond,qd] = cournot_markup(cd,q_d,par)
    sigma = par.sigma;
    sigma_prime = par.sigma_prime;
    X = par.X;
    P = par.P;
    constant = P^((sigma-1)/sigma)*X^(1/sigma);
    qd0 = 0.5;
    qd_sol = fsolve(@(qd) ...
        iter_qd(cd,qd,q_d,sigma,sigma_prime,constant),qd0);
    Q = Q_agg([qd_sol;q_d],sigma_prime);
    mkt_share = (qd_sol/Q)^((sigma_prime-1)/sigma_prime);
    epsilond = 1/(mkt_share/sigma+(1-mkt_share)/sigma_prime);
    qd  = qd_sol;
end

function val = iter_qd(cd,qd,q_d,sigma,sigma_prime,constant)
    omega = 1/sigma-1/sigma_prime;
    q = [qd;q_d];
    Q = Q_agg(q,sigma_prime);
    mkt_share = (qd/Q)^((sigma_prime-1)/sigma_prime);
    epsilond = 1/(mkt_share/sigma+(1-mkt_share)/sigma_prime);
    val = epsilond/(epsilond-1)*cd-qd^(-1/sigma_prime)...
            *Q^(-omega)*constant;
end

function [pm,p_m] = price_upstream(qum,q_um,q_m,par)
    eta = par.eta;
    n_u = par.n_u; % number of upstream firms
    p_old  = ones(n_u+1,1);
    qm  = sum([qum;q_um]);
    q   = [qm;q_m];
    Q   = Q_agg(q,eta);

    gap = 1;
    while gap>1e-6
        P = P_indx(p_old,eta);
        p_new = P*((q/Q).^(-1/eta));
        gap = p_new-p_old;
        p_old = p_new;
    end
    
    pm = p_old(1);
    p_m = p_old(2:end);
end


function cd = unit_cost(pm,p_m,par)
    alpha = par.alpha;
    eta = par.eta;
    w   = par.w;
    P   = P_indx([pm;p_m],eta);
    temp = ((1-alpha)/alpha)^alpha+(alpha/(1-alpha))^(1-alpha);
    cd = temp*w^(1-alpha)*P^alpha;
end

function P = P_indx(p,eta)
    % p is the vector of prices
    % eta is the elasticity of substitution
    P = sum(p.^(1-eta))^(1/(1-eta));
end

function Q = Q_agg(q,eta)
    % p is the vector of prices
    % eta is the elasticity of substitution
    Q = sum(q.^(eta/(eta-1)))^(eta/(eta-1));
end