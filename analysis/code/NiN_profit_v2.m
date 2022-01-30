function [pi_d,pi_u,qd_sol] = NiN_profit_v2(zd_H,zu_H,qumd_vec,par,u,m,d)
% find the upstream and downstream firm's profit as a function of
% intermediate input provision of the upstream firm, given other firms'
% input supply, other material supply, and competing downstream firms'
% supply to the market.
% qum: quantity supplied by upstream firm u on input m.
% q_um: quantity supplied on input m by unstream firms other thna suplier
% u.
% zd: productivity of the downstream firm.
% q_m: quantity of inputs other than input m supplied to the downstream
% firm.
% q_d: quantity of competing downstream firms supplied to the market 
% par: various paramters
    
    qum = exp(qumd_vec(:,:,d)); % for a positive number solution.
    q_um = qum(:,m);
    q_um(u) = [];
    q_m = qum;
    q_m(:,m) = [];
    qum = qum(u,m);
    zd = zd_H(d);
    zu = zu_H(u,m);
    qd_old = qd_vec(d);
    q_d = qd_vec;
    q_d(d) = [];
    gap = 1;
    while gap>1e-6
        [pm,p_m] = price_upstream(zd,qum,q_um,q_m,qd_old,par);
        cd = unit_cost(pm,p_m,par);        
        [~,qd_new] = cournot_markup(cd,q_d,par);
        gap = norm(qd_new-qd_old);
        qd_old = qd_new;
    end
    qd_sol = qd_new;
%     qd_sol = fsolve(@(qd) eq_qd(zd,qum,q_um,q_m,q_d,qd,par),qd_old);
    [pm,p_m] = price_upstream(zd,qum,q_um,q_m,qd_sol,par);
    cd = unit_cost(pm,p_m,par);    
    mktshare = qd_sol/sum([qd_sol;q_d]);
    pi_d = (mktshare/(par.sigma-mktshare))*(cd/zd)*qd_sol;
    w = par.w;
    pi_u = max([(pm-w/zu)*qum,0]);
end

function val = eq_qd(zd,qum,q_um,q_m,q_d,qd_old,par)
    [pm,p_m] = price_upstream(zd,qum,q_um,q_m,qd_old,par);
    cd = unit_cost(pm,p_m,par);
    [~,qd_new] = cournot_markup(cd,q_d,par);
    val = norm(qd_new-qd_old)/norm(qd_old);
end

function [epsilond,qd] = cournot_markup(cd,q_d,par)
    sigma = par.sigma;
    sigma_prime = par.sigma_prime;
    X = par.X;
    P = par.P;
    constant = P^((sigma-1)/sigma)*X^(1/sigma);
    omega = 1/sigma-1/sigma_prime;
    qd_old = 0.5;
    gap = 1;
    while gap >1e-6
        q = [qd_old;q_d];
        Q = Q_agg(q,sigma_prime);
        mkt_share = (qd_old/Q)^((sigma_prime-1)/sigma_prime);
        epsilond = 1/(mkt_share/sigma+(1-mkt_share)/sigma_prime);
        qd_new = (epsilond/(epsilond-1)*cd*(Q^(-omega)*constant))^(-1/sigma_prime);
        gap = norm(qd_old-qd_new)/norm(qd_old);
        qd_old = qd_new;
    end
    qd_sol = qd_new;
%     qd0 = 0.5;
%     qd_sol = fsolve(@(qd) ...
%         iter_qd(cd,qd,q_d,sigma,sigma_prime,constant),qd0);
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

function [pm,p_m] = price_upstream(zd,qum,q_um,q_m,qd_old,par)
    nm = par.n_m; % number of different materials
    
   
    qm  = sum([qum;q_um]);
    qmvec   = [qm;q_m];
    eta = par.eta;
    alpha =par.alpha;
    w = par.w;
    constant1 = (alpha/(1-alpha))^(1-alpha);
    constant2 = w^(1-alpha)*qd_old/zd;
    gap = 1;
    p_old = ones(nm,1);
    while gap >1e-6
        P = P_indx(p_old,eta);
        p_new = (constant1*constant2*P^(eta-(1-alpha))./qmvec).^(1/eta);
        gap = norm(p_old-p_new)/norm(p_old);
        p_old = p_new;
    end
    p_sol = p_new;
%     p0   = ones(nm,1);
%     p_sol = fsolve(@(p) factor_price(p,zd,qmvec,qd_old,par),p0);
    
    pm = p_sol(1);
    p_m = p_sol(2:end);
end

function val = factor_price(p,zd,qmvec,qd_old,par)
    eta = par.eta;
    alpha =par.alpha;
    w = par.w;
    constant1 = (alpha/(1-alpha))^(1-alpha);
    constant2 = w^(1-alpha)*qd_old/zd;
    P = P_indx(p,eta);
    val = p.^eta-constant1*constant2*P^(eta-(1-alpha))./qmvec;
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