% set parameters
par.alpha = 0.5;
par.sigma = 2;
par.sigma_prime = 3;
par.n_d = 2;
par.n_m = 2;
par.n_u = 2;
par.eta = 1.3;
par.gamma = 0.4;
par.w = 1;
par.X = 1;
par.P = 1;

% from 1 to 2 countries: a single simulation point.
nd = par.n_d;
nu = par.n_u;
nm = par.n_m;
theta = 4;
rng(1215)
Td      = 1+rand(2,1);
Tu      = 1+rand(2,1);
zd_H = (-log(rand([nd,1]))/Td(1)).^(1/theta);
zu_H = (-log(rand([nu,nm]))/Tu(1)).^(1/theta);

pqvec0 = ones(nu*nm*nd+nd+1,1);
options = optimset('Display','iter');
q_sol = fsolve(@(qvec) cournot_competition(zd_H,zu_H,qvec,par),pqvec0,options);

% qvec_old = ones(nu*nm*nd+nd,1)*.1;
% gap = 1;
% while gap>1e-6
%     qvec_new = cournot_competition(zd_H,zu_H,qvec_old,par);
%     gap = norm(qvec_new-qvec_old)/norm(qvec_old);
%     qvec_old = qvec_new;
% end

function val = cournot_competition(zd_H,zu_H,qvec,par)

nd = par.n_d;
nm = par.n_m;
nu = par.n_u;
qd_before = qvec(1:nd);
qum_before = reshape(qvec(nd+1:nu*nm*nd+nd),nu,nm,nd);
qd_sol = zeros(nd,1);
qum_sol = zeros(nu,nm,nd);
options = optimset('Display','off');
for d   = 1:nd
        % qd_sol as a function of qum_before
    qd_sol_temp = zeros(nu,nm);
    for u   = 1: nu
        % qu_sol as a function of qd_before
        for m   = 1:nm
            qd0 = qd_before(d);
            qum_sol(u,m,d) = fsolve(@(qd) tempfun_qum(qd,zd_H,zu_H,qum_before,qd_before,par,u,m,d),qd0,options);
            
            qum0 = qum_before(u,m,d);
            qd_sol_temp(u,m) = fsolve(@(qum) tempfun_qdm(qum,zd_H,zu_H,qum_before,qd_before,par,u,m,d), qum0,options);
        end
    end
    qd_sol(d) = mean(qd_sol_temp,'all');
end



q_sol(1:nd) = qd_sol;
q_sol(nd+1:nd+nm*nu*nd) = qum_sol(:);
val = q_sol;
end

function val = tempfun_qdm(qum,zd_H,zu_H,qum_before,qd_before,par,u,m,d)
    
    q_d_temp = qd_before;
    q_d_temp(d) = [];
    
    % Given qum, find the bargaining outcome of qdm.
    q_m_temp = qum_before(:,:,d);
    q_m_temp(:,m) = [];
    q_m_temp = Q_agg(q_m_temp,par.eta);
    q_um_temp = qum_before(:,m,d);
    q_um_temp(u) = [];
    zd_bargain = zd_H(d);
    zu_bargain = zu_H(u,m);
    [qum_sol,~] = NiN_bargaining(zd_bargain,zu_bargain,q_um_temp,q_m_temp,q_d_temp,par);
    val = qum-qum_sol;
    
end

function val = tempfun_qum(qd,zd_H,zu_H,qum_before,qd_before,par,u,m,d)
    % Given qum, find the bargaining outcome of qdm.
    q_d_temp = qd_before;
    q_d_temp(d) = [];    
    
    q_m_temp = qum_before(:,:,d);
    q_m_temp(:,m) = [];
    q_m_temp = Q_agg(q_m_temp,par.eta);
    
    q_um_temp = qum_before(:,m,d);
    q_um_temp(u) = [];
    
    zd_bargain = zd_H(d);
    zu_bargain = zu_H(u,m);
    [~,qd_sol] = NiN_bargaining(zd_bargain,zu_bargain,q_um_temp,q_m_temp,q_d_temp,par);
    val = qd-qd_sol;
    
end

function val = factor_price(p,zd,qmvec,par)
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
    Q = sum(q.^(eta/(eta-1))).^(eta/(eta-1));
end

