clc
clear

par.alpha = 0.5;
par.sigma = 2;
par.sigma_prime = 3;
par.nd_H = 2;
par.nm = 2;
par.nu_H = 2;
par.nd_F = 2;
par.nu_F = 2;
par.eta = 1.3;
par.gamma = 0.4;
par.w = 1;
par.X = 1;
par.P = 1;
nd = par.nd_H;
nu = par.nu_H;
nm = par.nm;
theta = 4;
rng(1)
Td      = 1+rand(2,1);
Tu      = 1+rand(2,1);
zd_H    = (-log(rand([nd,1]))/Td(1)).^(1/theta);
zu_H    = (-log(rand([nu,nm]))/Tu(1)).^(1/theta);
qd_vec0 = rand(nd,1)*0.9588;
options = optimset('Display','iter');
qd_sol  = fsolve( @(qd_vec) NiN_eq_v2(qd_vec,zd_H,zu_H,par),qd_vec0,options) ;
qvec_0  = [qd_sol; ones(nu*nm*nd,1)];
qumd_sol = fsolve(@(qvec) NiN_bargaining_v2(qvec,zd_H,zu_H,par),qvec_0, options);
function qd_sol = NiN_eq_v2(qd_vec,zd_H,zu_H,par)
    w   = par.w;
    nd = par.nd_H;
    nm = par.nm;
    nu = par.nu_H;
    eta     = par.eta;
    alpha   = par.alpha;
    qvec_0  = [qd_vec; ones(nu*nm*nd,1)];
    options = optimset('Display','off');
    qumd_sol = fsolve(@(qvec) NiN_bargaining_v2(qvec,zd_H,zu_H,par),qvec_0, options);
    pmd = zeros(nd,1);
    for d = 1: nd
        temp_qumd1 = reshape(qumd_sol(nd+1:end),nu,nm,nd);
        temp_qumd2 = temp_qumd1(:,:,d);
        p_input_vec = pmd_vec(temp_qumd2(1,1),temp_qumd1(:),qvec_0(1:nd),par,1,1,d);
        pmd(d)         = P_indx(p_input_vec,eta);        
    end
    constant1   = (((1-alpha)/alpha)^alpha+(alpha/(1-alpha))^(1-alpha));
    cd          = constant1*w^(1-alpha).*pmd.^alpha./zd_H;
    % find the new qd
    qd_old      = ones(nd,1);
    gap = 1;
    while gap >1e-6
        qd_new  = qd_foc(qd_old,cd,par);
        gap     = norm(qd_new-qd_old)/norm(qd_old);
        qd_old  = qd_new;
    end    
    qd_sol    = qd_old(:)-qd_vec;
end

function val = NiN_bargaining_v2(qvec,zd_H,zu_H,par)
w   = par.w;
nd  = par.nd_H;
nm  = par.nm;
nu  = par.nu_H;
eta = par.eta;
alpha = par.alpha;
qd_vec      = qvec(1:nd);
qumd_vec    = qvec(nd+1:nd+nu*nm*nd);
qumd        = zeros(nu,nm,nd,1);
pmd = zeros(nd,1);
for d   = 1:nd
            % qd_sol as a function of qum_before
    for u   = 1: nu
            % qu_sol as a function of qd_before
        for m   = 1:nm
            qumd0=.01;
            options = optimset('Display','off','MaxFunEvals',1e+4);
            qumd(u,m,d) = fminsearch(@(q) BargainingObjective(q,qumd_vec,qd_vec,zd_H,zu_H,par,u,m,d),qumd0,options);
        end
    end
    p_input_vec = pmd_vec(qumd(u,m,d),qumd,qd_vec,par,u,m,d);
    pmd(d)         = P_indx(p_input_vec,eta);
end

    constant1   = (((1-alpha)/alpha)^alpha+(alpha/(1-alpha))^(1-alpha));
    cd          = constant1*w^(1-alpha).*pmd.^alpha./zd_H;
    % find the new qd
    qd_old      = ones(nd,1);
    gap = 1;
    while gap >1e-6
        qd_new  = qd_foc(qd_old,cd,par);
        gap     = norm(qd_new-qd_old)/norm(qd_old);
        qd_old  = qd_new;
    end    
    qumd    = qumd(:);
    val     = qumd-qumd_vec;
end

function val = BargainingObjective(qumd,qumd_vec,qd_vec,zd_H,zu_H,par,u,m,d)
    p_vec = pmd_vec(qumd,qumd_vec,qd_vec,par,u,m,d);
    pi_u = max([piumd(qumd,zu_H,par,p_vec,u,m),1e-8]);
    pi_d = pid(zd_H,par,p_vec,d); 
    
    p_vec_second = pmd_vec(0,qumd_vec,qd_vec,par,u,m,d);
    pi_d_second = pid(zd_H,par,p_vec_second,d);
    gamma = par.gamma;
    val = -(pi_d-pi_d_second)^gamma*(pi_u)^(1-gamma);
end

function val = pmd_vec(qumd,qumd_vec,qd_vec,par,u,m,d)
    eta = par.eta;
    alpha = par.alpha;
    w = par.w;
    nd = par.nd_H;
    nu = par.nu_H;
    nm = par.nm;
    qumd_matrix = reshape(qumd_vec,nu,nm,nd);
    qumd_matrix(u,m,d) = qumd;
    qumd_d = qumd_matrix(:,:,d); % nu x nm matrix
    qmd = sum(qumd_d)'; % aggregate across nu firms
    
    pmd_vec_old = ones(nm,1);
    gap = 1;
    constant1 = (alpha/(1-alpha)*w)^(1-alpha);
    while gap >1e-6
        P = P_indx(pmd_vec_old,eta);
        pmd_vec_new = (constant1*qd_vec(d)*qmd.*P.^(eta-(1-alpha))).^(-1/eta);
        gap = norm(pmd_vec_new-pmd_vec_old);
        pmd_vec_old = pmd_vec_new;
    end
    val = pmd_vec_new;
end

function val = piumd(qumd,zu_H,par,p_vec,u,m)
    w = par.eta;
    val = (p_vec(m)-w/zu_H(u,m))*qumd;
end

function val= pid(zd_H,par,pmd_vec,d)
    nd = par.nd_H;
    sigma = par.sigma;
    eta = par.eta;
    w = par.w;
    alpha = par.alpha;
    pmd =  P_indx(pmd_vec,eta);
    constant1 = (((1-alpha)/alpha)^alpha+(alpha/(1-alpha))^(1-alpha));
    cd = constant1*w^(1-alpha)*pmd^alpha./zd_H;
    % find the new qd, normalized to one.
    qd_old = ones(nd,1);
    gap = 1;
    while gap >1e-6
        qd_new = qd_foc(qd_old,cd,par);
        gap = norm(qd_new-qd_old)/norm(qd_old);
        qd_old = qd_new;
    end
    mktshr_d = qd_new(d)/sum(qd_new);
    val = mktshr_d/(sigma-mktshr_d)*cd(d)*qd_new(d);
    
end

function qd_new = qd_foc(qd_old,cd,par)
    sigma = par.sigma;
    P = par.P;
    X = par.X;
    Q = sum(qd_old);
    temp1 = Q^(-1/sigma)*P^(sigma-1)*X^(1/sigma)-cd;
    temp2 = sigma\Q^(-1/sigma-1)*P^(sigma-1)*X^(1/sigma);
    qd_new = max(temp1,1e-3)./temp2;
end

function P = P_indx(p,eta)
    % p is the vector of prices
    % eta is the elasticity of substitution
    P = sum(p.^(1-eta))^(1/(1-eta));
end
