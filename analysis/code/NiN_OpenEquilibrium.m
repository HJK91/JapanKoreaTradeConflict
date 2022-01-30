clc
clear
rng(2)
N = 2;
par.N       = N;    % # of countries
par.alpha   = 0.1;  % cost share of intermediate input aggregation
par.sigma   = 2;    % final good substitution parameter
par.eta     = 1.3;  % intermediate input substitution parameter
par.theta   = 4;    % shape parameter of productivity distribution
par.gamma   = 0.4;  % bargaining power of downstream firm

par.tau_d   = 1+rand(N); % iceberg trade cost between downstream firms-absorption place
par.tau_u   = 1+rand(N); % iceberg trade cost between uptream and downstream firms
par.tau_d(1:N+1:end) = 1;
par.tau_u(1:N+1:end) = 1;

par.nd      = [1 1]; % 1xN vector, # of downstream producers
par.nm      = 1;     % scalar, # of intermediate inputs
par.nu      = [1 1]; % 1xN vector, # of upstream producers


par.Td      = 1+rand(N,1);
par.Tu      = 1+rand(N,1);
par.w       = ones(N,1);
par.X       = ones(N,1);
par.P       = ones(N,1);

nd          = par.nd;
nu          = par.nu;
nm          = par.nm;
Td          = par.Td;
Tu          = par.Tu;
theta       = par.theta;
rng(1)

 zd = zeros(max(nd),N);
 zu = zeros(max(nu),N);
 for i = 1: N
    zd(1:nd(i),i) = (-log(rand([nd(i),1]))/Td(i)).^(1/theta); % zd: max(nd_i)xN matrix.
     zu(1:nu(i),i) = (-log(rand([nu(i),1]))/Tu(i)).^(1/theta); % zu: max(nu_i)xN matrix.
 end
 zd          = zd(:);
 zu          = zu(:);
 zd(zd==0)   = [];
 zu(zu==0)   = [];
zd = ones(sum(nd),1)*6;
zu = ones(sum(nu),1)*2;
qd_vec0     = zeros(sum(nd*N),1); % qd_vec: downstream firm's quantity decisions over N countries.
qumd_vec0   = zeros(sum(nu)*nm*sum(nd),1); % qumd_vec: upstream firm's quantity decisions over sum(nd) downstream firms
options     = optimset('Display','iter');
qd_sol      = fsolve( @(qd_vec) NiN_eq_v2(qd_vec,zd,zu,par),qd_vec0,options) ;
qd_sol      = exp(qd_sol); % convert back to positive numbers
qumd_sol    = fsolve(@(qumd_vec) NiN_bargaining_v2(qd_sol,qumd_vec,zd,zu,par),qumd_vec0, options);
qumd_sol = exp(qumd_sol);


function qd_sol = NiN_eq_v2(qd_vec,zd,zu,par)
    N = par.N;
    qd_vec  = exp(qd_vec);
    nd      = par.nd;
    nm      = par.nm;
    nu      = par.nu;
    qumd_vec0  = zeros(sum(nd)*sum(nu)*nm,1);
    options = optimset('Display','off');
    qumd_sol = fsolve(@(qumd_vec) NiN_bargaining_v2(qd_vec,qumd_vec,zd,zu,par),qumd_vec0, options);
    qumd_sol = exp(qumd_sol); % convert back to positive numbers
    pmd = pmd_vec(qumd_sol(1,1,1),qumd_sol,qd_vec,par,1,1,1);
    pmd = P_indx(pmd,par.eta);
    cd = cost_dowmstream(zd,pmd,par);
    % find the new qd
    qd_new = zeros(N,sum(nd));
    qd0 = zeros(N,1);
    for d=1:sum(nd)
        qd_new(:,d) = fsolve(@(qd) qd_foc(qd,qd_vec,cd,par,d),qd0,options); 
    end
    qd_sol    = exp(qd_new(:))-qd_vec;
end

function val = NiN_bargaining_v2(qd_vec,qumd_vec,zd,zu,par)
nd  = par.nd;
nm  = par.nm;
nu  = par.nu;
qumd_vec    = exp(qumd_vec);
qumd        = zeros(sum(nu),nm,sum(nd)); % find the bargaining outcome quantity.
qumd0       = 0;
options     = optimset('Display','off','MaxFunEvals',1e+4);            
for d   = 1:sum(nd)
    for u   = 1: sum(nu)
        for m   = 1:nm
            [qumd(u,m,d),~,flag] = fminsearch(@(q) BargainingObjective(q,qumd_vec,qd_vec,zd,zu,par,u,m,d),qumd0,options);
            if flag == 1
                qumd(u,m,d) = exp(qumd(u,m,d));
            else
                qumd(u,m,d) = 0;
            end
            
            % updated qumd
        end
    end
end
    qumd    = qumd(:);
    val     = qumd-qumd_vec;
end

function val = BargainingObjective(qumd,qumd_vec,qd_vec,zd,zu,par,u,m,d)
    gamma   = par.gamma;
    qumd    = exp(qumd);
    
    p_vec   = pmd_vec(qumd,qumd_vec,qd_vec,par,u,m,d); % intermediate input prices that firm d faces (nmxsum(nd) matrix).
    pi_u    = piumd(qumd,qumd_vec,zu,par,p_vec,u,m,d);
    pi_d    = pid(zd,par,p_vec,qd_vec,d); 
    p_vec_second    = pmd_vec(0,qumd_vec,qd_vec,par,u,m,d);
    
    pi_d_second     = pid(zd,par,p_vec_second,qd_vec,d);
    pi_u_second    = piumd(0,qumd_vec,zu,par,p_vec_second,u,m,d);
    ud = max(pi_d-pi_d_second,0);
    uu = max(pi_u-pi_u_second,0);
    val = -(ud)^gamma*(uu)^(1-gamma);
end

function val = pmd_vec(qumd,qumd_vec,qd_vec,par,u,m,d)
% find the vector of intermediate inputs given quantity of upstream goods.
% qumd is the endogenous level of input supplied, while qumd_vec is the
% quantity of intermediate inputs other than qumd. 
    eta     = par.eta;
    alpha   = par.alpha;
    w       = par.w;
    nd      = par.nd;
    nu      = par.nu;
    nm      = par.nm;
    N       = par.N;
    qd_sum  = sum(reshape(qd_vec,N,sum(nd)),1);  % each downstream firm's total production (1xsum(nd))
    qumd_matrix         = reshape(qumd_vec,sum(nu),nm,sum(nd));
    qumd_matrix(u,m,d)  = qumd;
    qmd_matrix          = reshape(sum(qumd_matrix,1),nm,sum(nd)); % aggregate across nu firms, nmxsum(nd) matrix.
    pmd_vec_old         = ones(nm,sum(nd));
    
    gap         = 1;
    constant1   = ((alpha/(1-alpha).*repelem(w,nd)).^(1-alpha))';  % (1xsum(nd))

    while gap >1e-6
        P           = P_indx(pmd_vec_old,eta); % 1xsum(nd) vector
        pmd_vec_new = (constant1.*qd_sum.*qmd_matrix.*P.^(eta-(1-alpha))).^(-1/eta);
        gap = norm(pmd_vec_new-pmd_vec_old)/norm(pmd_vec_old);
        pmd_vec_old = pmd_vec_new;
    end
    val = pmd_vec_new; % output: nmxsum(nd) matrix.
end

function val = piumd(qumd,qumd_vec,zu,par,p_vec,u,m,d)
    w       = par.w;
    nd      = par.nd;
    nu      = par.nu;
    nm      = par.nm;
    tau_u   = par.tau_u;
    qumd_mat= reshape(qumd_vec,sum(nu),nm,sum(nd));
    qumd_mat(u,m,d) = qumd;
    qumd_sales = reshape(qumd_mat(u,m,:),sum(nd),1); % sum(nd)x1 vector
%     nd_sum  = cumsum(nd);
    nu_sum  = cumsum(nu);
%     d_n     = sum(d>nd_sum)+1; % country index for downstream firm d.
    u_n     = sum(u>nu_sum)+1;
    tau_ud  = tau_u(u_n,:)'; % iceberg trade costs toward firm d.
    tau_ud  = repelem(tau_ud,nd);
    oper_val= max(p_vec(m,:)'-tau_ud*w(u_n)/zu(u),0); % the upstream firm operates only if it is profitable 
    val     = sum(oper_val.*qumd_sales);
end

function val= pid(zd,par,pmd_vec,qd_vec,d)
    N       = par.N;
    nd      = par.nd;
    sigma   = par.sigma;
    eta     = par.eta;
    pmd     = P_indx(pmd_vec,eta); % 1xsum(nd) vector
    % cd: each element represents the marginal cost of d incurred to ship
    % its good from d's country to the destination. Its size is as large as
    % Nxsum(nd).
    cd      = cost_dowmstream(zd,pmd, par); % Nxsum(nd) matrix
    % find the new qd.
    qd_old  = zeros(N,1);
    qd_mat = reshape(qd_vec,N,sum(nd));
    options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8);
%      options = optimset('Display','iter');
    qd_new  = exp(fsolve(@(qd) qd_foc(qd,qd_vec,cd,par,d),qd_old,options));
    qd_mat(:,d)  = qd_new;
    cd_mat  = reshape(cd,N,sum(nd));
    mktshr_d = qd_new./sum(qd_mat,2);
    temp_pi2 = mktshr_d./(sigma-mktshr_d).*cd_mat(:,d).*qd_mat(:,d);
    val = sum(temp_pi2);
end

function val = qd_foc(qd,qd_vec,cd,par,d)
    N       = par.N;
    nd      = par.nd;
    sigma   = par.sigma;
    cd_mat = reshape(cd,N,sum(nd));
    qd_mat = reshape(qd_vec,N,sum(nd));
    qd_mat(:,d) = exp(qd);
    P       = par.P; % Nx1 macro variable vector
    X       = par.X; % Nx1 macro variable vector
    Q       = sum(qd_mat,2); % N by 1 supply
    temp1   = Q.^(-1/sigma).*P.^(sigma-1).*X.^(1/sigma);
    temp2   = temp1-cd_mat(:,d); % Nx1
    temp3   = (1/sigma)*Q.^(-1/sigma-1).*P.^(sigma-1).*X.^(1/sigma);
    qd_new  = temp2./temp3;
    val     = qd_new-exp(qd);
end

function cd = cost_dowmstream(zd,pmd,par)
    % pmd: 1xsum(nd) vector
    % output: (sum(nd)*N)x1 vector, including trade costs.
    nd          = par.nd;   % Nx1 vecor
    w           = par.w;    % Nx1 vector
    N           = par.N;    % scalar
    alpha       = par.alpha; % scalar
    tau_df      = par.tau_d; % NxN downstream trade costs
    temp_tau    = repelem(tau_df',1,nd); % Nxsum(nd) matrix
    % 
    constant1   = (((1-alpha)/alpha)^alpha+(alpha/(1-alpha))^(1-alpha));
    cd          = constant1*repelem(w',1,nd).^(1-alpha).*pmd.^alpha./zd';
    cd          = repmat(cd,N,1);
    cd          = cd.*temp_tau;
end
function P = P_indx(p,eta)
    % input p is nmxsum(nd) matrix.
    % ouput is intermediate price index for downstream firms
    % (1xsum(nd))
    % eta is the elasticity of substitution
    P = sum(p.^(1-eta),1).^(1/(1-eta));
end
