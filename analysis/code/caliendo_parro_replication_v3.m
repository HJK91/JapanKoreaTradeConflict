clc
clear
S = 6;
J = 5;
par.T       = rand(J*S,1)*0.2+1;
par.tau     = rand(S*J^2,1)*0.1+1;
par.theta   = 5+rand(S,1)*0.1;
par.sigma   = par.theta+1-rand(S,1);
par.gamma   = sort(rand(J,S+1,S),2);
par.gamma(:,1,:)   = 0;
par.gamma   = diff(sort(par.gamma,2),1,2);
par.gamma   = par.gamma(:);
par.gamma_sum   = sum(reshape(par.gamma,J,S,S),2);
par.gamma_sum   = par.gamma_sum(:);
par.upsilon = prod(reshape(par.gamma.^(-par.gamma),J,S,S),2);
par.upsilon = par.upsilon(:).*(1-par.gamma_sum).^(-(1-par.gamma_sum));
par.alpha   = diff([zeros(1,J); sort(rand(S-1,J),1) ; ones(1,J)])';
par.alpha   = par.alpha(:);
par.kappa   = gamma(1+(1-par.sigma)./par.theta).^(1./(1-par.sigma));

L_0     = ones(J,1);   % Jx1 data
D_0     = zeros(J,1);        % Jx1 data
options = optimoptions('fsolve','Display','iter','StepTolerance',1e-8,...
    'MaxFunctionEvaluations',1e+6,'MaxIterations',5e+3,'OptimalityTolerance',1e-10);
N = 5;
w_sol = zeros(J,N);
fval = zeros(J,N);
flag = zeros(N,1);
for i = 1: N
w_0 = rand(J,1);
% wX_0 = ones(J*S+J-1,1);
[w_sol(:,i), fval(:,i), flag(i)]   = fsolve(@(w) h(w,par,S,J,L_0,D_0), w_0, options);

end
w_sol = exp(w_sol);
w_sol = w_sol./w_sol(1,:);
function z = g(w,cP,par,S,J)
    
    z       = zeros(2*S*J,1);
    T       = par.T;
    tau     = par.tau;
    theta   = par.theta;
    gamma   = par.gamma;
    gamma_sum = par.gamma_sum;
    upsilon = par.upsilon;
    kappa   = par.kappa;
    c           = exp(cP(1:S*J));
    P           = exp(cP(S*J+1:2*S*J));
    
    c_temp = mididx(c,S,J,J);
    T_temp = mididx(T,S,J,J);
    THETA       = T_temp.*((c_temp.*tau).^(-kron(theta,ones(J^2,1))));
    THETA_sum   = sum(reshape(THETA,J,S*J),1)';
    P_rep       = mididx(P,S,S,J);
    P_gamma     = P_rep.^gamma;
    P_prod      = prod(reshape(P_gamma,J,S,S),2);
    P_prod      = P_prod(:);
    z(1:S*J) = upsilon.*(repmat(w,S,1).^(1-gamma_sum)).*P_prod; % c
    
    kappa_rep = kron(kappa,ones(J,1));
    theta_rep = kron(theta,ones(J,1));
    
    z(S*J+1:2*S*J) = kappa_rep.*(THETA_sum).^(-1./theta_rep); % P
    
   
end

function z = h(wX,par,S,J,L,D)
% w: J-1 vector (first country's wage is numeraire)
% X: SxJ vector
    T       = par.T;
    tau     = par.tau;
    theta   = par.theta;
    gamma   = par.gamma;
    %gamma_sum = par.gamma_sum;
    %upsilon = par.upsilon;
    alpha   = par.alpha;
    %kappa   = par.kappa;
    w = exp(wX(1:J));
    w = w/w(1);
%     X = exp(wX(J:S*J+J-1));
    X = exp(randn(S*J,1));
    options = optimoptions('fsolve','StepTolerance',1e-8,...
    'MaxFunctionEvaluations',1e+6,'MaxIterations',5e+3);
    cP_prev     = [randn(S*J,1); randn(S*J,1)*0.1];
    gap     = 1;
    tol     = 1e-8;
    while gap>tol
        cP_new = log(g(w,cP_prev,par,S,J));
        gap = norm(exp(cP_prev)-exp(cP_new));
        cP_prev = cP_new;
    end
    
    c_temp = mididx(exp(cP_new(1:S*J)),S,J,J);
    T_temp = mididx(T,S,J,J);
    THETA       = T_temp.*((c_temp.*tau).^(-kron(theta,ones(J^2,1))));
    THETA_sum   = sum(reshape(THETA,J,S*J),1)';
    pi_temp =  THETA./lastidx(THETA_sum,S,J,J);
    
%     X_sol   = fsolve(@(XX) TBandED(w,XX,pi_temp,par,S,J,L,D), X, options);
%     X_sol = exp(X_sol);
%      Y = sum(reshape(lastidx(X_sol,S,J,J).*pi_temp,J,J,S),2);
%      Y =Y(:);
     gap = 1;
     X_prev=X;
     while gap > 1e-8
        
        Y = sum(reshape(lastidx(X_prev,S,J,J).*pi_temp,J,J,S),2);
        Y =Y(:);
        Y_gamma = sum(reshape(mididx(Y,S,S,J).*gamma,J,S,S),3);
        Y_gamma = Y_gamma(:);
        X_new = (Y_gamma+alpha.*repmat(w.*L+D,S,1));
        gap = norm(X_prev-X_new);
        X_prev = X_new;
     end
%     Y_gamma = sum(reshape(mididx(Y,S,S,J).*gamma,J,S,S),3);
%     Y_gamma = Y_gamma(:);
%     z_temp = X-(Y_gamma+alpha.*repmat([1;w].*L+D,S,1));
%     z(1:S*J-1,1) = z_temp(1:end-1);
     X_sum = sum(reshape(X_new,J,S),2);
     Y_sum = sum(reshape(Y,J,S),2);
     z_temp = X_sum-(Y_sum+D);
    z(1:J) = z_temp;
end

function z=TBandED(w,X,pi,par,S,J,L,D)
    gamma = par.gamma;
    alpha = par.alpha;
    X = exp(X);
    Y = sum(reshape(lastidx(X,S,J,J).*pi,J,J,S),2);
    Y =Y(:);
    Y_gamma = sum(reshape(mididx(Y,S,S,J).*gamma,J,S,S),3);
    Y_gamma = Y_gamma(:);
    z= (Y_gamma+alpha.*repmat(w.*L+D,S,1));
    
%     z(1:S*J-1,1) = z_temp(1:end-1);
%     X_sum = sum(reshape(X,J,S),2);
%     Y_sum = sum(reshape(Y,J,S),2);
%     z_temp = X_sum-(Y_sum+D);
%     z(S*J) = z_temp(1);
end

function z = mididx(x,idx1,idx2,idx3)
% expand the middle index
x_temp = repmat(reshape(x,idx3,idx1),idx2,1);
z = x_temp(:);
end

function z = lastidx(x,idx1,idx2,idx3)
% expand the last index
% x is (idx1*idx2) by 1 vector
x_temp = repmat(reshape(x,1,idx1*idx2),idx3,1);
z = x_temp(:);
end
