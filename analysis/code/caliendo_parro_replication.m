% Set the number of equations and the number of unknowns
clc
clear
S = 4;
J = 3;

% Unknowns: pi_ij^s c_j^s P_j^s w_j X_j^s

%% Set parameters
par.T       = rand(J*S,1)*0.2+0.1;
par.tau     = rand(S*J^2,1)*0.1+1;
par.theta   = 1+rand(S,1)*0.1;
par.sigma   = par.theta+1-rand(S,1);
par.gamma   = sort(rand(J,S+1,S),2);
par.gamma(:,1,:)   = 0;
par.gamma   = diff(sort(par.gamma,2),1,2);
par.gamma   = par.gamma(:);

par.gamma_sum   = sum(reshape(par.gamma,J,S,S),2);
par.gamma_sum   = par.gamma_sum(:);

par.upsilon = prod(reshape(par.gamma.^(-par.gamma),J,S,S),2);
par.upsilon = par.upsilon(:).*par.gamma_sum.^(-par.gamma_sum);
par.alpha   = diff([zeros(1,J); sort(rand(S-1,J),1) ; ones(1,J)]);
par.alpha   = par.alpha(:);
par.kappa   = gamma(1+(1-par.sigma)./par.theta).^(1./(1-par.sigma));
%% convention: sector first, importing country next, exporting country next
% simultaneous equations
% Total # of equations
% clc
% 
% 
% N = 2;
% z_sol = zeros(length(Z_0),N);
% fval = zeros(length(Z_0),N);
% flag = zeros(N,1);
% iter = zeros(N,1);
% 
% for i = 1: N
% Z_0     = [pi_0; exp(randn(length(Z_0)-length(pi_0),1))/100];
% pi_0    = cat(1,cat(1,zeros(1,J,S),sort(rand(J-1,J,S),1)),ones(1,J,S));
% pi_0    = diff(pi_0,1,1);
% pi_0    = reshape(pi_0(1:end-1,:,:),S*J*(J-1),1); % Jx(J-1)xS unknowns
% 
% c_0     = randn(S*J,1);      % JxS unknowns
% w_0     = randn(J-1,1);   % J-1 unknowns
% P_0     = randn(S*J,1);   % JxS unknowns
% X_0     = randn(S*J,1);      % JxS unknowns
% L_0     = rand(J,1);   % Jx1 data
% D_0     = zeros(J,1);        % Jx1 data
% 
% Z_0      = [pi_0; c_0; w_0; P_0; X_0]; 
%  f1      = @(z) f(z,par,S,J,L_0,D_0);
% % [z_sol(:,i), fval(:,i), flag(i), iter(i), fjacinv]   = broyden(f1, Z_0);
% %
% clc
% options = optimoptions('fsolve','Display','iter','StepTolerance',1e-8,...
%     'MaxFunctionEvaluations',1e+6,'MaxIterations',5e+3);
% [z_sol(:,i), fval(:,i), flag(i), ~, fjacinv]   = fsolve(@(z) f(z,par,S,J,L_0,D_0), Z_0, options);
% end
%%
clc
L_0     = zeros(J,1);   % Jx1 data
D_0     = zeros(J,1);        % Jx1 data
%
options = optimoptions('fsolve','Display','iter','StepTolerance',1e-8,...
    'MaxFunctionEvaluations',1e+6,'MaxIterations',5e+3);
wX_0 = [zeros(J-1,1);rand(S*J,1)*1.5];
% wX_0 = ones(J*S+J-1,1);
[z_sol, fval, flag]   = fsolve(@(wX) h(wX,par,S,J,L_0,D_0), wX_0, options);
hh = @(wX) h(wX,par,S,J,L_0,D_0);
% [z_sol, fval, flag]   = broyden(hh, wX_0);
function z = f(Z,par,S,J,L,D)
    
    z       = zeros(S*J*(J+2)+J-1,1);
    T       = par.T;
    tau     = par.tau;
    theta   = par.theta;
    gamma   = par.gamma;
    gamma_sum = par.gamma_sum;
    upsilon = par.upsilon;
    alpha   = par.alpha;
    kappa   = par.kappa;
    
    pi = exp(Z(1:S*J*(J-1)))./(1+exp(Z(1:S*J*(J-1))));
    pi = reshape(pi,J-1,J,S);
    pi = cat(1,pi,1-sum(pi,1));
    pi = pi(:);
    
    c = exp(Z(S*J*(J-1)+1:S*J^2));
    w = exp([1; Z(S*J^2+1:S*J^2+J-1)]);
    P = exp(Z(S*J^2+J:S*J*(J+1)+J-1));
    X = exp(Z(S*J*(J+1)+J:S*J*(J+2)+J-1));
    
    THETA       = repmat(T,J,1).*((repmat(c,J,1).*tau).^(-kron(theta,ones(J^2,1))));
    THETA_sum   = sum(reshape(THETA,J,S*J),1)';
    z_temp      = pi-THETA./kron(THETA_sum,ones(J,1));
    z_temp(J:J:end)     = [];
    z(1:S*J*(J-1))      = z_temp; % pi

    P_rep       = repmat(P,S,1);
    P_gamma     = P_rep.^gamma;
    P_prod      = prod(reshape(P_gamma,J,S,S),2);
    P_prod      = P_prod(:);
    
    z(S*J*(J-1)+1:S*J^2) = c - upsilon.*(repmat(w,S,1).^(1-gamma_sum)).*P_prod; % c
    
    kappa_rep = kron(kappa,ones(J,1));
    theta_rep = kron(theta,ones(J,1));
    
    z(S*J^2+1:S*J*(J+1)) = P - kappa_rep.*(THETA_sum).^(-1./theta_rep); % P
    
    Y = sum(reshape(kron(X,ones(J,1)).*pi,J,J,S),2);
    Y =Y(:);
    Y_rep = repmat(Y,S,1);
    Y_gamma = sum(reshape(Y_rep.*gamma,J,S,S),3);
    Y_gamma = Y_gamma(:);
    z_temp  = X-(Y_gamma+alpha.*repmat(w.*L+D,S,1));
    z_temp(end)=[];
    z(S*J*(J+1)+1:S*J*(J+2)-1) = z_temp;
    
    X_sum = sum(reshape(X,J,S),2);
    Y_sum = sum(reshape(Y,J,S),2);
    z(S*J*(J+2):S*J*(J+2)+J-1) = X_sum-(Y_sum+D);
    
end

function z = g(w,cP,par,S,J)
    
    z       = zeros(2*S*J,1);
    T       = par.T;
    tau     = par.tau;
    theta   = par.theta;
    gamma   = par.gamma;
    gamma_sum = par.gamma_sum;
    upsilon = par.upsilon;
    %alpha   = par.alpha;
    kappa   = par.kappa;
    w           = [ 1; w];
    c           = exp(cP(1:S*J));
    P           = exp(cP(S*J+1:2*S*J));
    
    THETA       = repmat(T,J,1).*((repmat(c,J,1).*tau).^(-kron(theta,ones(J^2,1))));
    THETA_sum   = sum(reshape(THETA,J,S*J),1)';
    P_rep       = repmat(P,S,1);
    P_gamma     = P_rep.^gamma;
    P_prod      = prod(reshape(P_gamma,J,S,S),2);
    P_prod      = P_prod(:);
    
    z(1:S*J) = c - upsilon.*(repmat(w,S,1).^(1-gamma_sum)).*P_prod; % c
    
    kappa_rep = kron(kappa,ones(J,1));
    theta_rep = kron(theta,ones(J,1));
    
    z(S*J+1:2*S*J) = P - kappa_rep.*(THETA_sum).^(-1./theta_rep); % P
    
   
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
    w = exp(wX(1:J-1));
    X = exp(wX(J:S*J+J-1));
    options = optimoptions('fsolve','StepTolerance',1e-8,...
    'MaxFunctionEvaluations',1e+6,'MaxIterations',5e+3);
    cP0     = [randn(S*J,1); randn(S*J,1)*0.1];
    z_sol   = fsolve(@(cP) g(w,cP,par,S,J), cP0, options);
%     gg = @(cP) g(w,cP,par,S,J);
%     [z_sol,~,flag]   = broyden(gg, cP0);
    c_temp = exp(z_sol(1:S*J));
    %P_temp = exp(z_sol(S*J+1:end));
    THETA       = repmat(T,J,1).*((repmat(c_temp,J,1).*tau).^(-kron(theta,ones(J^2,1))));
    THETA_sum   = sum(reshape(THETA,J,S*J),1)';
    pi_temp =  THETA./kron(THETA_sum,ones(J,1));

    Y = sum(reshape(kron(X,ones(J,1)).*pi_temp,J,J,S),2);
    Y =Y(:);
    Y_rep = repmat(Y,S,1);
    Y_gamma = sum(reshape(Y_rep.*gamma,J,S,S),3);
    Y_gamma = Y_gamma(:);
    z_temp = X-(Y_gamma+alpha.*repmat([1;w].*L+D,S,1));
    z(1:S*J-1,1) = z_temp(2:end);
    X_sum = sum(reshape(X,J,S),2);
    Y_sum = sum(reshape(Y,J,S),2);
    z_temp = X_sum-(Y_sum+D);
    z(S*J:S*J+J-1,1) = z_temp;
end
