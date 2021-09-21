clc
clear
S = 6;      % Number of industries
J = 5;      % Number of countries

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Indexing convention
% s   |  r  |   n   |  i  | where
% s: importing industry
% r: exporting industry
% n: importing country
% i: exporting country
% For instance, a (SJ)^2 vector W^{rs}_{ni} goes as
% [W^{11}_{11}; W^{11}_{12}; ... ; W^{11}_{1J};
%   W^{11}_{21};W^{11}_{22},...; W^{11}_{JJ}; W^{12}_{11}; W^{12}_{12};...; 
%   W^{1J}_{JJ};...; W^{21}_{11}; ...; W^{JJ}_{JJ};].

par.T       = rand(J*S,1)*0.2+1;        % Technology parameter T_j^s
par.tau     = rand(S*J^2,1)+1;          % Trade cost parameters tau_ni^s
par.theta   = 5+rand(S,1)*0.1;          % CA parameter theta
par.sigma   = par.theta+1-rand(S,1);    % Substitution parameter sigma^s
par.gamma           = sort(rand(J,S+1,S),2);    % Industry share parameter gamma^sr_j
par.gamma(:,1,:)    = 0;                 % Make gamma to satisfy restriction
par.gamma           = diff(sort(par.gamma,2),1,2);  % 0< Sum of gamma across r <1
par.gamma           = par.gamma(:);                 
par.gamma_sum       = sum(reshape(par.gamma,J,S,S),2); % 1-labor share = gamma_sum
par.gamma_sum       = par.gamma_sum(:);
par.upsilon = prod(reshape(par.gamma.^(-par.gamma),J,S,S),2); 
par.upsilon = par.upsilon(:).*(1-par.gamma_sum).^(-(1-par.gamma_sum));  % As in Caliendo and Parro (2015)
par.alpha   = diff([zeros(1,J); sort(rand(S-1,J),1) ; ones(1,J)])';     % Restrictio: sum of alpha_j^s across s is 1
par.alpha   = par.alpha(:);
par.kappa   = gamma(1+(1-par.sigma)./par.theta).^(1./(1-par.sigma));    % As in Caliendo and Parro (2015)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Data and Equilibrium solution %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L_0         = ones(J,1);        % Jx1 data
D_0         = zeros(J,1);       % Jx1 data
options     = optimoptions('fsolve','Display','iter');
N           = 5;                % Number of initial points to try
w_sol = zeros(J,N);             % Wage vectors are contained here
fval = zeros(J,N);              % Function values that are meant to be close to zero
flag = zeros(N,1);              % Flags
for i = 1: N
    w_0 = randn(J,1);
    [w_sol(:,i), fval(:,i), flag(i)]   = fsolve(@(w) mktclearing(w,par,S,J,L_0,D_0), w_0, options);
    % Solve for equilibrium wage, using trade balance condition
end
solution.w = exp(w_sol(:,1))./exp(w_sol(1,1)); % Save the result in structure 'solution'
% Given the solution, obtain (c, P, pi, X) using while loop.
cP_prev     = [randn(S*J,1); randn(S*J,1)]; % This vector is converted to exp(cP_prev) in function g.
gap         = 1;
tol         = 1e-8;
while gap>tol
    cP_new      = log(cPeq(solution.w(:,1),cP_prev,par,S,J));
    gap         = norm(exp(cP_prev)-exp(cP_new));
    cP_prev     = cP_new;
end
solution.c  = exp(cP_new(1:S*J));
solution.P  = exp(cP_new(S*J+1:end));
% Using solution.c, obtain the bilateral trade shares.
solution.pi = pieq(cP_prev,par,S,J);
% Finally, using market clearing conditions compute expenditure X.
gap         = 1;
X_prev      = exp(randn(S*J,1));
while gap > 1e-8
    Y_prev  = sum(reshape(lastidx(X_prev,S,J,J).*solution.pi,J,J,S),2);
    X_new   = Xeq(solution.w,Y_prev,par,S,J,L_0,D_0);
    gap     = norm(X_prev-X_new);
    X_prev  = X_new;
end
solution.X = X_new;
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
function z = cPeq(w,cP,par,S,J)
    z           = zeros(2*S*J,1);
    T           = par.T;
    tau         = par.tau;
    theta       = par.theta;
    gamma       = par.gamma;
    gamma_sum   = par.gamma_sum;
    upsilon     = par.upsilon;
    kappa   = par.kappa;
    c       = exp(cP(1:S*J));           % SJx1 vector
    P       = exp(cP(S*J+1:2*S*J));     % SJx1 vector
    
    c_temp = mididx(c,S,J,J);   % expand the index to match SJ^2x1 vector tau.
    T_temp = mididx(T,S,J,J);
    
    THETA       = T_temp.*((c_temp.*tau).^(-kron(theta,ones(J^2,1))));  % numerators of pi.
    THETA_sum   = sum(reshape(THETA,J,S*J),1)';                         % denominators of pi.
    P_rep       = mididx(P,S,S,J);  % expand the index to match (S^2)J by 1 vector gamma.
    P_gamma     = P_rep.^gamma;
    P_prod      = prod(reshape(P_gamma,J,S,S),2);   % get products across absorbed sectors
    P_prod      = P_prod(:);
    z(1:S*J)    = upsilon.*(repmat(w,S,1).^(1-gamma_sum)).*P_prod; % get c
    kappa_rep   = kron(kappa,ones(J,1));
    theta_rep   = kron(theta,ones(J,1));    
    z(S*J+1:2*S*J) = kappa_rep.*(THETA_sum).^(-1./theta_rep); % get P
end

function z = Xeq(w,Y_prev,par,S,J,L,D)
    gamma   = par.gamma;
    alpha   = par.alpha;
    % production amount
    Y_temp  = mididx(Y_prev(:),S,S,J);
    Y_gamma = sum(reshape(Y_temp.*gamma,J,S,S),3);
    Y_gamma = Y_gamma(:);
    z       = Y_gamma+alpha.*repmat(w.*L+D,S,1);
end

function z = pieq(cP,par,S,J)
    T = par.T;
    tau = par.tau;
    theta = par.theta;
    c_temp      = mididx(exp(cP(1:S*J)),S,J,J);
    T_temp      = mididx(T,S,J,J);
    THETA       = T_temp.*((c_temp.*tau).^(-kron(theta,ones(J^2,1))));
    THETA_sum   = sum(reshape(THETA,J,S*J),1)';
    z     = THETA./lastidx(THETA_sum,S,J,J);
end
function z = mktclearing(w,par,S,J,L,D)
% w: J-1 vector (first country's wage is numeraire)
% X: SxJ vector
    w       = exp(w(1:J));
    w       = w/w(1);       % Normalize the wage vector
    X       = exp(randn(S*J,1));
    cP_prev = randn(2*S*J,1); % LOG of cost and price: this can be negative!
    gap     = 1;
    tol     = 1e-8;
    while gap>tol
        cP_new  = log(cPeq(w,cP_prev,par,S,J));
        gap     = norm(exp(cP_prev)-exp(cP_new));
        cP_prev = cP_new;
    end
    pi_temp     = pieq(cP_prev,par,S,J);
    gap         = 1;
    X_prev=X;
     while gap > 1e-8
        Y_prev  = sum(reshape(lastidx(X_prev,S,J,J).*pi_temp,J,J,S),2);
        X_new   = Xeq(w,Y_prev,par,S,J,L,D);
        gap     = norm(X_prev-X_new);
        X_prev  = X_new;
     end
     X_sum  = sum(reshape(X_prev,J,S),2);
     Y_sum  = sum(reshape(Y_prev,J,S),2);
     z_temp = X_sum-(Y_sum+D);
     z(1:J) = z_temp;
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
