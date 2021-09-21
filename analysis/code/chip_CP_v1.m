% Let's extend Caliendo and Parro (2015) to incorporate chip market
% environment. To consider the GVC through chip industry, consider S+3
% industries. There are S 'ususal' industries including electronic device.
% Electronic device sector source high and low quality chips and combine
% them in C-D fashion. These differentiated chips are classified as two
% distinguished sectors other than usual industries. In particular, high
% quality chips requires a mixture of high quality PR and chemical inputs
% in its recipe. Low quality chips, on the other hand, only requires
% chemical inputs.

clc
clear
S = 3;      % Number of industries
S = S+3;    % Accommodate H-chips, L-chips and H-RP sectors.
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
% For instance, a (SJ)^2 by 1 vector W^{rs}_{ni} goes as
% [W^{11}_{11}; W^{11}_{12}; ... ; W^{11}_{1J};
%   W^{11}_{21};W^{11}_{22},...; W^{11}_{JJ}; W^{12}_{11}; W^{12}_{12};...; 
%   W^{1J}_{JJ};...; W^{21}_{11}; ...; W^{JJ}_{JJ};].

par.T       = rand(J*S,1)*0.2+1;        % Technology parameter T_j^s
par.tau     = rand(S*J^2,1)+1;          % Trade cost parameters tau_ni^s
par.theta   = 5+rand(S,1)*0.1;          % CA parameter theta
par.sigma   = par.theta+1-rand(S,1);    % Substitution parameter sigma^s
par.gamma           = sort(rand(J,S+1,S),2);    % Industry share parameter gamma^sr_j
par.gamma(:,1,:)    = 0;                % Make gamma to satisfy restriction
par.gamma           = diff(sort(par.gamma,2),1,2);  % 0< Sum of gamma across r <1
par.gamma           = par.gamma(:);                 
par.gamma_sum       = sum(reshape(par.gamma,J,S,S),2); % 1-labor share = gamma_sum
par.gamma_sum       = par.gamma_sum(:);
par.delta           = rand(1);          % share of high-end chip within chip sector for electronic industry. 
par.upsilon = prod(reshape(par.gamma.^(-par.gamma),J,S,S),2); 
par.upsilon = par.upsilon(:).*(1-par.gamma_sum).^(-(1-par.gamma_sum));  % As in Caliendo and Parro (2015)
par.alpha   = diff([zeros(1,J); sort(rand(S-1,J),1) ; ones(1,J)])';     % Restrictio: sum of alpha_j^s across s is 1
par.alpha   = par.alpha(:);
par.kappa   = gamma(1+(1-par.sigma)./par.theta).^(1./(1-par.sigma));    % As in Caliendo and Parro (2015)
par.K       = 2; % the ratio of 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Data and Equilibrium solution %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L_0         = ones(J,1);        % Jx1 data
D_0         = zeros(J,1);       % Jx1 data
N           = 5;                % Number of initial points to try
w_sol       = zeros(J,N);       % Wage vectors are contained here
fval        = zeros(J,N);       % Function values that are meant to be close to zero
flag        = zeros(N,1);       % Flags
options     = optimoptions('fsolve','Display','iter');
for i = 1: N
    w_0 = randn(J,1);
    [w_sol(:,i), fval(:,i), flag(i)]   = fsolve(@(w) mktclearing(w,par,S,J,L_0,D_0), w_0, options);
    % Solve for equilibrium wage, using trade balance condition
end
%%
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
    X_new   = Xeq(solution.w,solution.P,Y_prev,par,S,J,L_0,D_0);
    gap     = norm(X_prev-X_new);
    X_prev  = X_new;
end
solution.X = X_new;
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
function z = cPeq(w,cP,par,S,J)
% For simplicity, let chip and PR sectors take last three indices in s.
% index: (1,...S-4 = chemical sector, S-3 = electronic device sector, L_chip, H_chip, PR).
    z           = zeros(J,2*S);
    T           = reshape(par.T,J,S); 
    tau         = reshape(par.tau,J,J,S);
    theta       = par.theta;
    gamma       = reshape(par.gamma,J,S,S);
    gamma_sum   = reshape(par.gamma_sum,J,S);
    upsilon     = reshape(par.upsilon,J,S);
    kappa   = par.kappa;
    c       = reshape(exp(cP(1:S*J)),J,S);           % SJx1 vector
    P       = reshape(exp(cP(S*J+1:2*S*J)),J,S);     % SJx1 vector
    
    c_temp = mididx(c,S,J,J);   % expand the index to match SJ^2x1 vector tau.
    T_temp = mididx(T,S,J,J);
    
    THETA       = T_temp.*((c_temp.*reshape(tau,S*J^2,1)).^(-kron(theta,ones(J^2,1))));  % numerators of pi.
    THETA_sum   = sum(reshape(THETA,J,S*J),1)';                         % denominators of pi.
    P_rep       = mididx(P,S,S,J);  % expand the index to match (S^2)J by 1 vector gamma.
    P_gamma     = reshape(P_rep,J,S,S).^gamma;
    P_prod      = reshape(P_gamma,J,S,S);   % get products across absorbed sectors
    P_prod      = reshape(prod(P_prod(:,1:S-3,:),2),J,S);   
    P_prod_nochem   = reshape(prod(P_gamma(:,(1:S-5),:),2).*P_gamma(:,S-3,:),J,S);  
    % This is the producs of price indices across usual sectors.
    P_prod      = reshape(P_prod(:),J,S); % SxJ by 1 vector
    z(:,1:S-4)  = upsilon(:,S-4).*(repmat(w,1,S-4).^(1-gamma_sum(:,S-4))).*P_prod(:,1:S-4); % get c for usual sectors
    z(:,S-3)    = upsilon(:,S-3).*(w.^(1-gamma_sum(:,S-3)-gamma(:,S-2,S-3))).*P_prod(:,S-3) ...
        .*P(:,S-2).^(gamma(:,S-2,S-3)*par.delta) ...  % contribution of unit price of high-tier chips
        .*P(:,S-1).^(gamma(:,S-2,S-3)*(1-par.delta)); % contribution of unit price of low-tier chips (same gamma as in high-end chips)
    z(:,S-2)    = upsilon(:,S-2).*(w.^(1-gamma_sum(:,S-2))).*P_prod_nochem(:,S-2) ... 
        .*(P(:,S-4)+par.K*P(:,S)); % contribution of unit price of PR-Chemical bundle
    z(:,S-1)  = upsilon(:,S-1).*(w.^(1-gamma_sum(:,S-1))).*P_prod(:,S-1); % get c for low-tier chip sectors
    z(:,S)      = w;
    kappa_rep   = kron(kappa,ones(J,1));
    theta_rep   = kron(theta,ones(J,1));    
    z(J*S+1:end) = reshape(kappa_rep.*(THETA_sum).^(-1./theta_rep),J,S); % get P
    z       = reshape(z,2*J*S,1);
end

function z = Xeq(w,P,Y_prev,par,S,J,L,D)
    P       = reshape(P,J,S);
    gamma   = reshape(par.gamma,J,S,S);
    alpha   = reshape(par.alpha,J,S);
    % production amount
    % S: PR sector
    % S-1: low-tier chip sector
    % S-2: high-tier chips sector
    % S-3: electronic device sector
    % S-5:
    Y_prev = reshape(Y_prev,J,S);
    Y_temp  = reshape(mididx(Y_prev(:),S,S,J),J,S,S);
    Y_gamma = sum(Y_temp.*gamma,3);
    Y_gamma = reshape(Y_gamma(:),J,S);
    z       = Y_gamma+alpha.*repmat(w.*L+D,1,S);
    z(:,S-4)= z(:,S-4)+gamma(:,S-2,S-4) ...
        +(par.K*P(:,end))./(par.K*P(:,end)+P(:,S-4)).*Y_prev(:,S-2);
    z(:,S-2)=gamma(:,S-2,S-3).*par.delta.*Y_prev(:,S-2);
    z(:,S-1)=gamma(:,S-1,S-3).*(1-par.delta).*Y_prev(:,S-1);
    z(:,S)=gamma(:,S-4,S-2).*(par.K*P(:,end))./(par.K*P(:,end)+P(:,S-4)).*Y_prev(:,S-2);
    z = z(:);
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
    p_temp      = cP_prev(S*J+1:end);
    gap         = 1;
    X_prev=X;
     while gap > 1e-8
        Y_prev  = sum(reshape(lastidx(X_prev,S,J,J).*pi_temp,J,J,S),2);
        X_new   = Xeq(w,p_temp,Y_prev,par,S,J,L,D);
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
