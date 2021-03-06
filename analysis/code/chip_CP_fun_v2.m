 function z = chip_CP_fun_v2(par, par_fixed,data)

% Let's extend Caliendo and Parro (2015) to incorporate chip market
% environment. To consider the GVC through chip industry, consider S+3
% industries. There are S 'ususal' industries including electronic device.
% Electronic device sector source high and low quality chips and combine
% them in C-D fashion. These differentiated chips are classified as two
% distinguished sectors other than usual industries. In particular, high
% quality chips requires a mixture of high quality PR and chemical inputs
% in its recipe. Low quality chips, on the other hand, only requires
% chemical inputs.

% par.S;      % Number of industries    % Accommodate H-chips and PR sectors.
% par.J;      % Number of countries

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Data and Equilibrium solution %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% par0 = [par.T; par.tau; par.sigma];

options     = optimoptions('fsolve','Display','off');
J= par_fixed.J;
S= par_fixed.S;
par_struc.T = par(1:S*J);
par_struc.tau = par(S*J+1:S*J+J^2*S);
par_struc.sigma = par(S*J+J^2*S+1:end);
w_0 = randn(J,1);
w_sol   = fsolve(@(w) mktclearing(w,par_struc,par_fixed), w_0, options);
    % Solve for equilibrium wage, using trade balance condition
solution.w = exp(w_sol(:,1))./exp(w_sol(1,1)); % Save the result in structure 'solution'
% Given the solution, obtain (c, P, pi, X) using while loop.
cP_prev     = [randn(S*J,1); randn(S*J,1)]; % This vector is converted to exp(cP_prev) in function g.
gap         = 1;
tol         = 1e-8; 
while gap>tol
    cP_new      = log(cPeq(solution.w(:,1),cP_prev,par_struc,par_fixed));
    gap         = norm(exp(cP_prev)-exp(cP_new));
    cP_prev     = cP_new;
end
solution.c  = exp(cP_new(1:S*J));
solution.P  = exp(cP_new(S*J+1:end));
% Using solution.c, obtain the bilateral trade shares.
solution.pi = pieq(cP_prev,par_struc,par_fixed); 
% Finally, using market clearing conditions compute expenditure X.
gap         = 1;
X_prev      = exp(randn(S*J,1));
while gap > 1e-8
    Y_prev  = sum(reshape(lastidx(X_prev,S,J,J).*solution.pi,J,J,S),2);
    X_new   = Xeq(solution.w,Y_prev,par_fixed,data); 
    gap     = norm(X_prev-X_new);
    X_prev  = X_new;
end
solution.X = X_new;

%% Map the model into the data

% From data, we get Z, Y, F, and VA.

% M = squeeze(sum(reshape(WIOD_Z,J,J,S,S),3));
% WIOD_X = squeeze(sum(M.*reshape(par.tau,J,J,S),1));
% gap_pi  = norm(reshape(solution.pi-WIOD_X./sum(WIOD_X),J^2*S,1));
model_Y = sum(reshape(lastidx(solution.X,S,J,J).*solution.pi,J,J,S),2);
temp  = reshape(mididx(model_Y(:),S,S,J),J,S,S);
Y_gamma = sum(temp(:,:,1:end-2).*par_fixed.mygamma(:,:,1:end-2),3);
model_F = solution.X-Y_gamma(:);
temp1 = lastidx(model_Y(:)-model_F,S,J,J).*solution.pi;
temp1 = mididx(temp1,J^2,S,S);
temp2 = lastidx(par_fixed.mygamma(:),S^2,J,J);
model_Z = reshape(temp1.*temp2,J,J,S,S);
model_VA= model_Y(:)-reshape((sum(model_Z,[1,3])),J*S,1);

z=norm([data.Z(:)-model_Z(:); data.F(:)-model_F(:); ...
    data.Y(:)-model_Y(:); data.VA(:)-model_VA(:)]);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
function z = cPeq(w,cP,par,par_fixed)
% For simplicity, let chip and PR sectors take last three indices in s.
% index: (1,...S-4 = chemical sector, S-3 = electronic device sector, L_chip, H_chip, PR).
    S = par_fixed.S;
    J = par_fixed.J;    
    z           = zeros(J,2*S);
    T           = reshape(par.T,J,S); 
    tau         = reshape(par.tau,J,J,S);
    theta       = par_fixed.theta;
    mygamma       = reshape(par_fixed.mygamma,J,S,S);
    mygamma_sum   = reshape(par_fixed.mygamma_sum,J,S);
    upsilon     = reshape(par_fixed.upsilon,J,S);
    kappa   = par_fixed.kappa;
    c       = reshape(exp(cP(1:S*J)),J,S);           % SJx1 vector
    P       = reshape(exp(cP(S*J+1:2*S*J)),J,S);     % SJx1 vector
    
    c_temp = mididx(c,S,J,J);   % expand the index to match SJ^2x1 vector tau.
    T_temp = mididx(T,S,J,J);
    
    THETA       = T_temp.*((c_temp.*reshape(tau,S*J^2,1)).^(-kron(theta,ones(J^2,1))));  % numerators of pi.
    THETA_sum   = sum(reshape(THETA,J,S*J),1)';                         % denominators of pi.
    P_rep       = mididx(P,S,S,J);  % expand the index to match (S^2)J by 1 vector mygamma.
    P_mygamma     = reshape(P_rep,J,S,S).^mygamma;
    P_prod      = reshape(P_mygamma,J,S,S);   % get products across absorbed sectors
    P_prod      = reshape(prod(P_prod(:,1:S-3,:),2),J,S);   
%     P_prod_nochem   = reshape(prod(P_mygamma(:,(1:S-4),:),2).*P_mygamma(:,S-2,:),J,S);  
    % This is the producs of price indices across usual sectors.
    P_prod      = reshape(P_prod(:),J,S); % SxJ by 1 vector
    z(:,1:S-3)  = upsilon(:,1:S-3).*(repmat(w,1,S-3).^(1-mygamma_sum(:,1:S-3))).*P_prod(:,1:S-3); % get c for usual sectors
    z(:,S-2)    = upsilon(:,S-2).*(w.^(1-mygamma_sum(:,S-2)-mygamma(:,S-1,S-2))).*P_prod(:,S-2) ...
        .*P(:,S-1).^(mygamma(:,S-1,S-2));  % contribution of unit price of high-tier chips
    z(:,S-1)    = upsilon(:,S-1).*(w.^(1-mygamma_sum(:,S-1))) ... 
        .*(P(:,S).^(mygamma(:,S,S-1))); % contribution of unit price of PR-Chemical bundle
    z(:,S)      = w; % PR production cost
    kappa_rep   = kron(kappa,ones(J,1));
    theta_rep   = kron(theta,ones(J,1));    
    z(J*S+1:end) = reshape(kappa_rep.*(THETA_sum).^(-1./theta_rep),J,S); % get P
    z       = reshape(z,2*J*S,1);
end

function z = Xeq(w,Y_prev,par_fixed,data)
    S = par_fixed.S;
    J = par_fixed.J;
%     P       = reshape(P,J,S);
    L = data.L;
    D = data.D;
    mygamma   = reshape(par_fixed.mygamma,J,S,S);
    alpha   = reshape(par_fixed.alpha,J,S);
    % production amount
    % S: PR sector
    % S-1: chips sector
    % S-2: electronic device sector
    % S-3: chemical sector
    Y_prev = reshape(Y_prev,J,S);
    Y_temp  = reshape(mididx(Y_prev(:),S,S,J),J,S,S);
    Y_mygamma = sum(Y_temp(:,:,1:end-2).*mygamma(:,:,1:end-2),3);
    Y_mygamma = reshape(Y_mygamma(:),J,S);
    z       = Y_mygamma+alpha.*repmat(w.*L+sum(D,2),1,S);
    z(:,S-1)=mygamma(:,S-2,S-1).*Y_prev(:,S-2);
    z(:,S)=mygamma(:,S-1,S).*Y_prev(:,S-1);
    z = z(:);
end

function z = pieq(cP,par,par_fixed)
    S = par_fixed.S;
    J = par_fixed.J;    
    T = par.T;
    tau = par.tau;
    theta = par_fixed.theta;
    c_temp      = mididx(exp(cP(1:S*J)),S,J,J);
    T_temp      = mididx(T,S,J,J);
    THETA       = T_temp.*((c_temp.*tau).^(-kron(theta,ones(J^2,1))));
    THETA_sum   = sum(reshape(THETA,J,S*J),1)';
    z     = THETA./lastidx(THETA_sum,S,J,J);
end

function z = mktclearing(w,par,par_fixed)
% w: J-1 vector (first country's wage is numeraire)
% X: SxJ vector
    S = par_fixed.S;
    J = par_fixed.J;
    D = data.D;
    w       = exp(w(1:J));
    w       = w/w(1);       % Normalize the wage vector
    X       = exp(randn(S*J,1));
    cP_prev = randn(2*S*J,1); % LOG of cost and price: this can be negative!
    gap     = 1;
    tol     = 1e-8;
    while gap>tol
        cP_new  = log(cPeq(w,cP_prev,par,par_fixed));
        gap     = norm(exp(cP_prev)-exp(cP_new));
        cP_prev = cP_new;
    end
    pi_temp     = pieq(cP_prev,par,par_fixed);
%     p_temp      = cP_prev(S*J+1:end);
    gap         = 1;
    X_prev=X;
     while gap > 1e-8
        Y_prev  = sum(reshape(lastidx(X_prev,S,J,J).*pi_temp,J,J,S),2);
        X_new   = Xeq(w,Y_prev,par_fixed,data);
        gap     = norm(X_prev-X_new);
        X_prev  = X_new;
     end
     X_sum  = sum(reshape(X_prev,J,S),2);
     Y_sum  = sum(reshape(Y_prev,J,S),2);
     z_temp = X_sum-(Y_sum+sum(D,2));
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
end