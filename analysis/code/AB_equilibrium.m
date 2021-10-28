

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

options     = optimoptions('fsolve','Display','iter');
J= par_fixed.J;
S= par_fixed.S;

w_0 = sort(randn(J,1),'descend');
w_sol   = fsolve(@(w) mktclearing(w,tau,par_fixed,data), w_0, options);
    % Solve for equilibrium wage, using trade balance condition
solution.w = exp(w_sol(:,1))./exp(w_sol(1,1)); % Save the result in structure 'solution'
% Given the solution, obtain (c, P, pi, X) using while loop.
cP_prev     = [randn(S*J,1); randn(S*J,1)]; % This vector is converted to exp(cP_prev) in function g.
gap         = 1;
tol         = 1e-8; 
while gap>tol
    cP_new      = log(cPeq(solution.w(:,1),cP_prev,tau,par_fixed)); % cPeq(w,cP,tau,par_fixed)
    gap         = norm(exp(cP_prev)-exp(cP_new));
    cP_prev     = cP_new;
end
solution.c  = exp(reshape(cP_new(1:S*J),J,S));
solution.P  = exp(reshape(cP_new(S*J+1:end),J,S));
solution.pi = reshape(pieq(cP_prev,tau,par_fixed).*par_fixed.mypi(:),J,J,S); %pieq(cP,tau,par_fixed)
% Finally, using market clearing conditions compute expenditure X.
gap         = 1;
X_prev      = exp(randn(S*J,1));
while gap > 1e-8
    Y_prev  = squeeze(sum(lastidx(X_prev,S,J,J).*solution.pi,2));
    X_new   = Xeq(solution.w,solution.P,Y_prev,par_fixed,data);  % Xeq(w,P,Y_prev,par_fixed,data)
    gap     = norm(X_prev-X_new);
    X_prev  = X_new;
end
solution.X = reshape(X_prev,J,S);
solution.Y = Y_prev;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
function z = cPeq(w,cP,tau,par_fixed)
% For simplicity, let chip and PR sectors take last three indices in s.
% index: (1,...S = electronic, S+1 = semiconductor, S+2 PR S+3 HF).
    S = par_fixed.S;
    J = par_fixed.J;    
    theta       = par_fixed.theta;
    delta       = par_fixed.delta;
    mypi        = par_fixed.mypi;
    mygamma     = par_fixed.mygamma;
    mygamma_sum = par_fixed.mygamma_sum;
    z           = zeros(J,2*S);
    c       = reshape(exp(cP(1:S*J)),J,S);           % SJx1 vector c_hat
    P       = reshape(exp(cP(S*J+1:2*S*J)),J,S);     % SJx1 vector P_hat
    c_temp = mididx(c,S,J,J);   % expand the index to match SJ^2x1 vector tau.
    P_rep       = repmat(P(:),S,1);  % expand the index to match (S^2)J by 1 vector mygamma.
    P_mygamma   = reshape(P_rep,J,S,S).^mygamma;
    P_prod      = prod(reshape(P_mygamma,J,S,S),2);   % get products across absorbed sectors
    P_prod      = reshape(P_prod,J,S);   
    % This is the producs of price indices across usual sectors.
    z(:,1:S)    = (repmat(w,1,S).^(1-mygamma_sum)).*P_prod; % get c for usual sectors
    z(:,S-2)    = (w.^(1-mygamma(:,S-1,S-2))).*(P(:,S-1).*delta+P(:,S).*(1-delta)).^(mygamma(:,S-1,S-2));  % delta: cost share of HF and PR composition.
    
    theta_rep   = kron(theta,ones(J,1));
    THETA       = mypi(:).*((c_temp(:).*reshape(tau,S*J^2,1)).^(-kron(theta,ones(J^2,1))));  % numerators of pi.
    THETA_sum   = sum(reshape(THETA,J,S*J),1)';        
    z(:,S+1:end) = reshape((THETA_sum).^(-1./theta_rep),J,S); % get P
    tau_PR = tau(2,3,S-1);
    tau_HF = tau(2,3,S);
    z(:,end-1) = AB_AtkesonBurstein(w,tau_PR);
    z(:,end) = AB_AtkesonBurstein(w,tau_HF);
    z       = reshape(z,2*J*S,1);
end

function z = Xeq(w,P,Y_prev,par_fixed,data)
    S = par_fixed.S;
    J = par_fixed.J;
    P = exp(reshape(P,J,S));
    D = data.D;
    wL = par_fixed.wL;
    
    mygamma   = reshape(par_fixed.mygamma,J,S,S);
    mygamma_sum   = reshape(par_fixed.mygamma_sum,J,S);
    delta       = par_fixed.delta;
    alpha   = reshape(par_fixed.alpha,J,S);
    % production amount
    % S: HF industry
    % S-1: PR industry
    % S-2: memory chip industry
    % S-3: electronic device industry
    Y_prev      = reshape(Y_prev,J,S);
    Y_temp      = mididx(Y_prev,S,S,J);
    Y_mygamma   = sum(Y_temp(:,:,1:end-3).*mygamma(:,:,1:end-3),3);
    Y_mygamma   = reshape(Y_mygamma(:),J,S);
%     z            = Y_mygamma+alpha.*repmat(w.*wL+sum(D,2),1,S);
    z            = Y_mygamma+alpha.*repmat(w.*wL+D,1,S);
    z(:,S-2)    = mygamma(:,S-2,S-3).*Y_prev(:,S-3); % expenditure of chip = cost share of electronic devices
    z(:,S-1)    = mygamma_sum(:,S-2).*((P(:,S-1).*delta)./(P(:,S-1).*delta+P(:,S).*(1-delta))).*Y_prev(:,S-2); % expenditure of PR = cost share of chips
    z(:,S)      = mygamma_sum(:,S-2).*((P(:,S).*(1-delta))./(P(:,S-1).*delta+P(:,S).*(1-delta))).*Y_prev(:,S-2);
    z = z(:);
end

function z = pieq(cP,tau,par_fixed)
    S       = par_fixed.S;
    J       = par_fixed.J;    
    theta   = par_fixed.theta;
    c = exp(mididx(cP(1:S*J),S,J,J));
    P = exp(lastidx(cP(S*J+1:end),S,J,J));
    z = (c(:).*tau(:)./P(:)).^(-kron(theta,ones(J^2,1))) ;
end

function z = mktclearing(w,tau,par_fixed,data)
% w: J-1 vector (first country's wage is numeraire)
% X: SxJ vector
    S = par_fixed.S;
    J = par_fixed.J;
    
    D = data.D;
    w       = exp(w(1:J));
    w       = w/w(1);       % Normalize the wage vector
    X       = data.X_ttl(:);
    cP_prev = randn(2*S*J,1); % LOG of cost and price: this can be negative!
    gap     = 1;
    tol     = 1e-8;
    while gap>tol
        cP_new  = log(cPeq(w,cP_prev,tau,par_fixed));
        gap     = norm(exp(cP_prev)-exp(cP_new));
        cP_prev = cP_new;
    end
    pi_temp     = reshape(pieq(cP_prev,tau,par_fixed),J,J,S).*par_fixed.mypi; % new trade shares after change 
    gap         = 1;
    X_prev=X;
    P_prev = cP_prev(S*J+1:end);
     while gap > 1e-8
        Y_prev  = sum(lastidx(X_prev,S,J,J).*pi_temp,2);
        X_new   = Xeq(w,P_prev,Y_prev,par_fixed,data);
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
    % idx3: last index in notation, first index in matrix
    z = reshape(repmat(reshape(x,idx3,idx1),idx2,1),idx3,idx2,idx1);
end


function z = lastidx(x,idx1,idx2,idx3)
    % expand the last index
    % x is (idx1*idx2) by 1 vector
    % idx3: last index in notation, first index in matrix
    z = reshape(repmat(reshape(x,1,idx1*idx2),idx3,1),idx3,idx2,idx1);
end