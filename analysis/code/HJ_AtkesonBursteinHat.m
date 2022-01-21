function z = HJ_AtkesonBursteinHat(w,S_prev,tau)
% w: w_hat percentage wage change vector (Jx1 vector)
% S_prev: vector of market share of HF and PR firms. (NKx1 vector)
% tau: trade barrier percentage change between Japan and Korea (JxJ matrix).
% Gotta solve for P_hat, which depends on variety level prices changes.

% parameters
J = 4;
K = length(S_prev); % # of varieties in each country
N = 20; % # of goods in each country
eta = 1.05; % substitution parameter across goods close to one
rho = 10;   % substitution parameter across varieties

% Unit cost percentage change = (wage change) X (trade cost change).
% let (NKJ)XJ matrix C denote the percentage change of cost in each
% country. Column is the importing country.
S_prev = repmat(S_prev,J,N);
C = kron(w.*tau,ones(K,1)); % Take the wage change first.

%% Now, solve for hat!
options = optimset('display','off');
X = 1;
Y0 = randn(K*N*J,1);
Y_sol = zeros(K*N*J,J);
P_sol = zeros(1,J);
j=1;
while j <= J
    C_j = repmat(C(:,j),1,N);
    [Y_sol(:,j),~,flag] = fsolve(@(Y) FOC(C_j,S_prev,Y,X,rho,eta,K,N,J),Y0,options);
    if flag == 1
        Y0 = Y_sol(:,j);
        Y_sol(:,j) = exp(Y_sol(:,j));
        P_sol(:,j) = P_agg(C_j,S_prev,reshape(Y_sol(:,j),K*J,N),rho,eta,N); 
        j=j+1;
    else
        Y0 = randn(K*N*J,1);
    end
end
z = P_sol'; % Jx1 vector

%% functions for hat variables
function z = Y_omega(Y_vrt,S_prev, rho)
    z = sum(S_prev.*Y_vrt.^((rho-1)/rho)).^(rho/(rho-1));
end

function z=S_k(Y_vrt,S_prev,rho)
    z = (Y_vrt./Y_omega(Y_vrt,S_prev, rho)).^(rho/(rho-1));
end

function z = endo_epsilon(Y_vrt,S_prev,rho,eta)
    mu = (1/eta)*S_prev./((1/eta)*S_prev+(1/rho)*(1-S_prev));
    z = 1./(mu.*S_k(Y_vrt,S_prev,rho)+(1-mu).*(1-S_prev.*S_k(Y_vrt,S_prev,rho))./(1-S_prev));
end

function z = endo_epsilon_prev(Y_vrt,S_prev,rho,eta)
    z = 1./((1/eta)*S_k(Y_vrt,S_prev,rho)+(1/rho)*(1-S_k(Y_vrt,S_prev,rho)));
end

function z = price_vrt(C_matrix,S_prev,Y_vrt,rho,eta)
   z= (endo_epsilon(Y_vrt,S_prev,rho,eta).*(endo_epsilon_prev(Y_vrt,S_prev,rho,eta)-1)./ ...
        (endo_epsilon(Y_vrt,S_prev,rho,eta).*endo_epsilon_prev(Y_vrt,S_prev,rho,eta)-1)).*C_matrix;
end

function z = P_omega(C_matrix,S_prev,Y_vrt,rho,eta)
    z = sum(S_k(Y_vrt,S_prev,rho).*price_vrt(C_matrix,S_prev,Y_vrt,rho,eta).^(1-rho)).^(1/(1-rho));
end

function z = P_agg(C_matrix,S_prev,Y_vrt,rho,eta,N)
    temp = 1/N*P_omega(C_matrix,S_prev,Y_vrt,rho,eta).^(1-eta);
    z = sum(real(temp),2).^(1/(1-eta));
end

function z = FOC(C_matrix,S_prev,Y,X,rho,eta,K,N,J)
    Omega = 1/eta-1/rho;
    Y_vrt = reshape(exp(Y),K*J,N);
    
    C_matrix = reshape(C_matrix,K*J,N);
    temp1 = (Y_vrt.^(-1/rho)).*(Y_omega(Y_vrt,S_prev, rho).^(-Omega)) ...
        .*P_agg(C_matrix,S_prev,Y_vrt,rho,eta,N).^((eta-1)/eta).*X.^(1/eta);
    temp2 = (rho-1)/rho-Omega.*S_prev.*S_k(Y_vrt,S_prev,rho);
    temp3 = (rho-1)/rho-Omega.*S_prev;
    z = temp1.*temp2./temp3-C_matrix;
    z = z(:);
end

end