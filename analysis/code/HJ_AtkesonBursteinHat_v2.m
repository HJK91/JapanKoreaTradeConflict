function z = HJ_AtkesonBursteinHat_v2(w,S_prev,tau)
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
% p_before = rand(K*J,N);
P_sol = zeros(1,J);
p_before0 = rand(K*J,N);
options = optimset('Display','none');
for j = 1:J
    C_j = repmat(C(:,j),1,N);
%     gap = 1;
%     while gap>1e-6
%         p_after = p_vrt_after(C_j,S_prev,p_before,rho,eta);
%         gap = norm(p_after(:)-p_before(:));
%         p_before = p_after;
%         disp(gap)
%     end
    p_after = fsolve(@(p_before) p_vrt_after(C_j,S_prev,p_before,rho,eta), p_before0,options);
    p_before0 = p_after;
    P_sol(:,j) = P_agg(S_prev,p_after,rho,eta,N); 
end
z = P_sol'; % Jx1 vector

%% functions for hat variables

function z=S_k(p_vrt_before,S_prev,rho) % S_hat
    z = (p_vrt_before./P_omega(p_vrt_before,S_prev,rho)).^(1-rho);
end

function z = epsilon_hat(p_vrt_before,S_prev,rho,eta)
    temp = 1./((1/eta)*S_prev.*S_k(p_vrt_before,S_prev,rho) ...
                +(1/rho)*(1-S_prev.*S_k(p_vrt_before,S_prev,rho)));
    z = temp./epsilon_prev(S_prev,rho,eta);
end

function z = epsilon_prev(S_prev,rho,eta)
    z = 1./((1/eta)*S_prev+(1/rho)*(1-S_prev));
end

function z = p_vrt_after(C_matrix,S_prev,p_vrt_before,rho,eta)
   z= (epsilon_hat(p_vrt_before,S_prev,rho,eta).*(epsilon_prev(S_prev,rho,eta)-1)./ ...
        (epsilon_hat(p_vrt_before,S_prev,rho,eta).*epsilon_prev(S_prev,rho,eta)-1)).*C_matrix;
    z = z-p_vrt_before;
end

function z = P_omega(p_vrt_before, S_prev,rho)
    z = sum(S_prev.*p_vrt_before.^(1-rho)).^(1/(1-rho));
end

function z = P_agg(S_prev,p_vrt_before,rho,eta,N)
    temp = 1/N*P_omega(p_vrt_before, S_prev,rho).^(1-eta);
    z = sum(real(temp),2).^(1/(1-eta));
end

% function z = FOC(C_matrix,S_prev,Y,X,rho,eta,K,N,J)
%     Omega = 1/eta-1/rho;
%     Y_vrt = reshape(exp(Y),K*J,N);
%     
%     C_matrix = reshape(C_matrix,K*J,N);
%     temp1 = (Y_vrt.^(-1/rho)).*(Y_omega(Y_vrt,S_prev, rho).^(-Omega)) ...
%         .*P_agg(C_matrix,S_prev,Y_vrt,rho,eta,N).^((eta-1)/eta).*X.^(1/eta);
%     temp2 = (rho-1)/rho-Omega.*S_prev.*S_k(Y_vrt,S_prev,rho);
%     temp3 = (rho-1)/rho-Omega.*S_prev;
%     z = temp1.*temp2./temp3;
%     z = z(:);
% end

end