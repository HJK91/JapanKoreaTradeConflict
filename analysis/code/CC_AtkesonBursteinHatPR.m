function [P_hat,pi_hat] = CC_AtkesonBursteinHatPR(prod_shr, firm_shr,tau)
% w: w_hat percentage wage change vector (Jx1 vector)
% S_prev: vector of market share of HF and PR firms. (NKx1 vector)
% tau: trade barrier percentage change between Japan and Korea (JxJ matrix).
% Gotta solve for P_hat, which depends on variety level prices changes.

% parameters
J = 4;
K = size(firm_shr,1); % # of varieties in each country
N = 20; % # of goods in each country
eta = 1.05; % substitution parameter across goods (close to one)
rho = 10;   % substitution parameter across varieties (large)

% Unit cost percentage change = (wage change) X (trade cost change).
% let (NKJ)XJ matrix C denote the percentage change of cost in each
% country. Column is the importing country.
firm_prevrep = repmat(firm_shr,1,1,N);
prod_prevrep = repmat(prod_shr,1,1,N);
C = repmat(kron(tau,ones(K/J,1)),1,1,N); % Take the wage change first.
pi_hat = zeros(J,J);
%% Now, solve for hat!
options = optimset('Display','none');
p_after = zeros(K,J,N);
%     gap = 1;
%     while gap>1e-6
%         p_after = p_vrt_after(C_j,S_prev,p_before,rho,eta);
%         gap = norm(p_after(:)-p_before(:));
%         p_before = p_after;
%         disp(gap)
%     end
for j = 1: J
    p_ini = ones(K,1,N);
    firm_ini = firm_prevrep(:,j,:);
    prod_ini = prod_prevrep(:,j,:);
    C_ini = C(:,j,:);
    S_prevpositive = firm_shr(:,j);
    temp_positive = 1:length(S_prevpositive);
    S_prevpositive = temp_positive(S_prevpositive>0);
    p_ini = p_ini(S_prevpositive,1,:);
    firm_ini = firm_ini(S_prevpositive,1,:); % prod_shr, firm_shr
    prod_ini = prod_ini(S_prevpositive,1,:);
    C_ini = C_ini(S_prevpositive,1,:);
    temp_price = fsolve(@(p_before) p_vrt_after(C_ini,firm_ini, prod_ini,p_before,rho,eta,N), p_ini,options);
    p_after(S_prevpositive,j,:) = temp_price;
    temp_hat1 = S_k(temp_price,prod_ini,rho,N);
    temp_hat2 = zeros(K,1);
    temp_hat2(S_prevpositive) = temp_hat1(:,:,1);
    s_hat = reshape(temp_hat2,K/J,J);
    temp_s1 = reshape(prod_shr(:,j,1),K/J,J);
    temp_s2 = temp_s1./sum(temp_s1,1);
    pi_hat(:,j) = sum(temp_s2.*s_hat,1)';
end
pi_hat(isnan(pi_hat)) = 1;
P_hat = P_agg(prod_prevrep,p_after,rho,eta,N,J); 

%% functions for hat variables
function z = P_omega(p_vrt_before, prod_shr,rho,N)
    p_vrt_before1 = reshape(p_vrt_before(p_vrt_before>0),[],1,N);
    S_prev1 = reshape(prod_shr(p_vrt_before>0),[],1,N);
    z = sum(S_prev1.*p_vrt_before1.^(1-rho),1).^(1/(1-rho));
end
function z=S_k(p_vrt_before,prod_shr,rho,N) % S_hat
    z = (p_vrt_before./P_omega(p_vrt_before,prod_shr,rho,N)).^(1-rho);
end

function z = epsilon_hat(p_vrt_before,prod_shr, firm_shr,rho,eta,N)
    temp = 1./((1/eta)*firm_shr.*S_k(p_vrt_before,prod_shr,rho,N) ...
                +(1/rho)*(1-firm_shr.*S_k(p_vrt_before,prod_shr,rho,N)));
    z = temp./epsilon_prev(firm_shr,rho,eta);
end

function z = epsilon_prev(firm_shr,rho,eta)
    z = 1./((1/eta)*firm_shr+(1/rho)*(1-firm_shr));
end

function z = p_vrt_after(C_matrix,prod_shr,firm_shr,p_vrt_before,rho,eta,N)
   z= (epsilon_hat(p_vrt_before,prod_shr,firm_shr,rho,eta,N).*(epsilon_prev(firm_shr,rho,eta)-1)./ ...
        (epsilon_hat(p_vrt_before,prod_shr,firm_shr,rho,eta,N).*epsilon_prev(firm_shr,rho,eta)-1)).*C_matrix;
    z = z-p_vrt_before;
    z(isnan(z)) = 0;
end



function z = P_agg(prod_shr,p_vrt_before,rho,eta,N,J)
    z = zeros(J,1);
    for k = 1 :J
        temp = P_omega(p_vrt_before(:,k,:), prod_shr(:,k,:),rho,N);
        temp1 = temp(temp>0 & temp<1e+10);
        temp2 = (1/length(temp1))*temp1.^(1-eta);
        z(k) = sum(real(temp2),3).^(1/(1-eta));
    end
    
end


end