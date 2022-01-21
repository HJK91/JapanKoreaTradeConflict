function [P_hat,pi_hat] = CC_chiphat(prod_shr, firm_shr,P_input,delta,mygamma)
% S_prev: matrix of market share of chip producers. (KxJ matrix)
% P_input: (Jx2) matrix of price index changes in industry HF and PR in country n.
% Gotta solve for P_hat, which depends on variety level prices changes.

% parameters
J = 4;
K = size(firm_shr,1); % # of varieties in each country
N = 20; % # of goods in each country
eta = 1.05; % substitution parameter across goods close to one
rho = 10;   % substitution parameter across varieties

% Unit cost percentage change = (wage change) X (trade cost change).
% adjust sizes to KxJxN matrices
% country. Column is the importing country.
prod_shr = repmat(prod_shr,1,1,N);
firm_shr = repmat(firm_shr,1,1,N);
C = (P_input(:,1).*delta+P_input(:,2).*(1-delta)).^mygamma; % Take the cost change due to input price change. 
% Convert Jx1 vector to (NKJ)xJ matrix.
C = repmat(C',K,1,N);

%% Now, solve for hat!
% p_before = rand(K*J,N);
p_before0 = ones(K,J,N);
options = optimset('Display','none');

%     gap = 1;
%     while gap>1e-6
%         p_after = p_vrt_after(C_j,S_prev,p_before,rho,eta);
%         gap = norm(p_after(:)-p_before(:));
%         p_before = p_after;
%         disp(gap)
%     end
    p_after = fsolve(@(p_before) p_vrt_after(C,prod_shr,firm_shr,p_before,rho,eta), p_before0,options);
    P_hat = P_agg(prod_shr,p_after,rho,eta,N)'; 
    temp_hat1 = S_k(p_after(:,:,1),prod_shr,rho);
    s_hat = reshape(temp_hat1(:,:,1),K/J,J,J);
    temp_s1 = reshape(prod_shr(:,:,1),K/J,J,J);
    temp_s2 = temp_s1./sum(temp_s1,1);
    pi_hat = squeeze(sum(temp_s2.*s_hat,1));
%% functions for hat variables
function z = P_omega(p_vrt_before, prod_shr, rho)
    z = sum(prod_shr.*p_vrt_before.^(1-rho),1).^(1/(1-rho));
end

function z=S_k(p_vrt_before,prod_shr,rho) % S_hat
    z = (p_vrt_before./P_omega(p_vrt_before,prod_shr,rho)).^(1-rho);
end

function z = epsilon_hat(p_vrt_before,prod_shr, firm_shr,rho,eta)
    temp = 1./((1/eta)*firm_shr.*S_k(p_vrt_before,prod_shr,rho) ...
                +(1/rho)*(1-firm_shr.*S_k(p_vrt_before,prod_shr,rho)));
    z = temp./epsilon_prev(firm_shr,rho,eta);
end

function z = epsilon_prev(firm_shr,rho,eta)
    z = 1./((1/eta)*firm_shr+(1/rho)*(1-firm_shr));
end

function z = p_vrt_after(C_matrix,prod_shr,firm_shr,p_vrt_before,rho,eta)
   z= (epsilon_hat(p_vrt_before,prod_shr, firm_shr,rho,eta).*(epsilon_prev(firm_shr,rho,eta)-1)./ ...
        (epsilon_hat(p_vrt_before,prod_shr, firm_shr,rho,eta).*epsilon_prev(firm_shr,rho,eta)-1)).*C_matrix;
    z = z-p_vrt_before;
    z(isnan(z)) = 0;
end


function z = P_agg(prod_shr,p_vrt_before,rho,eta,N)
    temp = 1/N*P_omega(p_vrt_before, prod_shr,rho).^(1-eta);
    z = sum(real(temp),3).^(1/(1-eta));
end

end