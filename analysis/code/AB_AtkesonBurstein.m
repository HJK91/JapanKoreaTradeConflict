function z = AB_AtkesonBurstein(w,tau)

% parameters
J = 4;
K = 10; % # of varieties in each country
N = 20; % # of goods in each country
eta = 1.05; % substitution parameter across goods close to one
rho = 10;   % substitution parameter across varieties
rng(1)
D = 1.1+rand(J,J)/5; % matrix of trade costs
D = D-diag(diag(D))+eye(J,J);
%%
D(2,3) = 1.4;
%%
% Frechet parameter
theta = 4.75; % trade elasticity of chemical sector
Omega = 1/eta-1/rho;
% draw firm productivities

w_before           = ones(J,1);
w_after         = w_before.*w;
q_underbar  = .1;
T           = K*q_underbar^theta;
myPhi1       = T.*w_before.^(-theta);
myPhi2       = T.*w_after.^(-theta);
invC_matrix = rand(K,N,J);
C_matrix1    = zeros(K,N,J);
C_matrix2    = zeros(K,N,J);
flag        = zeros(K,N,J);
C0          = 1;
options     = optimoptions('fsolve','Display','none');
for j = 1:4
k=1;
while k <= K
    n = 1;
    while n <=N
        q = invC_matrix(k,n,j);
        
        [C_matrix1(k,n,j),~,flag(k,n,j)] = fsolve(@(C) G(C,q,myPhi1(j),theta,K), C0, options);
        [C_matrix2(k,n,j),~,flag(k,n,j)] = fsolve(@(C) G(C,q,myPhi2(j),theta,K), C0, options);
        if flag(k,n,j)==1
            n=n+1;
        else
            C0          = exp(randn(1,1));
        end
    end
    k=k+1;
end

end

%%
X = 1;
Y0 = randn(K*N*J,1);
Y_sol1 = zeros(K*N*J,J);
Y_sol2 = zeros(K*N*J,J);
P_sol1 = zeros(1,J);
P_sol2 = zeros(1,J);
for j = 1:J
    C_mat_D1=zeros(K*J,N);
    C_mat_D2=zeros(K*J,N);
    for jj = 1:J
        C_mat_D1((jj-1)*K+1:jj*K,:) = C_matrix1(:,:,jj)*D(jj,j);
        if jj == 2 && j == 3
            C_mat_D2((jj-1)*K+1:jj*K,:) = C_matrix2(:,:,jj)*D(jj,j)*tau; % increase trade costs by tau (hat).
        else
            C_mat_D2((jj-1)*K+1:jj*K,:) = C_mat_D1((jj-1)*K+1:jj*K,:);
        end
    end
    Y_sol1(:,j) = fsolve(@(Y) FOC(C_mat_D1,Y,X,rho,eta,K,N,J),Y0,options);
    Y_sol2(:,j) = fsolve(@(Y) FOC(C_mat_D2,Y,X,rho,eta,K,N,J),Y0,options);
    Y_sol1(:,j) = exp(Y_sol1(:,j));
    Y_sol2(:,j) = exp(Y_sol2(:,j));
    P_sol1(:,j) = P_agg(C_mat_D1,reshape(Y_sol1(:,j),K*J,N),rho,eta);
    P_sol2(:,j) = P_agg(C_mat_D2,reshape(Y_sol2(:,j),K*J,N),rho,eta);
end
%
z = P_sol2./P_sol1;


function z = G(C,q,Phi,theta,K)
    temp = (Phi*C^theta).^(0:K-1)./factorial((0:K-1));
    z = q-(1-exp(-Phi*C^theta*sum(temp)));
end

function z = Y_omega(Y_vrt, rho)
    z = sum(Y_vrt.^((rho-1)/rho)).^(rho/(rho-1));
end

function z=S_k(Y_vrt,rho)
    z = (Y_vrt./Y_omega(Y_vrt,rho)).^(rho/(rho-1));
end

function z = endo_epsilon(Y_vrt,rho,eta)
    z = 1./((1/eta)*S_k(Y_vrt,rho)+(1/rho)*(1-S_k(Y_vrt,rho)));
end

function z = price_vrt(C_matrix,Y_vrt,rho,eta)
   z= (endo_epsilon(Y_vrt,rho,eta)./(endo_epsilon(Y_vrt,rho,eta)-1)).*C_matrix;
end

function z = P_omega(C_matrix,Y_vrt,rho,eta)
    z = sum(price_vrt(C_matrix,Y_vrt,rho,eta).^((rho-1)/rho)).^(rho/(rho-1));
end

function z = P_agg(C_matrix,Y_vrt,rho,eta)
    z = sum(P_omega(C_matrix,Y_vrt,rho,eta).^((eta-1)/eta),2).^(eta/(eta-1));
end

function z = FOC(C_matrix,Y,X,rho,eta,K,N,J)
    Omega = 1/eta-1/rho;
    Y_vrt = reshape(exp(Y),K*J,N);
    z = (Y_vrt.^(-1/rho)).*(Y_omega(Y_vrt, rho).^(-Omega))...
        .*(P_agg(C_matrix,Y_vrt,rho,eta).^((eta-1)/eta)) ...
        .*(X^(1/eta))-C_matrix;
end
end