clc
clear

%% Exercise 1. Match 2020 level of intermediate input trade
% PR import amount in dollar in 2020
calibration_CP_v3
par_fixed.mygamma   = squeeze(sum(Z_aug,1))./reshape(mididx(Y_aug,S_aug,S_aug,J_agg),J_agg,S_aug,S_aug);
par_fixed.mygamma(:,4,3) =  par_fixed.mygamma(:,4,3)+par_fixed.mygamma(:,5,3);% share of intermediate input expenditure in chip production
par_fixed.mygamma(:,5,3) = zeros(J_agg,1);
par_fixed.mygamma(isnan(par_fixed.mygamma )) = 0;
par_fixed.mygamma_sum   = squeeze(sum(reshape(par_fixed.mygamma,J_agg,S_aug,S_aug),2)); % 1-labor share = gamma_sum


% PR_import = [51103157 ; 328294898]; % from world; from japan
% HF_import = [63519010 ; 9375581]; % from world; from japan
% input_import_share_japan = [PR_import(2); HF_import(2)]./[sum(PR_import); sum(HF_import)];
input_import_share_japan = [0.8492; 0.07];
tau             = ones(4,4,5);
% tau(2,3,4)  = 1.00698;
tau(2,3,4)  = 1.262;
tau(2,3,5)  = 1.562;
chip_CP_fun_v3 % share of Japanese input in 2020

%
X_ni = reshape(lastidx(solution.X,S,J,J).*solution.pi(:),J,J,S);
input_share_japan_sol = squeeze(X_ni(2,3,[4 5]))./(solution.X(3,[4 5])'-squeeze(X_ni(3,3,[4 5])));
display([input_import_share_japan input_share_japan_sol])

% share of calibrated model  % share of calibrated model
real_wage_hat = solution.w./prod(solution.P.^par_fixed.alpha,2);
solution_Leontief = solution;
disp('Real wage change in %')
display((solution_Leontief.w./prod(solution_Leontief.P.^par_fixed.alpha,2)-1)*100)
%% Exercise 2. Compare it with uniform tariff increase
calibration_CP_v3

tau             = ones(4,4,5);
% tau(2,3,:)      = 1.000249;
 tau(2,3,:)      = 1.00110;
chip_CP_fun_v3
display(input_import_share_japan); % share of Japanese input in 2020
real_wage_hat_all = solution.w./prod(solution.P.^par_fixed.alpha,2);
display(real_wage_hat_all./real_wage_hat-1)
solution_unif = solution;
display(solution_Leontief.Y(:,3)./solution.Y(:,3)-1)
%% Exercise 3. Compare it with C-D production shape
tau             = ones(4,4,5);
tau(2,3,4)  = 1.262;
tau(2,3,5)  = 1.562;
% par_fixed.mygamma   = squeeze(sum(Z_aug,1))./reshape(mididx(Y_aug,S_aug,S_aug,J_agg),J_agg,S_aug,S_aug);
% par_fixed.mygamma(isnan(par_fixed.mygamma )) = 0;
% par_fixed.mygamma_sum   = squeeze(sum(reshape(par_fixed.mygamma,J_agg,S_aug,S_aug),2)); % 1-labor share = gamma_sum

chip_CD_fun_v1
solution_CD = solution;
% display(input_import_share_japan); % share of Japanese input in 2020
display(solution.pi(2,3,[4 5]));
real_wage_hat_CD = solution_CD.w./prod(solution_CD.P.^par_fixed.alpha,2);

display(real_wage_hat-1)
display(real_wage_hat_CD)
display((real_wage_hat_CD-real_wage_hat)./(real_wage_hat-1))

%% Exercise 4. What would have been the best trade policy for japanse governement?

N=30;
real_wage_hat_japan = zeros(N,1);
i = 1;
while i <= N
    tau             = ones(4,4,5);
    tau(2,3,3)      = 1.001+0.001*(i-1);
    chip_CP_fun_v3
    real_wage_hat_all = solution.w./prod(solution.P.^par_fixed.alpha,2);
    if solution.pi(2,3,3)>0.2 && real_wage_hat_all(2)>1
        real_wage_hat_japan(i) = real_wage_hat_all(2);
        i=i+1;
    end
end

plot((1.001:0.001:1.001+0.001*(N-1)),real_wage_hat_japan)

%% Exercise 5. Conditional on the current increase in tariff barrier, would korean gov's revenge do any help?

%% Exercise 6. What if there japanese components completely removed?



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