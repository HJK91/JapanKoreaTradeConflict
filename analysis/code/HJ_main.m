clc
clear
HJ_WIOD_aug
%% Exercise 1. Match 2020 level of intermediate input trade
% PR import amount in dollar in 2020

PR_import = [51103157 ; 328294898]; % from world; from japan
HF_import = [63519010 ; 9375581]; % from world; from japan
input_import_share_japan = [PR_import(2); HF_import(2)]./[sum(PR_import); sum(HF_import)];
tau             = ones(par_fixed.J,par_fixed.J,par_fixed.S);
tau(2,3,S_aug-1)  = 1.262;
tau(2,3,S_aug)  = 1.562;
HJ_equilibrium % share of Japanese input in 2020

%
X_ni = lastidx(solution.X,S,J,J).*solution.pi;
input_share_japan_sol = squeeze(X_ni(2,3,[S-1 S]))./(solution.X(3,[S-1 S])'-squeeze(X_ni(3,3,[S-1 S])));
display([input_import_share_japan input_share_japan_sol])

% share of calibrated model  % share of calibrated model
real_wage_hat = solution.w./prod(solution.P.^par_fixed.alpha,2);
solution_Leontief = solution;
display(solution_Leontief.w./prod(solution_Leontief.P.^par_fixed.alpha,2)-1)


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