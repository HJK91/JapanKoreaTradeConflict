clc
clear
HJ_WIOD_aug
%% Exercise 1. Match 2020 level of intermediate input trade
% PR import amount in dollar in 2020

% PR_import = [51103157 ; 328294898]; % from world; from japan
% HF_import = [63519010 ; 9375581]; % from world; from japan
% input_import_share_japan = [PR_import(2); HF_import(2)]./[sum(PR_import); sum(HF_import)];
input_import_share_japan = [0.8492; 0.07]; % Japan's PR and HF share in late 2019
tau             = ones(par_fixed.J,par_fixed.J,par_fixed.S);
tau(2,3,S_aug-1)  = 1.3375;
tau(2,3,S_aug)  = 2.0758;
% tau             = ones(par_fixed.J,par_fixed.J,par_fixed.S);
% tau(2,3,S_aug-1)  = 1.262;
% tau(2,3,S_aug)  = 1.562;
% Take the HF and PR market share of firms
% Firms are ONLY Japanese firms and their local factories 
% From 2014 TOK corp's financial statements, take the PR market share
% across firms

% No market share data for HF: use uniform distribution instead.

HJ_equilibrium % share of Japanese input in 2020

%
X_ni = lastidx(solution.X,S,J,J).*solution.pi;
input_share_japan_sol = squeeze(X_ni(2,3,[S-1 S]))./(solution.X(3,[S-1 S])'-squeeze(X_ni(3,3,[S-1 S])));
display([input_import_share_japan input_share_japan_sol])

% share of calibrated model  % share of calibrated model
real_wage_hat = solution.w./prod(solution.P.^par_fixed.alpha,2);
solution_HJ = solution;
disp((solution_HJ.w./prod(solution_HJ.P.^par_fixed.alpha,2)-1)*100)


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