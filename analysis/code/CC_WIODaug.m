clc
clear
% Call WIOD 
load('/Users/hj/Documents/GitHub/JapanKoreaTradeConflict/analysis/input/WIOD.mat')
S   = 56;     % # of manufacturing industries (HS 2-digit level)
J   = 44;     % # of countries

% Z: Current indexing:  n s i r (matrix index: r i s n; importing country, importing sector, exporting country, exporting sector)
% Want to make it to s r n i (matrix index: i n r s) => permute [2 4 1 3]. 
Z   = reshape(Z,S,J,S,J);
Z   = permute(Z, [2 4 1 3]);
F   = squeeze(sum(reshape(F, S, J, 5 ,J),3)); % sum over fiver different types of consumption
F   = permute(F,[2 3 1]);
Y  = reshape(GO,S,J)';
VA  = reshape(VA,S,J)';

S = 22;
Z = Z(:,:,1:S,1:S);
F = F(:,:,1:S);
Y = Y(:,1:S);
VA = VA(:,1:S);

% Let's aggregate the data over countries.
% Sectors: Others, Elec devices, (Chip, PR and HF).
% Countries: ROW, Japan, Korea, US.
% Indices-sector: Electronic devices (17, 18 ,19) 
% Indices-country: JPN(25) KOR (26) US(43) ROW(44)

S_agg = 20; % aggregate three electronic divice sector
J_agg = 4;
S_aug = S_agg+3; % including subindustries


idx_sector_others                   = (1:S)';
idx_sector_others([17, 18, 19])       = [];
idx_country_others                  = (1:J)';
idx_country_others([25, 26, 43])  = [];

Z_agg   = zeros(J_agg,J_agg,S_agg,S_agg);
F_agg   = zeros(J_agg,J_agg,S_agg);
VA_agg  = zeros(J_agg,S_agg);
Y_agg   = zeros(J_agg,S_agg);

% Z_agg_WIODform = zeros(S_agg,J_agg,S_agg,J_agg);
idx_sector_agg = [num2cell(idx_sector_others); [17, 18, 19]']; % Others, Electronic devices
idx_country_agg = {idx_country_others, 25, 26, 43}'; % ROW, JPN, KOR, US

for n=1:J_agg
    for s = 1:S_agg
        for i = 1:J_agg
            for r = 1:S_agg
                ss = idx_sector_agg{s}; % importing sector
                rr = idx_sector_agg{r}; % exporting sector
                nn = idx_country_agg{n}; % importing country
                ii = idx_country_agg{i}; % exporting country
                Z_agg(i,n,r,s) = sum(Z(ii,nn,rr,ss),'all');
            end
            F_agg(i,n,s)    = sum(F(ii,nn,ss),'all');
        end
        Y_agg(n,s)    = sum(Y(nn,ss),'all');
        VA_agg(n,s)    = sum(VA(nn,ss),'all');
    end
end

clear n s r i ss rr nn ii temp idx_sector_agg idx_country_agg ...
    idx_sector_others idx_country_others
% Now, add three sectors: chip, PR, and HF. 
% export country, importing country, exporting sector, importing sector
% Country 1: ROW 2: JPN 3: KOR 4: US
% Industry 1~S_agg-1(53): Others 54: Elec 55: Chip 56: PR 57: HF
% Replace value 2 with data from somewhere else later.
Z_aug = zeros(J_agg,J_agg,S_agg+3,S_agg+3);
Z_aug(:,:,1:S_agg,1:S_agg) = Z_agg;

HJ_memorychip
clc
Z_aug(:,:,S_agg+1,S_agg) = chip_trade; % S_agg+1 exporting sector-memory chip/ S_agg: importing sector-electronics
HJ_PR
Z_aug(:,:,S_agg+2,S_agg+1) = PR_trade;
HJ_HF
Z_aug(:,:,S_agg+3,S_agg+1) = HF_trade;

% Augment F, Y and VA from Z_aug.
F_aug                   = zeros(J_agg,J_agg,S_aug);
F_aug(:,:,1:S_agg)      = F_agg;
temp_Y                  = squeeze(sum(Z_aug,[2 4]));
temp_V                  = squeeze(sum(Z_aug,[1,3]));
Y_aug = temp_Y+squeeze(sum(F_aug,2));
VA_aug   = Y_aug - temp_V;
X_ind = squeeze(sum(Z_aug,4))+F_aug; % X_ni
X_ttl= squeeze(sum(Z_aug,[1,4]))+squeeze(sum(F_aug,1)); % X_n

clear temp_Y temp_V
% feed in the data into the function to be minimized
%% Take the labor and deficit data
data.Y      = Y_aug;
data.X_ind  = X_ind;
data.X_ttl  = X_ttl;
data.VA     = VA_aug;
data.F      = F_aug;
data.Z      = Z_aug;

%% Set parameters for counterfactuals: parameters are array forms
par_fixed.theta     = [4.55*ones(S_agg-1,1); 10.6; 10.6; 4.75; 4.75];    % CA parameter theta: change this to estimated one later
par_fixed.mygamma   = squeeze(sum(Z_aug,1))./mididx(Y_aug,S_aug,S_aug,J_agg);
par_fixed.mygamma(isnan(par_fixed.mygamma )) = 0;
par_fixed.mygamma_sum   = squeeze(sum(reshape(par_fixed.mygamma,J_agg,S_aug,S_aug),2)); % 1-labor share = gamma_sum
par_fixed.alpha     = squeeze(sum(F_aug,1))./squeeze(sum(F_aug,[1,3]))';
par_fixed.alpha(isnan(par_fixed.alpha)) = 0;
par_fixed.delta     = (sum(Z_aug(:,:,S_aug-1,S_aug-2),1)./(sum(Z_aug(:,:,S_aug-1,S_aug-2),1)+sum(Z_aug(:,:,S_aug,S_aug-2),1)))';
par_fixed.mypi = X_ind./lastidx(X_ttl,S_aug,J_agg,J_agg);
par_fixed.mypi=reshape(par_fixed.mypi,J_agg,J_agg,S_aug);
par_fixed.wL        = sum(data.VA,2);
par_fixed.S         = size(Z_aug,3);
par_fixed.J         = size(Z_aug,1);
data.D = sum(X_ttl,2)-sum(lastidx(data.X_ttl,S_aug,J_agg,J_agg).*par_fixed.mypi,[2 3]);

% Compute market shares: Take the HF and PR market share of firms
% Firms are ONLY Japanese firms and their local factories 
% From 2014 TOK corp's financial statements, take the PR market share
% across firms
S_prev =[25 ;18.5; 13.1; 13; 9.7; 6.5; 14.2]/100;
S_prev = [S_prev;S_prev.*0.15];
S_prev = S_prev/sum(S_prev);
data.S_prev(:,1) = S_prev;
S_prev =ones(length(S_prev(:,1))/2,1)/length(S_prev(:,1));
S_prev = [S_prev;S_prev.*0.1];
S_prev = S_prev/sum(S_prev);
% No market share data for HF: use uniform distribution instead.
data.S_prev(:,2) = S_prev ;

% functions that control index matching

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