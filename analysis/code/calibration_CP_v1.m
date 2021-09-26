clc
clear
load('/Users/hj/Documents/GitHub/JapanKoreaTradeConflict/analysis/input/WIOD.mat')
S = 56;
J = 44;

% Z: Current indexing:  n s i r (r i s n; importing country, importing sector, exporting country, exporting sector)
% Want to make it to s r n i (i n r s) => permute [2 4 1 3]. 
Z = reshape(Z,S,J,S,J);
Z = permute(Z, [2 4 1 3]);
Z = Z(:);

% F: current indexing: n i r (r i n; importing country, exporting country, exporting sector)
% Want to make it r n i (i n r)=> permute [2 3 1].
F = squeeze(sum(reshape(F, S, J, 5 ,J),3));
F = permute(F,[2 3 1]);
F = F(:);

% GO
GO = reshape(GO,S,J)';
GO = GO(:);

% VA
VA = reshape(VA,S,J)';
VA = VA(:);

% Let's make the data super simple.
% Sectors: Others, Chemical, Elec devices, (Chip-H, Chip-L and PR).
% Countries: ROW, Japan, Korea, US.
% Indices-sector: Chemical(11) Electronic devices (17, 18 ,19) 
% Indices-country: JPN(25) KOR (26) US(43) ROW(44)

S_agg = 3;
J_agg = 4;

idx_sector_others                   = (1:56)';
idx_sector_others([11, 17, 18, 19])       = [];
idx_country_others                  = (1:44)';
idx_country_others([25, 26, 43])  = [];
Z_mat   = reshape(Z,J,J,S,S);                            
F_mat   = reshape(F,J,J,S);
Y_mat   = reshape(GO,J,S);
VA_mat  = reshape(VA,J,S);

Z_agg   = zeros(J_agg,J_agg,S_agg,S_agg);
F_agg   = zeros(J_agg,J_agg,S_agg);
VA_agg  = zeros(J_agg,S_agg);
Y_agg   = zeros(J_agg,S_agg);

Z_agg_WIODform = zeros(S_agg,J_agg,S_agg,J_agg);
idx_sector_agg = {idx_sector_others, 11, [17, 18, 19]}'; % Others, Chemical, Electronic devices
idx_country_agg = {idx_country_others, 25, 26, 43}'; % ROW, JPN, KOR, US

for n=1:J_agg
    for s = 1:S_agg
        for i = 1:J_agg
            for r = 1:S_agg
                ss = idx_sector_agg{s}; % importing sector
                rr = idx_sector_agg{r}; % exporting sector
                nn = idx_country_agg{n}; % importing country
                ii = idx_country_agg{i}; % exporting country
                Z_agg(i,n,r,s) = sum(Z_mat(ii,nn,rr,ss),'all');
                %Z_agg_WIODform(r,i,s,n) = sum(Z_mat(ii,nn,rr,ss),'all');
                
                
            end
            F_agg(i,n,s)    = sum(F_mat(ii,nn,ss),'all');
        end
        Y_agg(n,s)    = sum(Y_mat(nn,ss),'all');
        VA_agg(n,s)    = sum(VA_mat(nn,ss),'all');
    end
end

Z_agg_vec = Z_agg(:); 
% To test, Z_agg_WIODform = reshape(Z_agg_WIODform,J_agg*S_agg,J_agg*S_agg);

%% Now, add two sectors: chip and PR.

Z_aug = ones(J_agg,J_agg,S_agg+2,S_agg+2);
% export country, importing country, exporting sector, importing sector
% Country 1: ROW 2: JPN 3: KOR 4: US
% Industry 1: Other 2: Chem 3: Elec 4: Chip 5: PR
% Replace value 2 with data from somewhere else later.
Z_aug(:,:,1:S_agg,1:S_agg) = Z_agg;
Z_aug(:,:,:,S_agg+2) = 0; % PR sector use labor only--no import

Z_aug(:,:,S_agg+2,:)    = 0;    % Prepare to enter PR trade data
Z_aug(2,[3 4],[2 5],4)  = 2;    % PR is exclusively used for chip-H in Korea and US.
Z_aug(:,[3 4],2,4)      = 2;    % Maybe other countries provided chemical inputs for chips
Z_aug(:,[1 2],2,4)      = 0;    % chip sectors of ROW and JPN do not import anything 
Z_aug(:,:,[1 3 4 5],4)  = 0;    % chip-H use PR and chemical only.

Z_aug(:,:,4,:) = 0;         % Prepare to enter chip exporting
Z_aug([3,4],:,4,3) = 2; % US and Korea are the only producers of chips.
% chips are only used to make electronic devices

%% Take the PR and chip global trade data
PR_table= readtable('/Users/hj/Documents/GitHub/JapanKoreaTradeConflict/build/input/PR_downstream.xlsx');
PR_table_importingcountries = PR_table{:,1:2};
temp = PR_table{:,end};
PR_table_importingamount = [temp(3); temp(4); temp(1); temp(2)+temp(5)]*1e-6;

chip_table =readtable('/Users/hj/Documents/GitHub/JapanKoreaTradeConflict/build/input/chip_trade_all.xlsx');
chip_korea_export = chip_table{1:51,end};
chip_us_export = chip_table{52:end,end};

korea_idx_other = (1:51)';
korea_idx_other([15 38 51]) = [];
us_idx_other = (1:length(chip_us_export))';
us_idx_other([70 115 132]) = [];
chip_korea_export = [chip_korea_export([38 51 15]); sum(chip_korea_export(korea_idx_other)) ]*1e-6;
chip_us_export = [chip_us_export([70 115 132]); sum(chip_us_export(us_idx_other))]*1e-6;

char_vec = {('도보로 이동 가능한 거리');('차로 30분 이상 1시간 이내 거리');...
    ('차로 1시간 이상 2시간 이내 거리'); ('차로 2시간 이상 거리')};
for i=1:size(char_vec,1)
    char = char_vec{i,:};
    idx = strcmp(dist,char);
    dist(idx) = {num2str(i-1)};
end

% Sector specific trade
Z_aug(2,:,6,4) = PR_table_importingamount;
Z_aug(1,3,6,4) = 1672113842*1e-6;

% Assume Micron inc. combine labor and japanese inputs

Z_aug(2,:,4,3) = chip_korea_export;
Z_aug(3,:,4,3) = chip_us_export;
 
clear chip_korea_export chip_us_export PR_table chip_table ...
    PR_table_importingamount temp
% Augment F, Y and VA from Z_aug.
F_aug   = zeros(J,J,S_agg+2);
Y_aug   = zeros(J,S_agg+2);
VA_aug  = zeros(J,S_agg+2);
F_aug(:,:,1:S_agg)  = F_agg;
Y_aug(:,1:S_agg)    = Y_agg;
VA_aug(:,1:S_agg)   = V_agg;
temp_Y = squeeze(sum(Z_aug,[2 4]));
temp_V = squeeze(sum(Z_aug,[1,3]));
Y_aug(:,S_agg+1:end) = temp_Y(:,S_agg+1:end);
VA_aug(:,S_agg+1:end) = Y_aug(:,S_agg+1:end) - temp_V(:,S_agg+1:end);

% feed in the data into the function to be minimized
%% Take the labor and deficit data
M = squeeze(sum(Z,3));
L = randn(J,1); % Change this to actual data
D_0     = squeeze(sum(M,1))-squeeze(sum(M,2)); 
WIOD_simple = [Z_aug(:); F_aug(:); Y_aug(:); VA_aug(:)];
data_simple = [L_0; D_0];
%% Set parameters for calibration
par.T       = rand(J*S,1)*0.2+1;        % Technology parameter T_j^s
par.tau     = rand(S*J^2,1)+1;          % Trade cost parameters tau_ni^s
par.theta   = 5+rand(S,1)*0.1;          % CA parameter theta
par.sigma   = par.theta+1-rand(S,1);    % Substitution parameter sigma^s
par.gamma   = squeeze(sum(Z_aug,1))./reshape(mididx(Y_aug,S,S,J),J,S,S);                 
par.gamma_sum       = 1-par.gamma; % 1-labor share = gamma_sum
par.upsilon = prod(reshape(par.gamma.^(-par.gamma),J,S,S),2); 
par.upsilon = par.upsilon(:).*(1-par.gamma_sum).^(-(1-par.gamma_sum));  % As in Caliendo and Parro (2015)
par.alpha   = squeeze(sum(F_aug,1))./squeeze(sum(F_aug,[1,3]));
par.kappa   = gamma(1+(1-par.sigma)./par.theta).^(1./(1-par.sigma));    % As in Caliendo and Parro (2015)
par.K       = 2; % the ratio between chemical and chips
par.S       = size(Z_aug,3)-2;
par.J       = size(Z_aug,1);





% Put the calibration work

function z = mididx(x,idx1,idx2,idx3)
    % expand the middle index
    x_temp = repmat(reshape(x,idx3,idx1),idx2,1);
    z = x_temp(:);
end