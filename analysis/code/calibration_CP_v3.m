
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

S_agg = 2;
J_agg = 4;

idx_sector_others                   = (1:56)';
idx_sector_others([17, 18, 19])       = [];
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

% Z_agg_WIODform = zeros(S_agg,J_agg,S_agg,J_agg);
idx_sector_agg = {idx_sector_others, [17, 18, 19]}'; % Others, Electronic devices
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

% Now, add three sectors: chip, PR, and HF.

Z_aug = zeros(J_agg,J_agg,S_agg+3,S_agg+3);
% export country, importing country, exporting sector, importing sector
% Country 1: ROW 2: JPN 3: KOR 4: US
% Industry 1: Other 2: Elec 3: Chip 4: PR 5: HF
% Replace value 2 with data from somewhere else later.
Z_aug(:,:,1:S_agg,1:S_agg) = Z_agg;
Z_aug(:,:,:,S_agg+2) = 0; % PR sector use labor only--no import

% chips are only used to make electronic devices

% Take the PR and chip global trade data
PR_table= readtable('/Users/hj/Documents/GitHub/JapanKoreaTradeConflict/build/input/PR_downstream.xlsx');
% PR_table_importingcountries = PR_table{:,1:2};
temp = PR_table{:,end};
PR_table_importingamount = [temp(2)+temp(5);temp(3); temp(4); temp(1) ]*1e-6;

chip_table =readtable('/Users/hj/Documents/GitHub/JapanKoreaTradeConflict/build/input/chip_trade_all.xlsx');
HF_table =readtable('/Users/hj/Documents/GitHub/JapanKoreaTradeConflict/build/input/HF_trade_all.xlsx');
chip_reporter   = chip_table{:,11};
HF_reporter     = HF_table{:,11};
chip_partner    = chip_table{:,14};
HF_partner      = HF_table{:,14};
chip_trade      = chip_table{:,8};
HF_trade        = HF_table{:,8};
chip_trade_val  = chip_table{:,32};
HF_trade_val    = HF_table{:,32};
% chip_korea_export = chip_table{1:51,end};
% chip_us_export = chip_table{52:end,end};
% 
% korea_idx_other = (1:51)';
% korea_idx_other([15 38 51]) = [];
% us_idx_other = (1:length(chip_us_export))';
% us_idx_other([70 115 132]) = [];
% chip_korea_export = [chip_korea_export([38 51 15]); sum(chip_korea_export(korea_idx_other)) ]*1e-6;
% chip_us_export = [chip_us_export([70 115 132]); sum(chip_us_export(us_idx_other))]*1e-6;

char_vec = {('WLD');('JPN');('KOR');('USA') };
char_trade= {'Export', 'Import'};
for i=1:size(char_vec,1)  % exporter
    for n = 1:size(char_vec,1) % importer
        if i == n 
            continue
        end
            char1   = char_vec{i,1}; % exporter name
            char2   = char_vec{n,1}; % importer country name
            char3   = char_trade{1,1};
            
            idx1_chip = strcmp(chip_reporter,char1);
            idx2_chip = strcmp(chip_partner,char2);
            idx3_chip = strcmp(chip_trade,char3);
            
            idx1_HF = strcmp(HF_reporter,char1);
            idx2_HF = strcmp(HF_partner,char2);
            idx3_HF = strcmp(HF_trade,char3);            
        if i ==1
            % if n = 2, world->Japan semiconductor trade
            idx1_chip   = strcmp(chip_partner,char1); % partner: world
            idx2_chip   = strcmp(chip_reporter,char2); % reporter: japan
            idx3_chip   = strcmp(chip_trade,char_trade{1,2}); % Japan's import
            idx_chip    = idx1_chip.*idx2_chip.*idx3_chip>0; % Japan's imports from all countries
            
            idx1_HF     = strcmp(HF_partner,char1); % partner: world
            idx2_HF     = strcmp(HF_reporter,char2); % reporter: japan
            idx3_HF     = strcmp(HF_trade,char_trade{1,2}); % Japan's import
            idx_HF      = idx1_HF.*idx2_HF.*idx3_HF>0; % Japan's imports from all countries
            % now, imports from US and Korea -- substracted from imoprts
            % from WLD => imports from ROW.
            char_vec_temp       = char_vec;
            char_vec_temp([1, n])= [];
            idx_others_chip     = zeros(length(idx_chip),length(char_vec_temp));
            idx_others_HF       = zeros(length(idx_HF),length(char_vec_temp));
            for j = 1: length(char_vec_temp)
                idx_others_chip(:,j) = strcmp(chip_partner,char_vec_temp{j,1});
                idx_others_HF(:,j) = strcmp(HF_partner,char_vec_temp{j,1});
            end
            idx_others_chip     = sum(idx_others_chip,2);
            idx_others_chip     = idx_others_chip.*idx2_chip.*idx3_chip>0;
            
            idx_others_HF       = sum(idx_others_HF,2);
            idx_others_HF       = idx_others_HF.*idx2_HF.*idx3_HF>0;
            
            Z_aug(i,n,3,2)      = (chip_trade_val(idx_chip)-sum(chip_trade_val(idx_others_chip)))*1e-6;
            Z_aug(i,n,5,3)      = (HF_trade_val(idx_HF)-sum(HF_trade_val(idx_others_HF)))*1e-6;
            
        else
            idx_chip = idx1_chip.*idx2_chip.*idx3_chip>0;
            idx_HF = idx1_HF.*idx2_HF.*idx3_HF>0;
            
            Z_aug(i,n,3,2) = chip_trade_val(idx_chip)*1e-6;
            Z_aug(i,n,5,3) = HF_trade_val(idx_HF)*1e-6;
            
        end
    end
end
%%

chip_table =readtable('/Users/hj/Documents/GitHub/JapanKoreaTradeConflict/build/input/dram_exports_all.xlsx');
Z_aug(3,3,3,2) = chip_table{51,end}*1e-6;
Z_aug(4,4,3,2) = chip_table{183,end}*1e-6;
Z_aug(2,2,3,2) = (2706.*1e+9/105)*1e-6*0.24; % Toshiba's domestic memory sales
Z_aug(2,:,4,3) = PR_table_importingamount;
Z_aug([1 4],3,4,3) = [700133 13352807]*1e-6; % use trade data from Korean custom 
Z_aug(2,[1 2 4],4,3) = Z_aug(2,[1 2 4],4,3)*0.3; % product memory chip industry share
HF_hq_korea = [35213805; 43640643; 0 ;361282]*1e-6;
for i = 1: 4
    pure_HF_factor = HF_hq_korea./Z_aug(:,3,5,3);
    pure_HF_factor(isnan(pure_HF_factor)) = 1;
    chip_factor =  squeeze(sum(Z_aug,[2 4]));
    chip_factor =  chip_factor(i,3)./chip_factor(3,3);
    Z_aug(:,i,5,3) = Z_aug(:,i,5,3).*(pure_HF_factor)*chip_factor;
    if i ==2
        continue
    else
        Z_aug(i,i,4,3) = Z_aug(2,i,4,3)*0.15; % accommodate the multinational production of PR
    end
end
Z_aug(isnan(Z_aug))=0;
 % 12,224 million yen: 2016 Stella cheifa's HF sales. 
%%
S_aug = S_agg+3;
clear chip_korea_export chip_us_export PR_table chip_table ...
     temp idx1_chip idx2_chip idx3_chip idx1_HF idx2_HF idx3_HF ...
     idx_chip idx_HP idx_others_chip idx_others_HF chip_partner HF_partner ...
     char_vec char_trade n i ii j s ss r rr idx_country_agg idx_country_others ...
     char1 char2 char3 char_vec_temp chip_reporter chip_trade chip_trade_val ...
     temp_V temp_Y HF_reporter HF_trade HF_trade_val idx_sector_agg ...
     idx_sector_others nn PR_table_importingamount PR_table HF_table ...
     idx_HF 

 
% Augment F, Y and VA from Z_aug.
F_aug                   = zeros(J_agg,J_agg,S_aug);
F_aug(:,:,1:S_agg)      = F_agg;
temp_Y                  = squeeze(sum(Z_aug,[2 4]));
temp_V                  = squeeze(sum(Z_aug,[1,3]));
Y_aug = temp_Y+squeeze(sum(F_aug,2));
VA_aug   = Y_aug - temp_V;
X_ind = squeeze(sum(Z_aug,4))+F_aug;
X_ttl= squeeze(sum(Z_aug,[1,4]))+squeeze(sum(F_aug,1));

% feed in the data into the function to be minimized
%% Take the labor and deficit data
M           = squeeze(sum(Z_aug,3));
L_0         = [5e+9; 1.2e+7; 5e+7; 3e+8].*1e-6; % Change this to actual data
D_0         = squeeze(sum(M,1))-squeeze(sum(M,2)); 
data.Y      = Y_aug;
data.X_ind  = X_ind;
data.X_ttl  = X_ttl;
data.VA     = VA_aug;
data.F      = F_aug;
data.Z      = Z_aug;
data.L      = L_0;
data.D      = D_0;

%% Set parameters for counterfactuals: parameters are array forms


par_fixed.theta     = [4.55; 10.6; 10.6; 4.75; 4.75];    % CA parameter theta: change this to estimated one later

par_fixed.mygamma   = squeeze(sum(Z_aug,1))./reshape(mididx(Y_aug,S_aug,S_aug,J_agg),J_agg,S_aug,S_aug);
% par_fixed.mygamma(:,4,3) = sum(Z_aug(:,:,4,3)+Z_aug(:,:,5,3),1)'./Y_aug(:,3); % share of intermediate input expenditure in chip production
par_fixed.mygamma(isnan(par_fixed.mygamma )) = 0;
par_fixed.mygamma_sum   = squeeze(sum(reshape(par_fixed.mygamma,J_agg,S_aug,S_aug),2)); % 1-labor share = gamma_sum

par_fixed.alpha     = squeeze(sum(F_aug,1))./squeeze(sum(F_aug,[1,3]))';
par_fixed.alpha(isnan(par_fixed.alpha)) = 0;

par_fixed.delta     = (sum(Z_aug(:,:,4,3),1)./(sum(Z_aug(:,:,4,3),1)+sum(Z_aug(:,:,5,3),1)))';
% par_fixed.mypi      = squeeze(sum(Z_aug,3))./reshape(lastidx(squeeze(sum(Z_aug,[1,3])),S_aug,J_agg,J_agg),J_agg,J_agg,S_aug);
% par_fixed.mypi(:,:,4) = eye(J_agg,J_agg);
% par_fixed.mypi(:,:,5) = eye(J_agg,J_agg);

par_fixed.mypi = X_ind./reshape(lastidx(X_ttl,S_aug,J_agg,J_agg),J_agg,J_agg,S_aug);
mypi = X_ind(:)./lastidx(X_ttl,S_aug,J_agg,J_agg);
par_fixed.mypi=reshape(mypi,J_agg,J_agg,S_aug);
data.D = sum(X_ttl,2)-sum(reshape(lastidx(data.X_ttl,S_aug,J_agg,J_agg).*par_fixed.mypi(:),J_agg,J_agg,S_aug),[2 3]);
par_fixed.wL        = sum(data.VA,2);
par_fixed.S         = size(Z_aug,3);
par_fixed.J         = size(Z_aug,1);


% Put the calibration work
%% check the validity of parameters
Y_accounting = squeeze(sum(reshape(par_fixed.mypi(:).*lastidx(data.X_ttl(:),S_aug,J_agg,J_agg),J_agg,J_agg,S_aug),2));
Y_gap = (data.Y(:)-Y_accounting(:))./data.Y(:);
X_gap = data.X_ttl-squeeze(sum(reshape(par_fixed.mygamma(:).*mididx(data.Y,S_aug,S_aug,J_agg),J_agg,S_aug,S_aug),3)) ...
    -par_fixed.alpha.*(par_fixed.wL+sum(data.D,2));
TB_gap1 = sum(data.X_ttl,2)...
    -sum(reshape(lastidx(data.X_ttl,S_aug,J_agg,J_agg).*par_fixed.mypi(:),J_agg,J_agg,S_aug),[1 3])';
TB_gap2 = sum(reshape(lastidx(data.X_ttl,S_aug,J_agg,J_agg).*par_fixed.mypi(:),J_agg,J_agg,S_aug),[2 3])+sum(data.D,2);
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