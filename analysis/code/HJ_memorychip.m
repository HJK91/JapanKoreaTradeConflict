chip_trade =readtable('/Users/hj/Documents/GitHub/JapanKoreaTradeConflict/build/input/memorychip_trade.csv');
% country code 392 = Japan 410 = Korea, 842 = US, 
% reporter code column: 9 partner code: 12
% trade flow code 1: import 2: export
TradeFlw = chip_trade{:,7};
reporter = chip_trade{:,9};
partner = chip_trade{:,12};
TradeVal = chip_trade{:,32};
%country order: ROW, JPN, KOR, US.
jpnkorus = [392 410 842];
reporter_id = unique(reporter);

chip_trade = zeros(J_agg,J_agg);

% Trade of ROW--export to ROW from each country
for i = 1:size(reporter_id)
    idx1 = reporter ==reporter_id(i); % Trade flow from country i...
    idx2 = TradeFlw ==2; % exported to...
    idx3 = sum(partner == jpnkorus,2)==1; 
    idx4 = sum(partner == 0,2)==1; % somewhere other than Japan, Korea and US.
    chip_trade(i+1,1)=sum(TradeVal(idx1&idx2&idx4))-sum(TradeVal(idx1&idx2&idx3)); 
    % export to all countries - export to Japan, Korea and US
end

% Trade of ROW--imports from ROW into each country
for i = 1:size(reporter_id)
    
    idx1 = reporter ==reporter_id(i); % Trade flow into country i...
    idx2 = TradeFlw ==1; % imported from...
    idx3 = sum(partner == jpnkorus,2)==1; 
    idx4 = sum(partner == 0,2)==1; % somewhere other than Japan, Korea and US.
    chip_trade(1,i+1)=sum(TradeVal(idx1&idx2&idx4))-sum(TradeVal(idx1&idx2&idx3)); 
    % export to all countries - export to Japan, Korea and U
end

for i = 2:size(reporter_id)+1
    for j = 2:size(reporter_id)+1
        if i == j
            continue
        else
            idx1 = reporter==reporter_id(i-1); % Country i exports to
            idx2 = partner==reporter_id(j-1); % Country j
            idx3 = TradeFlw==2;
            chip_trade(i,j) = sum(TradeVal(idx1&idx2&idx3));
        end
    end
end

% deflate chip trade by million dollar
chip_trade = chip_trade*1e-6;

% Domestic use of memory chip: first, take the global memory chip
% production from US congress report
chip_prod = 106*1e+9*1e-6; % 106 billion dollars = 106000 million dollars
% Get the regional production share from fabriacation capacity shares
fab_share = [43-24 18 26 13]'*(100/(100-24)); % reflect the Taiwan's capacity

fab_share([1 3 4]) = fab_share([1 3 4])*(100/(100-(fab_share(2) -9)));
fab_share(2) = 9;% fix Japan's capcity to 9%.
fab_share = fab_share*.01;
ttl_prod = chip_prod*fab_share;
chip_trade = chip_trade+diag(ttl_prod-sum(chip_trade,2));

clear ttl_prod fab_share chip_prod idx1 idx2 idx3 idx4 i j jpnkorus ...
    reporter partner TradeVal TradeFlw