HF_table =readtable('/Users/hj/Documents/GitHub/JapanKoreaTradeConflict/build/input/HF_trade_all.xlsx');
% country code 392 = Japan 410 = Korea, 842 = US, 
% reporter code column: 9 partner code: 12
% trade flow code 1: import 2: export
TradeFlw = HF_table{:,7};
reporter = HF_table{:,9};
partner = HF_table{:,12};
TradeVal = HF_table{:,32};
%country order: ROW, JPN, KOR, US.
jpnkorus = [392 410 842];
reporter_id = unique(reporter);

HF_trade = zeros(J_agg,J_agg);

% Trade of ROW--export to ROW from each country
for i = 1:size(reporter_id)
    idx1 = reporter ==reporter_id(i); % Trade flow from country i...
    idx2 = TradeFlw ==2; % exported to...
    idx3 = sum(partner == jpnkorus,2)==1; 
    idx4 = sum(partner == 0,2)==1; % somewhere other than Japan, Korea and US.
    HF_trade(i+1,1)=sum(TradeVal(idx1&idx2&idx4))-sum(TradeVal(idx1&idx2&idx3)); 
    % export to all countries - export to Japan, Korea and US
end

% Trade of ROW--imports from ROW into each country
for i = 1:size(reporter_id)
    
    idx1 = reporter ==reporter_id(i); % Trade flow into country i...
    idx2 = TradeFlw ==1; % imported from...
    idx3 = sum(partner == jpnkorus,2)==1; 
    idx4 = sum(partner == 0,2)==1; % somewhere other than Japan, Korea and US.
    HF_trade(1,i+1)=sum(TradeVal(idx1&idx2&idx4))-sum(TradeVal(idx1&idx2&idx3)); 
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
            HF_trade(i,j) = sum(TradeVal(idx1&idx2&idx3));
        end
    end
end
HF_hq_korea = [35213805; 43640643; 0 ;361282]*1e-6; % From BandTrass
% deflate chip trade by million dollar
HF_trade = HF_trade*1e-6;
ratio   = HF_hq_korea([1 2 4])./HF_trade([1 2 4],3); % Korea's 10-6digit import ratio
for i = 1:4
    idx1 = 1:4;
    idx1(i) = [];
    HF_trade(idx1,i) = HF_trade(idx1,i).*ratio;
end

%allow for 10% domestic production, compared to imported inputs
HF_trade = HF_trade+ diag(sum(HF_trade)*0.1);
clear HF_hq_korea ratio idx1 idx2 idx3 idx4 i j jpnkorus ...
    reporter partner reporter_id HF_table