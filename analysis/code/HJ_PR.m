PR_table= readtable('/Users/hj/Documents/GitHub/JapanKoreaTradeConflict/build/input/PR_downstream.xlsx');
% Total sales_dollar is computed from (PR sale of TOK corp X 1/(market
% share of TOK corop)). Sales share is directly sourced from TOK corp's
% sales proportions across countries. Sales_dollar is computed by multiply
% total PR sales to market shares.
JPN_export = PR_table{:,end};
PR_trade = zeros(J_agg,J_agg);
% Country : ROW, Japan, Korea and US
PR_trade(2,:) = [JPN_export(2)+JPN_export(5) JPN_export(3) JPN_export(4) ...
    JPN_export(2)]; 

for i = [1 3 4]
    PR_trade(i,i)=PR_trade(2,i)*0.15;
end

% Use BandTrass data to fill out Korea's PR trade with ROW and US
PR_table =readtable('/Users/hj/Documents/GitHub/JapanKoreaTradeConflict/build/input/PR_trade2014.xlsx');
PR_impval=PR_table{:,6};
PR_trade(1,3) = sum(PR_impval)-PR_impval(10)-PR_impval(5);
PR_trade(4,3) = PR_impval(5);
PR_trade = PR_trade*1e-6;
clear i JPN_export PR_table PR_impval