clc
clear
% Raw file
% From row 7-62 AUS -2514 US - ROW 2470 2471 Total intermediate
% consumption 2476 VA 2478 GO
% From col E-BH AUS -CNU USA -CPX ROW -CPY Final consumption, AUS, last
% column: total output.
table= readtable('/Users/hj/Documents/GitHub/JapanKoreaTradeConflict/build/input/WIOT2014_Nov16_ROW.xlsx');
%%
% S =56, J = 44.
% Abosorbing sector starts from column 5. 
% 1~4: row items, 5~2468: intermediate input trade 2469~2688: final
% consumption, 2689: GO
% Final consumption: 5 entries per country.

% Intermediate input trade table
Z = table{1:2464,5:2468};
F = table{1:2464,2469:2688};
GO = table{1:2464,end};
VA = table{2470,5:2468}';


save('/Users/hj/Documents/GitHub/JapanKoreaTradeConflict/analysis/input/WIOD.mat', 'Z', 'F', 'GO', 'VA')