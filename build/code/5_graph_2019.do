cd "/Users/hj/Documents/GitHub/JapanKoreaTradeConflict/build/input"


import excel "PR_trade_2019.xlsx", firstrow clear
reshape wide Export_dollar Export_kg Import_dollar Import_kg, i(Month) j(Country) string
egen Import_dollarAll = rowtotal(Import_dollar*)
gen Import_dollarFrac=Import_dollarJapan/Import_dollarAll
label variable Import_dollarJapan "PR import from Japan"
label variable Import_dollarAll "PR import from the world"
label variable Import_dollarFrac "Japan's import share"
twoway (bar Import_dollarAll Month, fcolor(%80) lcolor(none%0)) ///
		(bar Import_dollarJapan Month, lcolor(none%0) lwidth(none%0)) ///
		(line Import_dollarFrac Month, yaxis(2)), ///
 		ytitle("Import amount(Dollar)") yscale(range(0 1)) ylabel(1e+7(2e+7)6e+7) ///
 		ylabel(0.6(0.1)1, axis(2)) ytitle(Import share of Japan, axis(2)) xtitle("Month, 2019") ///
 		xlabel(1(2)12) xline(7, lwidth(thick) lpattern(dash) lcolor(black) noextend) title("Monthly PR Import of Korea")


cd "/Users/hj/Documents/GitHub/JapanKoreaTradeConflict/analysis/output"
graph export PR_Japan_Share_2019.png, replace

cd "/Users/hj/Documents/GitHub/JapanKoreaTradeConflict/build/input"

import excel "HF_trade_2019.xlsx", firstrow clear
drop if missing(Country)
reshape wide Export_dollar Export_kg Import_dollar Import_kg, i(Month) j(Country) string
egen Import_dollarAll = rowtotal(Import_dollar*)
gen Import_dollarFrac=Import_dollarJapan/Import_dollarAll
label variable Import_dollarJapan "HF import from Japan"
label variable Import_dollarAll "HF import from the world"
label variable Import_dollarFrac "Japan's import share"
twoway (bar Import_dollarAll Month, fcolor(%80) lcolor(none%0)) ///
		(bar Import_dollarJapan Month, lcolor(none%0) lwidth(none%0)) ///
		(line Import_dollarFrac Month, yaxis(2)), ///
 		ytitle("Import amount(Dollar)") yscale(range(0 1)) ylabel(1e+7) ///
 		ylabel(0(0.1)0.5, axis(2)) ytitle(Import share of Japan, axis(2)) xtitle("Month, 2019") ///
 		xlabel(1(2)12) xline(7, lwidth(thick) lpattern(dash) lcolor(black) noextend) title("Monthly HF Import of Korea")


cd "/Users/hj/Documents/GitHub/JapanKoreaTradeConflict/analysis/output"
graph export HF_Japan_Share_2019.png, replace
