
// plots for imported goods
cd ..
cd output

foreach x in "HF" "PR" "PI" {
	twoway (bar KW_`x'_ImportQuantity_Dollar Year) (bar KJ_`x'_ImportQuantity_Dollar Year), ytitle(`x' Import(Dollar)) ylabel(#3)
	graph export historic_`x'_trade.png, replace
	twoway (bar Japan_`x'_share Year) , ytitle(Japan's import share of `x') ylabel(#3)
	graph export historic_`x'_Japan_Share.png, replace
}

// Bar plots for exported goods
foreach x in "Dram" "Flash" {
	twoway (bar `x'_ExportQuantity_Dollar Year), ytitle(`x' Export(Dollar)) ylabel(#3)
	graph export historic_`x'_export_dollars.png,replace
}

foreach x in "Dram" "Flash" {
	twoway (bar `x'_ExportQuantity_Kg Year), ytitle(`x' Export(Kg)) ylabel(#3)
	graph export historic_`x'_trade_kg.png,replace
}


//HF graph

gen HF_Japan_share = KJ_HF_ImportQuantity_Dollar/KW_HF_ImportQuantity_Dollar
gen PR_Japan_share = KJ_PR_ImportQuantity_Dollar/KW_PR_ImportQuantity_Dollar
label variable HF_Japan_share "Japan's import share in Korea"
label variable PR_Japan_share "Japan's import share in Korea"
twoway (bar KW_HF_ImportQuantity_Dollar Year, fcolor(%80) lcolor(none%0)) (bar KJ_HF_ImportQuantity_Dollar Year, lcolor(none%0) lwidth(none)) (line HF_Japan_share Year, yaxis(2)) if Year>2010, ytitle(Import amount) yscale(range(0 1)) ylabel(1e+8(.5e+8)1.5e+8) ytitle(Import share of Japan, axis(2)) xtitle(Year) xlabel(2011(2)2020 2020) xline(2019, lwidth(thick) lpattern(dash) lcolor(black) noextend)
graph export HF_Japan_Share.png, replace
twoway (bar KW_PR_ImportQuantity_Dollar Year, fcolor(%80) lcolor(none%0)) (bar KJ_PR_ImportQuantity_Dollar Year, lcolor(none%0) lwidth(none)) (line PR_Japan_share Year, yaxis(2)) if Year>2010, ytitle(Import amount) yscale(range(0 1)) ylabel(1e+8(2e+8)5.5e+8) xtitle(Year) xline(2019, lwidth(thick) lpattern(dash) lcolor(black) extend) xlabel(2011(2)2020 2020) ylabel(0.8(0.05).99, ax(2))
graph export PR_Japan_Share.png, replace
// cd ..
// cd input
// use DramNANDprices, clear
// twoway (tsline DramDD48GbDram) (tsline NAND128GbMLC)
// cd ..
// cd output
// graph export historic_DramFlash_prices.png, replace
