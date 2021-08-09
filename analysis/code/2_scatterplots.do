
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

cd ..
cd input
use DramNANDprices, clear
twoway (tsline DramDD48GbDram) (tsline NAND128GbMLC)
cd ..
cd output
graph export historic_DramFlash_prices.png, replace
