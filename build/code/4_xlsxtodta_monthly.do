// Set the current directory as the folder where all data are contained.

cd "/Users/hj/Documents/GitHub/JapanKoreaTradeConflict/build/input"
// Convert .xlsx files to .dta files

foreach x in "HF" "PR"{
	foreach y in "NorthAmerica" "Asia" "Europe" {
		import excel "Korea`y'`x'tradeMonthly.xlsx", firstrow clear
		rename ImportDollar KW_`x'_ImportDollar
		rename ImportKg KW_`x'_ImportKg
		generate MonthYear = ym(Year, Month)
		format MonthYear %tmMonYY
		tsset  MonthYear, monthly
		save "Korea`y'`x'tradeMonthly.dta", replace
	}	
}



// Trade with Japan
foreach x in "HF" "PR" {
	import excel "KoreaJapan`x'tradeMonthly.xlsx", firstrow clear
	rename ImportDollar KJ_`x'_ImportDollar
	rename ImportKg KJ_`x'_ImportKg
	generate MonthYear = ym(Year, Month)
	format MonthYear %tmMonYY
	tsset  MonthYear, monthly
	save "KoreaJapan`x'tradeMonthly.dta", replace
}


// Trade of Dram and NAND Flash
foreach x in "NAND" "Dram" {
	foreach y in "NorthAmerica" "SouthAmerica" "Asia" "MiddleEast" "Oceania" "Africa" "Europe" {
		import excel "Korea`y'`x'tradeMonthly.xlsx", firstrow clear
		rename ExportDollar `x'_ExportDollar
		rename ExportKg `x'_ExportKg
		generate MonthYear = ym(Year, Month)
		format MonthYear %tmMonYY
		tsset  MonthYear, monthly
		save "Korea`y'`x'tradeMonthly.dta", replace
	}
}

// Prices of Dram and NAND
import excel "DramNANDprices.xlsx", firstrow clear
drop if missing(Year)
generate MonthYear = ym(Year, Month)
format MonthYear %tmMonYY
tsset  MonthYear, monthly
save "DramNANDprices.dta", replace

// Merge across continents
foreach x in "HF" "PR" {
	use "KoreaNorthAmerica`x'tradeMonthly.dta", clear
	foreach y in "Asia" "Europe" {
		append using "Korea`y'`x'tradeMonthly.dta"
	}
	collapse (sum) KW_`x'_ImportDollar (sum) KW_`x'_ImportKg, by(MonthYear)
	merge 1:1 MonthYear using "KoreaJapan`x'tradeMonthly.dta", nogenerate
	destring MonthYear, replace
	save "Korea`x'tradeMonthly.dta", replace
}

foreach x in "Dram" "NAND" {
	use "KoreaNorthAmerica`x'tradeMonthly.dta", clear
	foreach y in "SouthAmerica" "Asia" "MiddleEast" "Oceania" "Africa" "Europe"  {
		append using "Korea`y'`x'tradeMonthly.dta"
	}
	collapse (sum) `x'_ExportDollar (sum) `x'_ExportKg, by(MonthYear)
	destring MonthYear, replace
	save "Korea`x'tradeMonthly.dta", replace
}

use "KoreaPRtradeMonthly.dta", clear
foreach x in "HF" "NAND" "Dram"{
	merge 1:1 MonthYear using "Korea`x'tradeMonthly.dta", nogenerate
}

merge 1:1 MonthYear using "DramNANDprices.dta", nogenerate

drop Year Month
rename NAND128GbMLC price_NAND
rename DramDD48GbDram price_Dram
cd ..
cd ..
cd analysis/input
save "TradeData", replace
