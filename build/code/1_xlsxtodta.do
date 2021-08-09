// Set the current directory as the folder where all data are contained.

cd "/Users/hj/Documents/GitHub/JapanKoreaTradeConflict/build/input"
// Convert .xlsx files to .dta files

foreach x in "HF" "PR" "PI" {
	foreach y in "NorthAmerica" "Asia" "Europe" {
		import excel "Korea`y'`x'trade.xlsx", firstrow clear
		drop if missing(HS10)
		rename ExportQuantityDollar KW_`x'_ExportQuantity_Dollar
		rename ExportQuantityKg KW_`x'_ExportQuantity_Kg
		rename ImportQuantityDollar KW_`x'_ImportQuantity_Dollar
		rename ImportQuantityKg KW_`x'_ImportQuantity_Kg
		save "Korea`y'`x'trade.dta", replace
	}	
}



// Trade with Japan
foreach x in "HF" "PR" "PI" {
	import excel "KoreaJapan`x'trade.xlsx", firstrow clear
	drop if missing(HS10)
	rename ExportQuantityDollar KJ_`x'_ExportQuantity_Dollar
	rename ExportQuantityKg KJ_`x'_ExportQuantity_Kg
	rename ImportQuantityDollar KJ_`x'_ImportQuantity_Dollar
	rename ImportQuantityKg KJ_`x'_ImportQuantity_Kg
	save "KoreaJapan`x'trade.dta", replace
}


// Trade of Dram and NAND Flash
foreach x in "Flash" "Dram" {
	foreach y in "Others" "NorthAmerica" "SouthAmerica" "Asia" "MiddleEast" "Oceania" "Africa" "Europe" {
		import excel "Korea`y'`x'trade.xlsx", firstrow clear
		drop if missing(HS10)
		rename ExportQuantityDollar `x'_ExportQuantity_Dollar
		rename ExportQuantityKg `x'_ExportQuantity_Kg
		rename ImportQuantityDollar `x'_ImportQuantity_Dollar
		rename ImportQuantityKg `x'_ImportQuantity_Kg
		save "Korea`y'`x'trade.dta", replace
	}
}

// Priced of Dram and NAND
import excel "DramNANDprices.xlsx", firstrow clear
drop if missing(Year)
generate MonthYear = ym(Year, Month)
format MonthYear %tmMonYY
tsset  MonthYear, monthly
cd ..
cd ..
cd analysis/input
save "DramNANDprices.dta", replace
