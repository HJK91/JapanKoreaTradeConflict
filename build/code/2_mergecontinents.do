cd "/Users/hj/Documents/GitHub/JapanKoreaTradeConflict/build/input"

// Merge tade across continents
foreach x in "HF" "PR" "PI" {
	use "KoreaNorthAmerica`x'trade.dta", clear
	foreach y in "Asia" "Europe" {
		append using "Korea`y'`x'trade.dta"
		drop if missing(HS10)
	}
	collapse (sum) KW_`x'_ExportQuantity_Dollar (sum) KW_`x'_ExportQuantity_Kg (sum) KW_`x'_ImportQuantity_Dollar (sum) KW_`x'_ImportQuantity_Kg, by(Year)
	merge 1:1 Year using "KoreaJapan`x'trade.dta", nogenerate
	destring Year, replace
	drop if Year>2020
	cd ..
	cd output
	save "Korea`x'trade.dta", replace
	cd .. 
	cd ..
	cd analysis/input
	save "Korea`x'trade.dta", replace
	cd ..
	cd ..
	cd build/input
}


//Trade of chips


foreach x in "Flash" "Dram" {
	use "KoreaOthers`x'trade.dta", clear
	foreach y in "NorthAmerica" "SouthAmerica" "Asia" "MiddleEast" "Oceania" "Africa" {
		append using "Korea`y'`x'trade.dta"
	}
	collapse (sum) `x'_ExportQuantity_Dollar (sum) `x'_ExportQuantity_Kg (sum) `x'_ImportQuantity_Dollar (sum) `x'_ImportQuantity_Kg, by(Year)
	destring Year, replace
	drop if Year>2020
	save "Korea`x'trade.dta", replace
	cd ..
	cd output
	save "Korea`x'trade.dta", replace
	cd .. 
	cd ..
	cd analysis/input
	save "Korea`x'trade.dta", replace
	cd ..
	cd ..
	cd build/input
}
