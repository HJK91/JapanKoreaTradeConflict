// Set the current directory as the folder where all data are contained.
cd "/Users/hj/Documents/GitHub/JapanKoreaTradeConflict/analysis/input"
use KoreaHFtrade, clear

// merge across trade
foreach x in "PR" "PI" "Dram" "Flash" {
	merge 1:1 Year using "Korea`x'trade.dta", nogenerate noreport
}

drop HS10

//Change labels
foreach x of varlist KW_* {
	foreach y in "HF" "PR" "PI" {
		label variable KW_`y'_ExportQuantity_Dollar "`y' exports to all countries"
		label variable KW_`y'_ImportQuantity_Dollar "`y' imports from all countries"
		label variable KW_`y'_ExportQuantity_Kg "`y' exports to all countries (in kg)"
		label variable KW_`y'_ImportQuantity_Kg "`y' imports from all countries (in kg)"
	}
}

foreach x of varlist KJ_* {
	foreach y in "HF" "PR" "PI" {
		label variable KJ_`y'_ExportQuantity_Dollar "`y' exports to Japan"
		label variable KJ_`y'_ImportQuantity_Dollar "`y' imports from Japan"
		label variable KJ_`y'_ExportQuantity_Kg "`y' exports to Japan (in kg)"
		label variable KJ_`y'_ImportQuantity_Kg "`y' imports from  Japan (in kg)"
	}
}

foreach x in "Dram" "Flash" {
	foreach y in "Export" "Import"{
		foreach z in "Dollar" "Kg"{
			label variable `x'_`y'Quantity_`z' "`x' `y' from Korea measured in `z'"
		}
	}
}


// generate share variables
foreach x in "HF" "PR" "PI" {
	generate Japan_`x'_share = KJ_`x'_ImportQuantity_Dollar/KW_`x'_ImportQuantity_Dollar
	label variable Japan_`x'_share "Japan's share of `x' importing to Korea"
}
