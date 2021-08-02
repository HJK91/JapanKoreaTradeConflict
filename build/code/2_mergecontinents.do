// Merge tade across continents
use "KoreaNorthAmericaHFtrade.dta", clear
append using "KoreaAsiaHFtrade.dta"
append using "KoreaEuropeHFtrade.dta"
drop if missing(HS10)
collapse (sum) KW_ExportQuantity_Dollar (sum) KW_ExportQuantity_Kg (sum) KW_ImportQuantity_Dollar (sum) KW_ImportQuantity_Kg, by(Year)
merge 1:1 Year using "KoreaJapanHFtrade.dta", nogenerate
save "KoreaHFtrade.dta", replace

use "KoreaNorthAmericaPRtrade.dta", clear
append using "KoreaAsiaPRtrade.dta"
append using "KoreaEuropePRtrade.dta"
drop if missing(HS10)
collapse (sum) KW_ExportQuantity_Dollar (sum) KW_ExportQuantity_Kg (sum) KW_ImportQuantity_Dollar (sum) KW_ImportQuantity_Kg, by(Year)
merge 1:1 Year using "KoreaJapanPRtrade.dta", nogenerate
save "KoreaPRtrade.dta", replace

use "KoreaNorthAmericaPItrade.dta", clear
append using "KoreaAsiaPItrade.dta"
append using "KoreaEuropePItrade.dta"
drop if missing(HS10)
collapse (sum) KW_ExportQuantity_Dollar (sum) KW_ExportQuantity_Kg (sum) KW_ImportQuantity_Dollar (sum) KW_ImportQuantity_Kg, by(Year)
merge 1:1 Year using "KoreaJapanPItrade.dta", nogenerate
save "KoreaPItrade.dta", replace

// Merge datasets with Japanese trade

