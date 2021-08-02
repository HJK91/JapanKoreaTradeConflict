// Set the current directory as the folder where all data are contained.

cd "/Users/hj/Documents/GitHub/JapanKoreaTradeConflict/build/input"
// Convert .xlsx files to .dta files

import excel "KoreaNorthAmericaHFtrade.xlsx", firstrow clear
rename ExportQuantityDollar KW_ExportQuantity_Dollar
rename ExportQuantityKg KW_ExportQuantity_Kg
rename ImportQuantityDollar KW_ImportQuantity_Dollar
rename ImportQuantityKg KW_ImportQuantity_Kg
save "KoreaNorthAmericaHFtrade.dta", replace

import excel "KoreaAsiaHFtrade.xlsx", firstrow clear
rename ExportQuantityDollar KW_ExportQuantity_Dollar
rename ExportQuantityKg KW_ExportQuantity_Kg
rename ImportQuantityDollar KW_ImportQuantity_Dollar
rename ImportQuantityKg KW_ImportQuantity_Kg
save "KoreaAsiaHFtrade.dta", replace

import excel "KoreaEuropeHFtrade.xlsx",  firstrow clear
rename ExportQuantityDollar KW_ExportQuantity_Dollar
rename ExportQuantityKg KW_ExportQuantity_Kg
rename ImportQuantityDollar KW_ImportQuantity_Dollar
rename ImportQuantityKg KW_ImportQuantity_Kg
save "KoreaEuropeHFtrade.dta", replace

import excel "KoreaNorthAmericaPRtrade.xlsx",  firstrow clear
rename ExportQuantityDollar KW_ExportQuantity_Dollar
rename ExportQuantityKg KW_ExportQuantity_Kg
rename ImportQuantityDollar KW_ImportQuantity_Dollar
rename ImportQuantityKg KW_ImportQuantity_Kg
save "KoreaNorthAmericaPRtrade.dta", replace

import excel "KoreaAsiaPRtrade.xlsx", firstrow clear
rename ExportQuantityDollar KW_ExportQuantity_Dollar
rename ExportQuantityKg KW_ExportQuantity_Kg
rename ImportQuantityDollar KW_ImportQuantity_Dollar
rename ImportQuantityKg KW_ImportQuantity_Kg
save "KoreaAsiaPRtrade.dta", replace

import excel "KoreaEuropePRtrade.xlsx",  firstrow clear
rename ExportQuantityDollar KW_ExportQuantity_Dollar
rename ExportQuantityKg KW_ExportQuantity_Kg
rename ImportQuantityDollar KW_ImportQuantity_Dollar
rename ImportQuantityKg KW_ImportQuantity_Kg
save "KoreaEuropePRtrade.dta", replace

import excel "KoreaNorthAmericaPItrade.xlsx",  firstrow clear
rename ExportQuantityDollar KW_ExportQuantity_Dollar
rename ExportQuantityKg KW_ExportQuantity_Kg
rename ImportQuantityDollar KW_ImportQuantity_Dollar
rename ImportQuantityKg KW_ImportQuantity_Kg
save "KoreaNorthAmericaPItrade.dta", replace

import excel "KoreaAsiaPItrade.xlsx", firstrow clear
rename ExportQuantityDollar KW_ExportQuantity_Dollar
rename ExportQuantityKg KW_ExportQuantity_Kg
rename ImportQuantityDollar KW_ImportQuantity_Dollar
rename ImportQuantityKg KW_ImportQuantity_Kg
save "KoreaAsiaPItrade.dta", replace

import excel "KoreaEuropePItrade.xlsx",  firstrow clear
rename ExportQuantityDollar KW_ExportQuantity_Dollar
rename ExportQuantityKg KW_ExportQuantity_Kg
rename ImportQuantityDollar KW_ImportQuantity_Dollar
rename ImportQuantityKg KW_ImportQuantity_Kg
save "KoreaEuropePItrade.dta", replace

// Trade with Japan

import excel "KoreaJapanHFtrade.xlsx", firstrow clear
drop if missing(HS10)
rename ExportQuantityDollar KJ_ExportQuantity_Dollar
rename ExportQuantityKg KJ_ExportQuantity_Kg
rename ImportQuantityDollar KJ_ImportQuantity_Dollar
rename ImportQuantityKg KJ_ImportQuantity_Kg
save "KoreaJapanHFtrade.dta", replace

import excel "KoreaJapanPRtrade.xlsx", firstrow clear
drop if missing(HS10)
rename ExportQuantityDollar KJ_ExportQuantity_Dollar
rename ExportQuantityKg KJ_ExportQuantity_Kg
rename ImportQuantityDollar KJ_ImportQuantity_Dollar
rename ImportQuantityKg KJ_ImportQuantity_Kg
save "KoreaJapanPRtrade.dta", replace

import excel "KoreaJapanPItrade.xlsx", firstrow clear
drop if missing(HS10)
rename ExportQuantityDollar KJ_ExportQuantity_Dollar
rename ExportQuantityKg KJ_ExportQuantity_Kg
rename ImportQuantityDollar KJ_ImportQuantity_Dollar
rename ImportQuantityKg KJ_ImportQuantity_Kg
save "KoreaJapanPItrade.dta", replace
