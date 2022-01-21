cd "/Users/hj/Documents/GitHub/JapanKoreaTradeConflict/build/input"
// plot for HF
use Trade_all.dta, clear
keep if Item == 2811111000
keep Country Year monthyear importshare_monthly Import_dollar 
reshape wide Import_dollar importshare_monthly, i(monthyear) j(Country) string
twoway (line importshare_monthlyTaiwan monthyear, sort) (line importshare_monthlyJapan monthyear, sort) ///
		(line importshare_monthlyChina monthyear, sort) (line importshare_monthlyUSA monthyear, sort) if  Year >2015, xline(714)
		
// plot for PR
use Trade_all.dta, clear
keep if Item == 3707901010
keep Country Year monthyear importshare_monthly Import_dollar 
reshape wide Import_dollar importshare_monthly, i(monthyear) j(Country) string
twoway (line importshare_monthlyTaiwan monthyear, sort) (line importshare_monthlyJapan monthyear, sort) ///
		(line importshare_monthlyChina monthyear, sort) (line importshare_monthlyUSA monthyear, sort) ///
		(line importshare_monthlyBelgium monthyear, sort) (line importshare_monthlyCanada monthyear, sort) ///
		(line importshare_monthlyGermany monthyear, sort)  (line importshare_monthlySingapore monthyear, sort) if  Year >2015, xline(714)
		
// plot for Dram
use Trade_all.dta, clear
keep if Item == 8542321010
keep Country Year monthyear exportshare_monthly Export_dollar 
reshape wide Export_dollar exportshare_monthly, i(monthyear) j(Country) string
twoway (line Export_dollarTaiwan monthyear, sort) (line Export_dollarJapan monthyear, sort) ///
		(line Export_dollarChina monthyear, sort) (line Export_dollarUSA monthyear, sort) if  Year >2015, xline(714) //  display ym(2019,7) = 714

// plot for NAND
use Trade_all.dta, clear
keep if Item == 8542321030
keep Country Year monthyear exportshare_monthly Export_dollar 
reshape wide Export_dollar exportshare_monthly, i(monthyear) j(Country) string
twoway (line Export_dollarTaiwan monthyear, sort) (line Export_dollarJapan monthyear, sort) ///
		(line Export_dollarChina monthyear, sort) (line Export_dollarUSA monthyear, sort) ///
		(line Export_dollarHongKong monthyear, sort) (line Export_dollarSingapore monthyear, sort) ///
		(line Export_dollarThailand monthyear, sort) if  Year >2015, xline(714) 
		
// plot for Silicon wafer
use Trade_all.dta, clear
keep if Item == 3818001000
keep Country Year monthyear importshare_monthly Import_dollar 
reshape wide Import_dollar importshare_monthly, i(monthyear) j(Country) string
twoway (line importshare_monthlyTaiwan monthyear, sort) (line importshare_monthlyJapan monthyear, sort) ///
		(line importshare_monthlyChina monthyear, sort) (line importshare_monthlyUSA monthyear, sort) ///
		(line importshare_monthlyBelgium monthyear, sort) (line importshare_monthlyCanada monthyear, sort) ///
		(line importshare_monthlyGermany monthyear, sort)  (line importshare_monthlySingapore monthyear, sort) if  Year >2015, xline(714)
		
		
// plot for export quantity: Dram
use Trade_all.dta, clear
keep if Item == 8542321010 & Year>2016
keep Country Year monthyear export_quantity
reshape wide export_quantity, i(monthyear) j(Country) string
twoway (line export_quantityTaiwan monthyear, sort) (line export_quantityJapan monthyear, sort) ///
		(line export_quantityChina monthyear, sort) (line export_quantityUSA monthyear, sort) ///
		(line export_quantityBelgium monthyear, sort) (line export_quantityCanada monthyear, sort) ///
		(line export_quantityGermany monthyear, sort)  (line export_quantitySingapore monthyear, sort) if  Year >2015, xline(714)
