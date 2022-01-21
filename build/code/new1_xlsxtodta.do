cd "/Users/hj/Documents/GitHub/JapanKoreaTradeConflict/build/input"
import excel "Trade_all.xlsx", firstrow clear

egen yearly_import=sum(Import_dollar), by (Item Country Year)
egen yearly_export=sum(Export_dollar), by (Item Country Year)
egen ttl_monthly_import=sum(Import_dollar), by (Item Year Month)
egen ttl_monthly_export=sum(Export_dollar), by (Item Year Month)
egen ttl_yearly_import=sum(Import_dollar), by (Item Year)
egen ttl_yearly_export=sum(Export_dollar), by (Item Year)
gen importshare_monthly = Import_dollar/ttl_monthly_import
gen exportshare_monthly = Export_dollar/ttl_monthly_export
keep if Year>2007
gen monthyear = ym(Year,Month)
format monthyear %tm
merge m:1 Item monthyear using "price_long.dta", nogen
gen export_quantity = Export_dollar/price
save "Trade_all.dta", replace

