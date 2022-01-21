use DramNANDprices.dta, clear
gen price1 = DramDD48GbDram
gen price2 = NAND128GbMLC
keep MonthYear price1 price2
reshape long price, i(MonthYear) j(item)
gen monthyear = MonthYear
drop MonthYear
gen Item = .
recast double Item
replace Item = 8542321010 if item ==1
replace Item = 8542321030 if item ==2
drop item
save price_long.dta, replace
