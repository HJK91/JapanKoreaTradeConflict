// input-output quantity relationship
cd "/Users/hj/Documents/GitHub/JapanKoreaTradeConflict/build/input"
use Trade_all.dta, clear
keep if Item == 2811111000 | Item == 3707901010 |  Item == 8542321010 | Item == 8542321030
collapse (sum) Export_dollar Export_weight Import_dollar Import_weight export_quantity, by(Item monthyear)
reshape wide Export_dollar Export_weight Import_dollar Import_weight export_quantity, i(monthyear) j(Item)

