cd "/Users/hj/Documents/GitHub/JapanKoreaTradeConflict/analysis/input"
use TradeData, clear


gen Z_Dram  = log(Dram_ExportKg * price_Dram / Dram_ExportDollar)
gen Z_NAND = log(NAND_ExportKg * price_NAND / NAND_ExportDollar)
foreach x in "7" "8" "9" "10" "11" "12"{
	gen dispute_dummy`x' = MonthYear==tm(2019m`x')
}

foreach x in "1" "2" "3" "4" "5"{
	gen dispute_dummy`x' = MonthYear==tm(2020m`x')
}
gen t = _n
gen t_sq = t^2
gen dispute_dummy`x' = MonthYear>=tm(2019m8) 
gen covid_dummy = MonthYear>=tm(2020m5)
gen Z_Dram_diff = Z_Dram[_n]-Z_Dram[_n-1]
gen Z_NAND_diff = Z_NAND[_n]-Z_NAND[_n-1]
reg Z_Dram_diff t t_sq covid_dummy dispute_dummy
reg Z_NAND_diff covid_dummy dispute_dummy
