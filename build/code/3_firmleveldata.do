// Set the current directory as the folder where all data are contained.

cd "/Users/hj/Documents/GitHub/JapanKoreaTradeConflict/build/input"
import delimited CE_MiningManufacturing_2015.csv, clear 
destring workers payroll production_shipment_amount operating_cost material_cost operating_profit production value_added, replace force
egen total_VA = sum(value_added)
keep if industry_large == "C" // manufacturing sector
keep if industry_medium == 26 //manufacturing_electric equipment
keep if industry_small == 1 // manufacturing_electric equipment_semiconductor
egen semiconductor_VA = sum(value_added)
egen semiconductor_matrial_cost = sum(material_cost)
