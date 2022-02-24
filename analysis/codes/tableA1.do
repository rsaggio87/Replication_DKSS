
local dir_tmp_aux = "../build/src/RAKIM/"
*********************************************************************************

*			SUMMARIZE DIFFERENCES ACROSS SAMPLES // FIRM LEVEL

*********************************************************************************
if 0 == 1{
*load baseline file that has sector codes
	use ateco07_2_calc matrx MATRX_MADRE CF_PERSG ID_PERSF codice_contri DIP* year if year>=2005 & year<=2015 using "../build/src/matr_bitalia_1990_2013NEW_INVIND",clear	
	gen mother = 0
	replace mother = 1 if matrx==MATRX_MADRE
	gen id_impresa = CF_PERSG
	tostring ID_PERSF,replace
	replace id_impresa= ID_PERSF if CF_PERSG==""
	destring ateco07_2_calc,replace
	gen industry_code_aux=ateco07_2_calc if mother == 1
	bys id_impresa: egen industry_code=mode(industry_code_aux),maxmode //sometimes firm might switch ateco code, keep modal. around 5% do that.
	egen EMPL_aux = rowmax(DIP1 DIP2 DIP3 DIP4 DIP5 DIP6 DIP7 DIP8 DIP9 DIP10 DIP11 DIP12)
	bys id_impresa year: egen EMPL = total(EMPL_aux) // sum across plants
	gen log_EMPL = log(EMPL)
	collapse (mean) EMPL log_EMPL (mean) industry_code, by(id_impresa) // average across years.
	gen log_avg_EMPL=log(EMPL)
 	save "`dir_tmp_aux'/firms_universe",replace 
}
	
if 1 == 1 {	
	use "`dir_tmp_aux'/firms_universe",replace 
	gen settore = 0
	replace settore = 1 if industry_code == 47 | industry_code == 46 | industry_code == 45  // Retail
	replace settore = 2 if industry_code == 41 | industry_code== 43  // Construction
	replace settore = 3 if industry_code == 56 | industry_code == 55 // Restaurants - Hotels
	replace settore = 4 if industry_code == 96 // Hairdressing/Care Center/
	replace settore = 5 if industry_code == 69 // Lawyers
	replace settore = 6 if industry_code == 33 | industry_code == 10 | industry_code==14 | industry_code==28 | industry_code==25 | industry_code==23 | industry_code==27 // Manufacturing/Reparaining
	replace settore = 7 if industry_code == 49 | industry_code == 50 | industry_code == 51  | industry_code == 52 // Transportation
	replace settore = 8 if industry_code == 81 | industry_code == 80 // Cleaning+Security+
	replace settore = 9 if industry_code == 78 // temp agencies
	replace settore = 10 if industry_code == 70 | industry_code == 62 | industry_code == 82   // Management/Consulting/Software
	replace settore = 11 if industry_code == 64 // finance
	replace settore = 12 if industry_code == 85  | industry_code == 88 | industry_code == 86 // Education and Health
	label define industrie 1 "Retail" 2 "Construction" 3 "Restaurants" 4 "Hairdressing" 5 "Lawyers" 6 "Manufacturing" 7 "Transportation" 8 "Cleaning" 9 "Temp" 10 "Management" 11 "Finance " 12 "Education and Health"
	label values settore industrie
	keep settore id_impresa EMPL log_EMPL log_avg_EMPL
	sum EMPL,d
	sum log_EMPL,d
	sum log_avg_EMPL,d
	tab settore [aw=EMPL]
	save "`dir_tmp_aux'/firms_universe_upd",replace
}	 

*now match with UNIVERSE INPS-INVIND
if 1 == 1{
	use year id_impresa if year>=2005 using "../build/src/RAKIM/Q_CLEANED_SPELL_GIANT_INVIND_data_1990_2016",clear
	collapse year, by(id_impresa)
	drop year
	save "`dir_tmp_aux'/INPS_INVIND_firms",replace 
}
if 1 == 1{
	use "`dir_tmp_aux'/INPS_INVIND_firms",replace 	
	merge 1:1 id_impresa using "`dir_tmp_aux'/firms_universe_upd", gen(match_universe_INPS_INVIND)
	keep if match_universe_INPS_INVIND==3
	sum EMPL,d
	sum log_EMPL,d
	sum log_avg_EMPL,d
	tab settore [aw=EMPL]
	save "`dir_tmp_aux'/firms_INPS_INVIND",replace
}
	
*now match with my estimation sample
if 1 == 1{
	import delimited "`dir_tmp_aux'/POOLED_SAMPLE__INVIND_data_1990_2016_2005_2016_RS_build.csv", encoding(ISO-8859-1) clear
	rename v1 log_dailywages
	rename v3 firmid
	collapse log_dailywages, by(firmid)
	drop  log_dailywages
	merge 1:m firmid using "`dir_tmp_aux'/DAVE_BUILD_RAKIM_APPENDED100_INVIND_data_1990_2016_2005_2016_RS_build", keepusing(id_impresa) gen(cross_walk)
	keep if cross_walk == 3
	merge m:1 id_impresa using "`dir_tmp_aux'/firms_universe_upd", gen(my_sample)
	keep if my_sample == 3
	collapse EMPL log_EMPL log_avg_EMPL settore, by(id_impresa)
	sum EMPL,d
	sum log_EMPL,d
	sum log_avg_EMPL,d
	tab settore [aw=EMPL]
	save "`dir_tmp_aux'/firms_INPS_INVIND",replace
}
	
