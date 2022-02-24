********************************************************************************

						*DECLARE AND INSTALL* 	
					
********************************************************************************
set seed 1234
set more off
local small_sample		= 0 //<---------------------------Instead of working with the enormous file take only,say, 1000 workers from each file
local ss				= 2
local total_files_app	= 100
if `small_sample' == 1 {
	local sample_load = "small"
}
local type_data 		"INVIND_data_1990_2016`sample_load'"
local dir_tmp_aux.      "$dirVP/RAKIM/"
********************************************************************************

*	IMPOSE RESTRICTIONS AND EXPORT TO MATLAB
	
*******************************************************************************												 	    
*LOAD
		use "$dirVP/RAKIM/RAKIM_APPENDED`total_files_app'_`type_data'", replace 
		
*DROP JOBS FOR WHICH I HAVE A RIGHT CENSORING PROBLEM (CAN'T REALLY SAY WHETHER LAGFIRM = 0 or not)
		drop if right_cen == 1 // notice that by imposing this condition, there is going to be some individuals (those who had a job in 1990) for which gap will always be non-missing. 		
		
*DROP DUDES OBSERVED WITH ONLY ONE JOB
		distinct id
		bys id: gen TTT = _N
		drop if TTT == 1
		distinct id
		drop TTT
		
*count j-j transitions
		gen 		jj_trans = 0
		replace 	jj_trans = 1 		if 	lagfirmid>0
		
		sum jj_trans if conteggio>1		
		
*redefine lagfirmid
		drop lagfirmid conteggio
		gen lagfirmid = 0
		bys id: gen pre_firm = firmid[_n-1]
		bys id:  gen conteggio = _n
		replace lagfirmid 	  =	-1 				if conteggio==1
		replace lagfirmid 	  =  pre_firm		if reason == 2 & conteggio > 1 | gap==1 & reason == . & conteggio>1
		
		 		
		
*count j-j transitions
		drop jj_trans
		gen 		jj_trans = 0
		replace 	jj_trans = 1 		if 	lagfirmid>0
		
		sum jj_trans if conteggio>1 		
		
**ss==2, only from 2005.
	if `ss' == 1{
			local 	type_data = "`type_data'_RS_BUILD"
	}
	if `ss' == 2 {
			keep if year>=2005
			local type_data "`type_data'_2005_2016_RS_build"
	}
	
*investigate missing in separations
		tab reason, missing
		tab reason if conteggio>1, missing
		tab gap if reason == . & conteggio > 1
		
		save "$dirVP/RAKIM/DAVE_BUILD_RAKIM_APPENDED`total_files_app'_`type_data'", replace

*EXPORT
		export delimited id firmid lagfirmid year log_dailywages age female firmsize using "$dirVP/RAKIM/INVIND_RAKIM_`type_data'", replace novarnames nolabel
	

*EXPORT FIRM-LEVEL FILE WITH INFO ON SIZE AND SECTOR TO BE READ BY MATLAB
		collapse firmid DIP10, by(id_impresa) // this is going to average firm size over the years.
		save "$dirVP/RAKIM/firm_size_aux",replace
		use ateco07_2_calc matrx MATRX_MADRE CF_PERSG ID_PERSF codice_contri DIP* year if year>=2005 & year<=2015 using 	"src/matr_bitalia_1990_2013NEW_INVIND",clear	
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
		save "$dirVP/RAKIM/firms_universe_upd",replace
		use  f"$dirVP/RAKIM/firm_size_aux",replace
		replace DIP10 = ceil(DIP10)
		sum DIP10,d
		merge 1:1 id_impresa using "`dir_tmp_aux'/firms_universe_upd", gen(my_sample) keepusing(settore)
		keep if my_sample == 3
		export delimited firmid DIP10 settore using "`dir_tmp_aux'/FIRM_SIZE.csv", replace novarnames nolabel			
