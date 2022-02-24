local ss				= $it
local dir_tmp_aux 		= "src/RAKIM/"	
**************************************************************************

*	SECTION 1: START BY TAKING WAGES IN YEAR 3 OF THE MATCH
	
**************************************************************************
if 1 == 1 {	
*LOAD THE MONTHLY PANEL
		use id date id_impresa id_azienda log_dailywages year working using "/scratch/public/board/RAKIM/CLEANED_MONTHLY_PANEL_GIANT_`ss'_INVIND_data_1990_2016", replace 
			
*CREATE THE TENURE VARIABLE	
		xtset id date
		bys id id_impresa working (date): gen tenure =_n
		replace tenure = 0 if working == 0			

*SHOW ME SPELLS AT MONTH 36
		bys id id_impresa: egen maxTenure = max(tenure)
		keep if tenure == 36 & working == 1
		tostring id_azienda, replace
		keep id id_impresa id_azienda log_dailywages year
		keep if year>=2005
		drop year
		drop id_azienda
		rename log_dailywages log_dailywages_3
*SAVE				
		save "/scratch/public/leave_out/Third_CLEANED_MONTHLY_PANEL_GIANT_`ss'_INVIND_data_1990_2016", replace 
}
**************************************************************************

*	SECTION 2: APPEND THE FILES CREATED ABOVE
	
**************************************************************************
if  1 == 1 {	
*LOAD THE FIRST ONE
		use "/scratch/public/leave_out/Third_CLEANED_MONTHLY_PANEL_GIANT_1_INVIND_data_1990_2016", replace 
		
		forval pp = 2(1)100 {
		append using "/scratch/public/leave_out/Third_CLEANED_MONTHLY_PANEL_GIANT_`pp'_INVIND_data_1990_2016",
		}
		
*SAVE				
		save "/scratch/public/leave_out/Third_INVIND_data_1990_2016", replace 		
}		
**************************************************************************

*	SECTION 3: NOW CALL THE BASELINE FILE AND MERGE IN 

*	(i) 	THE ESTIMATED EFFECTS FROM DWL
*	(ii)	THE WAGE IN YEAR 3
	
**************************************************************************	
if 1 == 1 {
*firm effects first	
	import delimited "`dir_tmp_aux'/POOLED_SAMPLE__INVIND_data_1990_2016_2005_2016_RS_build.csv", encoding(ISO-8859-1) clear
	rename v3 firmid 
	rename v4 lagfirmid
	rename v6 fe
	rename v7 lag_fe
	rename v8 controls
	rename v19 age2 
	rename v20 age3
	gen year = 2004
	local aux = 1
	forval pp = 9(1)18 {
	replace year = 2004 + `aux' if v`pp' == 1
	local aux = `aux' +1
	}
	replace year = 2015 if year == 2004
	sum year
	collapse fe lag_fe controls, by(firmid lagfirmid age2 age3 year)
	save "`dir_tmp_aux'/identified_firms",replace

*merge in effects	
	use "`dir_tmp_aux'/DAVE_BUILD_RAKIM_APPENDED100_INVIND_data_1990_2016_2005_2016_RS_build",replace	
	gen age2= ((age-40)/40)^2
	gen age3= ((age-40)/40)^3
	merge m:1 firmid lagfirmid age2 age3 year using "`dir_tmp_aux'/identified_firms", gen(identified)
	keep if identified ==3
	
*merge in VA/L
	merge m:1 firmid using "`dir_tmp_aux'/va_file", gen(VA) keepusing(log_VA_L)
	keep if VA == 3	

*reconstruct worker effect
	gen res = log_dailywages-fe-lag_fe-controls
	bys id: egen pe = mean(res)	
	
*now bring wage at time 3	
	merge 1:1 id id_impresa using "/scratch/public/leave_out/Third_INVIND_data_1990_2016", gen(merge_3)
	
*compute quantiles
	drop qqq
	xtile qqq_aux = log_VA_L if fe!=. & merge_3==3, nquantile(100)	
	
*save
	save "`dir_tmp_aux'/Third_eye",replace
	 	
}
		
