local J_J_STRINGA = "$JJSTRINGA"
local file_type = "COMBINED_SAMPLE_INVIND_data_`J_J_STRINGA'"
local dir_tmp_aux = "../build/src/RAKIM/"	
local coorte	  = 1

if `coorte' == 1{
local CONDA "id>0"
loca  nameRESULTS "start_in_2005"
}
if `coorte' == 2{
local CONDA "birth>=1978 & birth<=1980"
loca  nameRESULTS "start_in_2005_1978_1982"
}
*********************************************************************************

*			CCK TABLE IN RAKIM

*********************************************************************************
if 1  == 1 {
use pe prediction id year date log_dailywages female psi* lambda* age firmid lagfirmid gap reason* tipo_contratto birth real_age_at_entry if `CONDA' using "`dir_tmp_aux'/oaxaca_sample_`file_type'",replace

*just focus on people that entered in 2005
		gen 		sel 		= 0 
		replace		sel			= 1									if age == real_age_at_entry & year == 2005
		bys id: egen SEL 		= mean(sel)
		keep if SEL>0	
		
*calculate the counterfactuals now
		gen 		female_wage = log_dailywages 					if female 	== 1 
		gen 		male_wage 	= log_dailywages 					if female   == 0
		
		gen 		female_pe	= pe 					    		if female 	== 1 
		gen 		male_pe	 	= pe 								if female   == 0
		
		gen 		female_pred	= prediction 					    if female 	== 1 
		gen 		male_pred	= prediction 						if female   == 0
	
		gen			female_psi	= psi_female						if female 	== 1
		gen 		male_psi	= psi_male							if female   == 0
		gen 		psi_male_fem= psi_female						if female 	== 0	
				 					 		
		gen			female_lam	= lambda_female						if female 	== 1
		gen 		male_lam	= lambda_male						if female   == 0
		gen 		lam_male_fem= lambda_female						if female 	== 0
	
		gen		    barga_psi	= psi_male-psi_female				if female   == 0
		gen			barga_lam	= lambda_male-lambda_female			if female 	== 0	

*auxiliaries to run different type of decompositions
		gen 			unoo	  = 1
		gen 			age_group = 0 if  age>=18 & age<=21 
		replace 		age_group = 1 if  age>21  & age<=24
		replace 		age_group = 2 if  age>24  & age<=27
		replace 		age_group = 3 if  age>27  & age<=30
		replace 		age_group = 4 if  age>30  & age<=33
		replace 		age_group = 5 if  age>33  & age<=36
		replace 		age_group = 6 if  age>36  & age<=39
		replace 		age_group = 7 if  age>39  & age<=42
		replace 		age_group = 8 if  age>42  & age<=45
		replace 		age_group = 9 if  age>45  & age<=49
		replace 		age_group = 10 if  age>50
		
		gen 			age_group_short = 0 if  age<=30
		replace 		age_group_short = 1 if  age>30  & age<=40
		replace 		age_group_short = 2 if  age>41 
		
		
		gen 			J_to_J 	  = .
		replace 		J_to_J    = 0 if lagfirmid ==0
		replace 		J_to_J    = 1 if lagfirmid >0
		gen temp = .
		replace temp = 1 if tipo_contratto == "D"
		replace temp = 0 if tipo_contratto == "I"
		sort id date
		
		bys id (date): gen NUMERO_LAVORI = _N
		bys id (date): gen n_job_aux = _n
		
*localize
		forval ss = 1(1)7{
		preserve
		if `ss' == 1 {
		local var_coll = "unoo"
		}
		
		if `ss' == 2 {
		local var_coll = "age_group_short"
		}
		
		if `ss' == 3 {
		local var_coll = "J_to_J"
		
		}
		
		if `ss' == 4 {
		local var_coll = "temp"
		}
		
		if `ss' == 5 {
		gen first_job = 0
		replace first_job = 1 if n_job_aux == 1 & lagfirmid<0
		bys id: egen NEW_COHORT = mean(first_job)
		keep if NEW_COHORT > 0
		xtset id date
		bys id: egen MIN_AGE = min(age)
		gen pot_experience = age-MIN_AGE+1
		local var_coll = "pot_experience"
		}
		
		if `ss' == 6 {
		local var_coll = "age"
		}
		
		if `ss' == 7 {
		local var_coll = "year"
		}
		
	
*construct the gaps
		
		collapse female_pred male_pred female_pe male_pe female_wage male_wage female_psi male_psi psi_male_fem female_lam male_lam lam_male_fem barga_psi barga_lam, by(`var_coll')
		drop if `var_coll' == .
		gen predi_gap = 		male_pe-female_pe
*		gen predi_gap =			male_pred-female_pred
		gen wage_gap  = 		male_wage-female_wage
		gen psi_gap	  = 		male_psi-female_psi
		gen lam_gap	  = 		male_lam-female_lam
		gen sort_psi  = 	    psi_male_fem-female_psi
		gen sort_lam  = 	    lam_male_fem-female_lam
		
if `ss' == 7 {
		gen pred_gap_2005 = predi_gap if year == 2005
		gen UNO = 1
		bys UNO: egen pred_gap_norm = mean(pred_gap_2005)
		replace wage_gap = wage_gap + predi_gap -pred_gap_norm
}	
		
*construct the table		
		keep  predi_gap  male_wage female_wage wage_gap male_psi female_psi male_lam female_lam psi_gap lam_gap sort_psi barga_psi sort_lam barga_lam `var_coll'
		order wage_gap male_psi female_psi male_lam female_lam psi_gap lam_gap sort_psi barga_psi sort_lam barga_lam `var_coll' male_wage female_wage
		gen collapse_var = "`var_coll'"
		tab collapse_var if _n == 1
		tab collapse_var if _n == 2
		sort `var_coll'
		
		gen valori = "" 
		
		if `ss' == 1{
		replace valori = "All Jobs"
		}
		
		if `ss' == 2{
		replace valori = "Age <=30" 		if _n == 1
		replace valori = "Age 31-40" 	    if _n == 2
		replace valori = "Age >40" 		    if _n == 3
		}
		
		if `ss' == 3{
		replace valori = "J-to-J transition" 		if _n == 2
		replace valori = "J-to-U-to-J transition"   if _n == 1
		}
		
		if `ss' == 4{
		replace valori = "Permanent Contract"     if _n == 1
		replace valori = "Temporary Contract"     if _n == 2
		}
		
		if `ss' == 5{
				tostring pot_experience, replace
				replace valori = pot_experience
		}
		
		if `ss' == 6{
		tostring age, replace
				replace valori = age
		}
		
		if `ss' == 7{
		tostring year, replace
				replace valori = year
		}
		
		global vars_to_shape male_psi female_psi male_lam female_lam psi_gap lam_gap sort_psi barga_psi sort_lam barga_lam male_wage female_wage
		foreach var in $vars_to_shape{
		rename  `var' `var'_1
		}
		save "tables/oaxaca_decom`ss'",replace
		
		global vars_to_shape male_psi female_psi male_lam female_lam psi_gap lam_gap sort_psi barga_psi sort_lam barga_lam male_wage female_wage
		foreach var in $vars_to_shape{
		replace `var' = (`var'/wage_gap)*100
		rename  `var' `var'_2
		}
		merge 1:1  valori using "tables/oaxaca_decom`ss'", nogen
		reshape long female_wage male_wage male_psi female_psi male_lam female_lam psi_gap lam_gap sort_psi barga_psi sort_lam barga_lam, i(valori) j(rows, string)
		sort `var_coll' rows
		save "tables/oaxaca_decom`ss'",replace
		sum male_psi
		local NNNNN=r(N)+1
		set obs `NNNNN'
		restore
		}

*combine them
		use "tables/oaxaca_decom1",replace
		forval ss = 2(1)7{
		append using "tables/oaxaca_decom`ss'"
		}
		list
		save "tables/oaxaca_decom`J_J_STRINGA'`nameRESULTS'",replace	
		
*shorter version
		use "tables/oaxaca_decom1",replace
		forval ss = 7(1)7{
		append using "tables/oaxaca_decom`ss'"
		}
		list
		destring year,replace
		save "tables/short_oaxaca_decom`J_J_STRINGA'`nameRESULTS'",replace			
}			
*********************************************************************************

*			GRAPH BY # YEAR OF HIRE

*********************************************************************************
if 1 == 1 {
*load sample
use 	"tables/oaxaca_decom`J_J_STRINGA'`nameRESULTS'",replace
		save "figures/full_oaxaca_sample_INVIND`J_J_STRINGA'`nameRESULTS'",replace
		keep if collapse_var == "year"  & rows == "_1"
		gen n_job = year
		destring n_job,replace
		sort n_job
		
*Do the Baseline Diffs
		twoway (connected wage_gap n_job,    lcolor(ltblue)  msize(medlarge) msymbol(square)   mfcolor(ltblue)   legend(label(1 "Adjusted-Hiring Wage Gender Gap"))) ///
			   (connected psi_gap n_job,    lcolor(gold) 	 msize(medlarge) msymbol(circle) mfcolor(gold)     legend(label(2 "Gap in Destination Effects -- {&psi} "))) ///
			   (connected lam_gap n_job ,   lcolor(pink)    msize(medlarge) msymbol(triangle) mfcolor(pink)     legend(label(3 "Gap in Origin Effects -- {&lambda}"))), ///		 		
			   legend(on) xtitle("Year of Hire") legend(order(1 2 3)) ylabel(-0.05 0 0.1 0.2 0.3 0.4) ytitle("Log Hiring Wage") graphregion(color(white))   bgcolor(white)  legend(ring(3) position(6) rows(3))
		graph export "figures/`J_J_STRINGA'_year_of_hire_`nameRESULTS'.pdf",replace
		graph save  "figures/`J_J_STRINGA'_year_of_hire_`nameRESULTS'",replace
			
}
*********************************************************************************

*			GRAPH BY # AGE

*********************************************************************************
if 1 == 1 {

*load sample
use 	"tables/oaxaca_decom`J_J_STRINGA'`nameRESULTS'",replace
		save "figures/full_oaxaca_sample_INVIND`J_J_STRINGA'`nameRESULTS'",replace
		keep if collapse_var == "age" & rows == "_1"
		destring age,replace
		keep if age <=35
		
*Do the Baseline Diffs
		twoway (connected wage_gap age,    lcolor(ltblue)  msize(medlarge) msymbol(square)   mfcolor(ltblue)   legend(label(1 "Adjusted-Hiring Wage Gender Gap"))) ///
			   (connected psi_gap age,    lcolor(gold) 	 msize(medlarge) msymbol(circle) mfcolor(gold)     legend(label(2 "Gap in Destination Effects -- {&psi} "))) ///
			   (connected lam_gap age ,   lcolor(pink)    msize(medlarge) msymbol(triangle) mfcolor(pink)     legend(label(3 "Gap in Origin Effects -- {&lambda}"))), ///	 		
			   legend(on) legend(order(1 2 3 4)) ylabel(-0.05 0 0.1 0.2 0.3 0.4 0.5) xtitle("Age at Hiring") ytitle("Log Hiring Wage") graphregion(color(white))   bgcolor(white)  legend(ring(0) position(11) rows(3)) 	
		graph export "figures/AGE`J_J_STRINGA'`nameRESULTS'.pdf",replace
		graph save  "figures/AGE`J_J_STRINGA'`nameRESULTS'",replace

}
