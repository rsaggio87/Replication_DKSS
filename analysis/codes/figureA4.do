*read the global
local dir_tmp_aux 		= "../build/src/RAKIM/"
local J_J_STRINGA 		= "$JJSTRINGA"	
local coorte	  		= 2

if `coorte' == 1{
local CONDA "id>0"
loca  nameRESULTS "start_in_2005"
local var_using "year"
local LABELLAMI "Year"
}
if `coorte' == 2{
local CONDA "birth>=1978 & birth<=1980"
loca  nameRESULTS "start_in_2005_1978_1982"
local var_using "age"
local LABELLAMI "Age"
}
	
**************************************************************************

*	SECTION 1: LOAD THE SPELL LEVEL DATA
	
**************************************************************************
if 1 == 1{	
*LOAD CLEANED SPELL DATA
		use id year id_impresa id_azienda imponibile log_dailywages age sesso qqq real_age_at_entry if year>=2005 using "../build/src/RAKIM/Q_CLEANED_SPELL_GIANT_INVIND_data_1990_2016", replace
		
*COLLAPSE TO PERSON-YEAR
		gsort id year -imponibile
		by id year: gen top_job=_n
		keep if top_job == 1
		drop top_job		

*QUANTILES OF THE DATA
		gen 	   female 			 = 0
		replace    female 			 = 1 		 	 if sesso == "F"
		
*SAVE THIS AUXILIARY FILE
		save "`dir_tmp_aux'/Y_CLEANED_SPELL_GIANT_INVIND_data_1990_2016", replace		
}		
**************************************************************************

*	SECTION 2: AGE PROFILE
	
**************************************************************************	
if 1 == 1 {
*LOAD CLEANED SPELL DATA
		use "`dir_tmp_aux'/Y_CLEANED_SPELL_GIANT_INVIND_data_1990_2016", replace

*ENTERED IN 2005		
		gen 		sel 		= 0 
		replace		sel			= 1									if age == real_age_at_entry & year == 2005
		bys id: egen SEL 		= mean(sel)
		keep if SEL>0			
		
*BIRTH
		capture drop birth
		gen birth = year-age
		keep if `CONDA'		

*GAPS
		gen 		female_wage = log_dailywages 					if female 	== 1 
		gen 		male_wage 	= log_dailywages 					if female   == 0
		
		sum male_wage female_wage		

*COLLAPSE
		collapse female_wage male_wage, by(`var_using')
		gen wage_gap_py = 		male_wage-female_wage
		
*BRING INFO ON HIRING WAGES
		preserve
		use "tables/oaxaca_decom`J_J_STRINGA'`nameRESULTS'",replace
		keep if collapse_var == "`var_using'" & rows == "_1"
		list
		keep wage_gap male_wage female_wage `var_using'
		destring `var_using',replace
		rename wage_gap wage_gap_hiring
		rename male_wage 	 male_wage_hiring
		rename female_wage   female_wage_hiring
		save "figures/hiring_wage`nameRESULTS'",replace
		restore		
		merge 1:1 `var_using' using "figures/hiring_wage`nameRESULTS'"
		
*Do the Baseline Diffs
if `coorte' == 2{
		drop if age >35 
		}
		twoway (connected wage_gap_py `var_using',    lcolor(ltblue)  msize(medlarge) msymbol(square)   mfcolor(ltblue)   legend(label(1 "Log Wage Gender Gap"))) /// 	
		 (connected wage_gap_hiring `var_using',    lcolor(gold)  msize(medlarge) msymbol(square)   mfcolor(gold)   legend(label(2 "Log Hiring Wage Gender Gap"))), /// 		
			   legend(on) ylabel(0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.80) xtitle("Age") ytitle("Log Wage") graphregion(color(white))   bgcolor(white)  legend(ring(3) position(6) rows(3))
		graph export "figures/py_AGE_BASELINE`nameRESULTS'.pdf",replace
		
*Do Levels Diffs
		twoway (connected male_wage `var_using',    lcolor(dknavy)   msize(medlarge)  msymbol(square)   mfcolor(dknavy)   legend(label(1 "Log Wage Men"))) /// 	
		 (connected female_wage `var_using',    lcolor(pink)  msize(medlarge)  msymbol(square)   mfcolor(pink)   legend(label(2 "Log Wage Women"))) /// 
		(connected male_wage_hiring `var_using',    lcolor(dknavy)  msize(medlarge)   msymbol(circle)   mfcolor(dknavy)   legend(label(3 "Log Hiring Wage Men"))) /// 	
		 (connected female_wage_hiring `var_using',    lcolor(pink)  msize(medlarge)  msymbol(circle)   mfcolor(pink)   legend(label(4 "Log Hiring Wage Women"))), /// 		
		legend(on) xtitle("`LABELLAMI'") ytitle("Log Wage") graphregion(color(white))   bgcolor(white)  legend(ring(3) position(6) rows(3)) legend(ring(0) position(11) rows(3)) ylabel(3.8(0.2)4.8)
		graph export "figures/py_levels_AGE_BASELINE`nameRESULTS'.pdf",replace					

}								 	    
**************************************************************************

*	DO THE PLOT THAT NICOLE WAS SUGGESTING
	
**************************************************************************	
if 1 == 1 {
*LOAD CLEANED SPELL DATA
		use if birth>=1965 & birth<=1970 using "`dir_tmp_aux'/Y_CLEANED_SPELL_GIANT_INVIND_data_1990_2016_really", replace

*age at entry
		bys id: egen age_at_entry = min(age)
		keep if age_at_entry == 20
		

*Count experience
		xtset id year
		by id (year): gen experience = _n
		gen experience_male=experience if female == 0
		gen experience_female=experience if female == 1
		
		sum experience if age == 18

*COLLAPSE
		collapse female_wage* male_wage* experience*, by(age_group)
		replace age_group = age_group+1
		drop if age_group == .
			
*Do Levels Diffs
		twoway (connected male_wage age_group,    lcolor(dknavy)   msize(medlarge)  msymbol(square)   mfcolor(dknavy)   legend(label(1 "Log Wage Men"))) /// 	
		 (connected female_wage age_group,    lcolor(pink)  msize(medlarge)  msymbol(square)   mfcolor(pink)   legend(label(2 "Log Wage Women"))) /// 
		(connected male_wage_hiring age_group,    lcolor(dknavy)  msize(medlarge)   msymbol(circle)   mfcolor(dknavy)   legend(label(3 "Log Hiring Wage Men"))) /// 	
		 (connected female_wage_hiring age_group,    lcolor(pink)  msize(medlarge)  msymbol(circle)   mfcolor(pink)   legend(label(4 "Log Hiring Wage Women"))), /// 		
			   legend(on) xtitle("Age") ytitle("Log Wage") graphregion(color(white))   bgcolor(white)  legend(ring(3) position(6) rows(3)) xlabel(1 "18-21" 2 "22-24" 3 "25-27" 4 "28-30" 5 "31-33" 6 "34-36" 7 "37-39" 8 "40-42" 9 "43-45" 10 "45-49" 11 "50+") legend(ring(0) position(11) rows(3))
		graph export "figures/py_levels_AGE_BASELINE_borth.pdf",replace		
		
*Do Experience
		twoway (connected experience_male age_group,    lcolor(dknavy)   msize(medlarge)  msymbol(square)   mfcolor(dknavy)   legend(label(1 "Men"))) /// 	
		 (connected experience_female age_group,    lcolor(pink)  msize(medlarge)  msymbol(square)   mfcolor(pink)   legend(label(2 "Women"))), /// 
  		 legend(on) xtitle("Age") ytitle("Years of Experience in the Labor Market") graphregion(color(white))   bgcolor(white)  legend(ring(3) position(6) rows(3)) xlabel(1 "18-21" 2 "22-24" 3 "25-27" 4 "28-30" 5 "31-33" 6 "34-36" 7 "37-39" 8 "40-42" 9 "43-45" 10 "45-49" 11 "50+") legend(ring(0) position(11) rows(3))
		graph export "figures/py_levels_AGE_BASELINE_borth_exp.pdf",replace							

}	
