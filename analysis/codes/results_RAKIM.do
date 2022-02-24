********************************************************************************
*INPUT:  FILE CREATED FROM ESTIM_RAKIM + FILE WITH PSI AND LAMBDA
********************************************************************************
set more off
local dir_tmp_aux "/scratch/public/leave_out"	
cd 	  "/accounts/projects/decomposition/RAKIM"
set scheme plotplain
*local dir_tmp_aux "scratch"	
capture log close
log using "logs/results_RAKIM.log", replace

********************************************************************************

*					TELL ME WHETHER YOU WANT JJ or JUJ

********************************************************************************

local jj_flag = 1

if `jj_flag' == 1{
	local J_J_STRINGA = "1990_2016_2005_2016"	
}
********************************************************************************

*					BUILD AUXILIARY FILE FOR ANALYSIS

*********************************************************************************
if 0 == 1{
*run the auxiliary do file
	global file_type = "POOLED_SAMPLE_INVIND_data_`J_J_STRINGA'"
	do "codes/firm_effects_PVR"
	local file_type = "$file_type"
}	
********************************************************************************

*					FIRM EFFECTS ONTO FIRM SIZE

*********************************************************************************
*Load the file created by firm_effects_PVR
	use "`dir_tmp_aux'/fe_`file_type'",replace		    
   
*Estimate the slope
	reg psi log_firm_size
	local slope_psi = round(_b[log_firm_size],0.001)
	
	reg lambda log_firm_size
	local slope_lambda = round(_b[log_firm_size],0.001)
	
*visualize the relationship
	xtile qqq=log_firm_size, nquantiles(100)
	
*collapse
	collapse psi lambda log_firm_size, by(qqq)
	
*plot
	twoway (scatter psi log_firm_size, yaxis(1) msymbol(square) mcolor(ltblue) legend(label(1 "{&psi} -- Current Firm Effects"))) (lfit psi log_firm_size, yaxis(1) lcolor(ltblue) lpattern(solid) legend(label(2 "Regression Line"))) (scatter lambda log_firm_size, yaxis(1) msymbol(square) mcolor(gold) legend(label(3 "{&lambda} -- Lag Firm Effects"))) (lfit lambda log_firm_size, yaxis(1) lcolor(gold) lpattern(solid) legend(label(4 "Regression Line"))), ytitle("Firm Effects") xtitle("Log Firm Size") note("Regression slope for {&psi}: `slope_psi' " "Regression slope for {&lambda}: `slope_lambda'",  pos(11) place(s)) plotregion(margin(zero)) legend(order(1 3)) legend(ring(2) position(6) rows(1))
	graph export "figures/figure_size`J_J_STRINGA'.pdf",replace
	graph save  "figures/figure_size`J_J_STRINGA'",replace
	
*********************************************************************************

*							DAKM VS AKM

*********************************************************************************
*Load the file created by firm_effects_PVR
	use "`dir_tmp_aux'/fe_`file_type'",replace	

*evaluate
	reg psi psi_AKM	
	local slopEE = round(_b[psi_AKM],0.001)
	 	
*find the quantiles
	xtile qqq=psi_AKM, nquantiles(100)

*collapse and plot
	collapse psi psi_AKM, by(qqq)

*plot
	twoway (scatter psi psi_AKM, yaxis(1) msymbol(square) mcolor(cranberry) legend(label(1 "{&psi} -- Current Firm Effects"))) (lfit psi psi_AKM, yaxis(1) lcolor(cranberry) lpattern(solid) legend(label(2 "Regression Line"))), legend(off)  ytitle("Firm Effects - DAKM") xtitle("Firm Effects - AKM") note("Regression slope: `slopEE'",  pos(11) place(s)) xlabel(-0.6(0.2)0.6) ylabel(-0.6(0.2)0.6) plotregion(margin(zero))
	graph export "figures/figure_AKM_DAKM_`J_J_STRINGA'.pdf",replace
	graph save  "figures/figure_AKM_DAKM_`J_J_STRINGA'",replace
}	































*********************************************************************************
*	CALL THE POOLED ESTIMATES AND SPLIT THEM ACROSS PSI AND LAMBDA
*********************************************************************************
	use "`dir_tmp_aux'/PVR_POOLED_SAMPLE_`J_J_STRINGA'_ESTIMATES_WITH_CROSS_WALK",replace 
	
	preserve
	keep psi* identified_psi* firmid
	save "`dir_tmp_aux'/PSI_PVR_POOLED_SAMPLE_`J_J_STRINGA'_ESTIMATES_WITH_CROSS_WALK",replace 
	restore
	
	preserve
	keep lambda* identified_lambda* firmid
	rename firmid firmid_lag
	save "`dir_tmp_aux'/LAMBDA_PVR_POOLED_SAMPLE_`J_J_STRINGA'_ESTIMATES_WITH_CROSS_WALK",replace 
	restore
	
*********************************************************************************
*	CONSTRUCT THE T = 4 AUXILIARY DATASET NEEDED FOR THE LAST TWO GRAPHS
*********************************************************************************
if 0 == 1 {
*save the ids and group identifiers
	use "`dir_tmp_aux'/RAKIM_GROUP1",replace
	gen group_type = 1
	append using "`dir_tmp_aux'/RAKIM_GROUP2"
	replace group_type = 2 if group_type==.
	append using "`dir_tmp_aux'/RAKIM_GROUP3"
	replace group_type = 3 if group_type==.
	tab group_type
	collapse group_type, by(id)
	tab group_type
	save "`dir_tmp_aux'/RAKIM_GROUP_collapsed",replace

*start with the micro-data for everyone
	 use "`dir_tmp_aux'/RAKIM_start",replace
	 
*work with the right ids
	 merge m:1 id using "`dir_tmp_aux'/RAKIM_GROUP_collapsed", gen(rakim)
	 keep if rakim == 3
	 
*tsfill, balanced the data
	bys id: egen age_max = max(age)
	append using "`dir_tmp_aux'/aux_file"	 	 
	xtset id year
	drop age type_worker
	
	tsfill, full
	replace unemployed			= 1 			if unemployed==.
	replace log_dailywages		= 0 			if unemployed == 1
	replace firmidnew 			= 0 			if unemployed == 1
	replace provN 				= 0 			if unemployed == 1
	replace firmid				= 0 			if unemployed == 1
	replace firmsize			= 0 			if unemployed == 1
	
	bysort id: egen aux =  mean(birth)
	replace birth=aux
	drop aux
	
	bysort id: egen aux =  mean(age_max)
	replace age_max=aux
	drop aux
	
	bysort id: egen aux =  mean(group_type)
	replace group_type=aux
	drop aux
	
	bysort id: egen aux =  mean(female)
	replace female=aux
	drop aux
	
	bysort id: egen aux =  mean(real_age_at_entry)
	replace real_age_at_entry=aux
	drop aux
	
	drop if id == 0
	gen     age 					= year-birth
	egen    type_worker				= group(year real_age_at_entry female)
	gen	    differenze 				= age-real_age_at_entry
	keep if differenze 				>= -1 // keep the year before entry in the labor market
	replace firmidnew 				= -1  													if differenze == -1
	replace firmid 					= -1  													if differenze == -1
	replace MIN_YEAR 				= year 													if differenze == -1
	
*I don't need your observations once you have retired
	drop if age > age_max
	
*now the very important stuff	
	keep if cum_tag <=3 | unemployed == 1 // keep the first three jobs only and the unemployment state
	bys id: egen min_year_unemployed = min(year) if firmid == 0
	keep if year == MIN_YEAR | year == min_year_unemployed & group_type == 2
	bys id: gen TTT=_N
	tab TTT group_type // makes sense!!!
	
*now count properly the jobs
	bys id: gen n_job = _n
	replace n_job = n_job - 1 // first row is always before school
	replace n_job = . if firmid == 0 // first row is always before school
	replace n_job = 2 if n_job == 3 & group_type == 2
	replace n_job = 3 if n_job == 4 & group_type == 2
	tab n_job group_type // makes sense!!!
	
*save
	xtset id year
	save "`dir_tmp_aux'/OAXACA_POOLED_`J_J_STRINGA'",replace
}	
*********************************************************************************
*	MERGE FIRM EFFECTS INFORMATION
*********************************************************************************
if 0 == 1 {
*start with the micro-data
	 use id year firmid female age female log_dailywages group_type n_job using "`dir_tmp_aux'/OAXACA_POOLED_`J_J_STRINGA'",replace
	 
*set up information to residualize wage
	preserve
	sum year
	local min_year = r(min)
	local max_year = r(max)
	import delimited "tables/controls_effects_POOLED_SAMPLE_`J_J_STRINGA'.csv", encoding(ISO-8859-1) clear // this includes the year effects and the last two numbers are the coefficients in the age polynomial
	list
	
	gen    obs_count = _n
	sum    v1
	scalar NNN=r(N)
	sum    v1 					if obs_count==NNN
	scalar age3_coeff=r(mean)
	sum    v1 					if obs_count==NNN-1
	scalar age2_coeff=r(mean)
	
	gen     year = `min_year'   if obs_count == 1
	replace year = year[_n-1]+1 if obs_count>1 
	drop if year > `max_year'
	
	rename v1 year_effect
	list
	save "`dir_tmp_aux'/controls_effects_POOLED_SAMPLE_`J_J_STRINGA'",replace
	restore
	
*residualize wages
	merge m:1 year using "`dir_tmp_aux'/controls_effects_POOLED_SAMPLE_`J_J_STRINGA'", nogen
	gen age2 = ((age-40)/40)^2
	gen age3 = ((age-40)/40)^3
	bys female: sum log_dailywages
	replace log_dailywages = log_dailywages - age2*age2_coeff - age3*age3_coeff - year_effect
			 
*bring in the psi firm effects (pooled across gender)
	merge m:1 firmid using "`dir_tmp_aux'/PSI_PVR_POOLED_SAMPLE_`J_J_STRINGA'_ESTIMATES_WITH_CROSS_WALK", gen(merge_psi)
		
*bring in the lambdas (pooled across gender)
	xtset id year
	bys id: gen double firmid_lag = firmid[_n-1]
	merge m:1 firmid_lag using "`dir_tmp_aux'/LAMBDA_PVR_POOLED_SAMPLE_`J_J_STRINGA'_ESTIMATES_WITH_CROSS_WALK", gen(merge_lambda)
	
*the ones that I cannot merge are firms that show up only with a contemporaneous effect
	sum lambda if merge_lambda == 2
	drop if merge_lambda == 2
	
*Normalize the lambdas
	sum lambda      if n_job ==1		
	replace lambda=lambda - r(mean)

*Normalize the psi's
	preserve
	use "`dir_tmp_aux'/PVR_POOLED_SAMPLE_group_1_2_3_ESTIMATES_WITH_CROSS_WALK_AND_FIRM_SIZE", replace //calling this file where I have all the info where I can collect the scalar that I need to run the normalization
	sum psi [aw=log_firm_size] if vinG == 1
	scalar psi_bar = r(mean)
	restore
	replace psi    = psi-psi_bar
	
*save the corresponding file
	save "`dir_tmp_aux'/PVR_AGE_GENDER_PROFILE_MEE",replace	
}	
*********************************************************************************
*	DO THE FIRST GRAPH
*********************************************************************************
if 0 == 1 {
*load the big file
	use "`dir_tmp_aux'/PVR_AGE_GENDER_PROFILE_MEE",replace

*do the graph by group type	
	forval gg_coll = 1(1)4 {
		preserve

		if `gg_coll' == 3{
			keep if group_type == 3
			local title_frame = ""
		}

		if `gg_coll' == 2{
			keep if group_type == 2
			local title_frame = ""
		}

		if `gg_coll' == 1{
			keep if group_type == 1
			local title_frame = ""
		}

		if `gg_coll' == 4{
			keep if id > 0
			local title_frame = ""
			
		}
		
	
*Read to collapse
		collapse log_dailywages lambda psi, by(n_job)
		drop if n_job == 0
		save "`dir_tmp_aux'/simple_graph",replace
	
*Do the graph
		twoway (connected psi n_job,       lcolor(gold)    msize(medlarge)  msymbol(triangle) mfcolor(gold)   legend(label(1 "Current firm assignments -- {&psi} "))) ///
			   (connected lambda n_job ,   lcolor(pink)    msize(medlarge)  msymbol(triangle) mfcolor(pink)   legend(label(2 "Lagged  firm assignments -- {&lambda} "))), ///		 		
			   legend(on) legend(order(1 2)) xtitle("") ytitle("Mean") graphregion(color(white))   bgcolor(white)  legend(ring(2) position(6) rows(1)) xlabel(1 "1st Job" 2 "2nd Job" 3 "3rd Job") title(`title_frame')
		graph export "figures/simple_jobs`J_J_STRINGA'_GROUPTYPE_`gg_coll'.pdf",replace
		graph save  "figures/simple_jobs`J_J_STRINGA'_GROUPTYPE_`gg_coll'",replace

		twoway (scatter psi n_job,       lcolor(gold)   msize(medlarge)  msymbol(triangle) mfcolor(gold)   legend(label(1 "Current firm assignments -- {&psi} "))) ///
			   (scatter lambda n_job ,   lcolor(pink)   msize(medlarge) msymbol(triangle) mfcolor(pink)   legend(label(2 "Lagged  firm assignments -- {&lambda} "))), ///		 		
			   legend(on) legend(order(1 2)) xtitle("") ytitle("Mean") graphregion(color(white))   bgcolor(white)  legend(ring(2) position(6) rows(1)) xlabel(1 "1st Job" 2 "2nd Job" 3 "3rd Job") title(`title_frame')
		graph export "figures/type2simple_jobs`J_J_STRINGA'_GROUPTYPE_`gg_coll'.pdf",replace
		graph save  "figures/type2simple_jobs`J_J_STRINGA'_GROUPTYPE_`gg_coll'",replace
	
	
		restore	
	}
}		
*********************************************************************************
*	DO THE SIMPLE OAXACA GRAPH
*********************************************************************************	
if 0  == 1 {
*load the file
	use "`dir_tmp_aux'/PVR_AGE_GENDER_PROFILE_MEE",replace			
	
*Wages 	- Female vs. Male
	gen 		female_wage = log_dailywages if female 	 == 1 & log_dailywages > 0 
	gen 		male_wage 	= log_dailywages if female   == 0 & log_dailywages > 0 
	
*Psi  	- Female vs. Male	
	gen			female_psi	= psi			if female 	== 1
	gen 		male_psi	= psi			if female   == 0
				 		
*Lambda  - Female vs. Male					 		
	gen			female_lam	= lambda		if female 	== 1
	gen 		male_lam	= lambda		if female   == 0
	
	
*Define the grouping/collapse structure	
	
	forval gg_coll = 1(1)4 {
		preserve

		if `gg_coll' == 3{
			keep if group_type == 3
			local title_frame = ""
		}

		if `gg_coll' == 2{
			keep if group_type == 2
			local title_frame = ""
		}

		if `gg_coll' == 1{
			keep if group_type == 1
			local title_frame = ""
		}

		if `gg_coll' == 4{
			keep if id > 0
			local title_frame = ""
			
		}
	
*Read to collapse
		collapse female_wage male_wage female_psi male_psi female_lam male_lam, by(n_job)
		drop if n_job == 0
		save "`dir_tmp_aux'/oaxaca_sample",replace
		
*Create the gaps
		gen		   wage_gap	 	= 	male_wage-female_wage
		gen		   psi_gap		= 	male_psi-female_psi
		gen 	   lam_gap		= 	male_lam-female_lam
		gen 	   pe_gap		=   wage_gap - psi_gap - lam_gap 		if n_job>1
		replace    pe_gap		=   wage_gap - psi_gap 			  		if n_job==1
	
*Do the graph
		twoway (connected wage_gap n_job,   lcolor(ltblue)  msize(medlarge) msymbol(square)   mfcolor(ltblue) legend(label(1 "Gender Wage Gap (Age adjusted)"))) ///
			   (connected psi_gap n_job,    lcolor(gold) 	msize(medlarge) msymbol(triangle) mfcolor(gold)   legend(label(2 "Gender difference in current firm assignments -- {&psi} "))) ///
			   (connected lam_gap n_job ,   lcolor(pink)    msize(medlarge) msymbol(triangle) mfcolor(pink)   legend(label(3 "Gender difference in lagged  firm assignments -- {&lambda} "))) ///
			   (connected pe_gap n_job ,    lcolor(orange)  msize(medlarge) msymbol(diamond) mfcolor(orange)   legend(label(4 "Gender difference in person effects -- {&alpha} "))), ///		 		
			   legend(on) legend(order(1 2 3 4)) xtitle("") ytitle("Mean") graphregion(color(white))   bgcolor(white)  legend(ring(2) position(6) rows(2)) xlabel(1 "1st Job" 2 "2nd Job" 3 "3rd Job") title(`title_frame')
		graph export "figures/Oaxaca_jobs`J_J_STRINGA'_GROUPTYPE_`gg_coll'.pdf",replace
		graph save  "figures/Oaxaca_jobs`J_J_STRINGA'_GROUPTYPE_`gg_coll'",replace
		
		twoway (scatter wage_gap n_job,   lcolor(ltblue)  msize(medlarge) msymbol(square)   mfcolor(ltblue) legend(label(1 "Gender Wage Gap (Age adjusted)"))) ///
			   (scatter psi_gap n_job,    lcolor(gold) 	  msize(medlarge) msymbol(triangle) mfcolor(gold)   legend(label(2 "Gender difference in current firm assignments -- {&psi} "))) ///
			   (scatter lam_gap n_job ,   lcolor(pink)    msize(medlarge) msymbol(triangle) mfcolor(pink)   legend(label(3 "Gender difference in lagged  firm assignments -- {&lambda} "))) ///
			   (scatter pe_gap n_job ,    lcolor(orange)  msize(medlarge) msymbol(diamond) mfcolor(orange)   legend(label(4 "Gender difference in person effects -- {&alpha} "))), ///		 		
			   legend(on) legend(order(1 2 3 4)) xtitle("") ytitle("Mean") graphregion(color(white))   bgcolor(white)  legend(ring(2) position(6) rows(2)) xlabel(1 "1st Job" 2 "2nd Job" 3 "3rd Job") title(`title_frame')
		graph export "figures/type2_Oaxaca_jobs`J_J_STRINGA'_GROUPTYPE_`gg_coll'.pdf",replace
		graph save  "figures/type2_Oaxaca_jobs`J_J_STRINGA'_GROUPTYPE_`gg_coll'",replace
	
		restore	
	}
}
*********************************************************************************
*	BEAUDRY - Di Nardo
*********************************************************************************	
if 0  == 1 {

*Create the file of firm effects with cross walk	
	global file_type = "POOLED_SAMPLE_group_1_2_3_BEAUDRY_DI_NARDO"
	local file_type = "$file_type"
	
*import the file with the pooled estimates, set it in wide form
	import delimited "tables/RAKIM_LAMBDA_`file_type'.csv", encoding(ISO-8859-1) clear
	rename v1 firmid
	rename v2 lambda
	rename v3 identified_lambda
	save "`dir_tmp_aux'/lambda",replace
	

	import delimited "tables/RAKIM_PSI_`file_type'.csv", encoding(ISO-8859-1) clear
	rename v1 firmid
	rename v2 psi
	rename v3 identified_psi
	rename v4 psi_AKM

	merge 1:1 firmid using "`dir_tmp_aux'/lambda", nogen
	rename firmid firmidnew
	distinct firmidnew
	save "`dir_tmp_aux'/PVR_`file_type'_ESTIMATES",replace	

*merge the cross-walk MATLAB
	import delimited "tables/cross_walk_firmids`file_type'.csv", encoding(ISO-8859-1) clear
	rename v1 firmidnew
	rename v2 firmid
	distinct firmidnew
	merge 1:1 firmidnew using "`dir_tmp_aux'/PVR_`file_type'_ESTIMATES"
	keep if _merge ==3
	keep psi* lambda identified_lambda identified_psi firmid
	rename firmid firmidnew
		
*merge the cross walk STATA
	merge 1:1 firmidnew using "`dir_tmp_aux'/dictionary_BEAUDRY_DI_NARDO"
	keep if _merge == 3
	drop _merge
	
*keep only identified effects
	keep if identified_lambda == 1 | identified_psi == 1
	
*normalize lambdas, save file
	sum lambda if firmid == -1 & year == 1984 
	replace lambda = lambda - r(mean)
	sum lambda if firmid > 0
	local mybar = round(r(mean),0.01)
	save "`dir_tmp_aux'/psi_lambda_effects_BEAUDRY_DI_NARDO",replace
	
*reshape, keep just U and N
	keep if firmid == 0 | firmid == -1
	keep lambda year firmid
	replace  firmid  = 2 if firmid == 0
	replace  firmid  = 1 if firmid == -1
	reshape wide lambda, i(year) j(firmid)
	rename lambda1 lambda_N
	rename lambda2 lambda_U
	 
*merge unemployment rate
	merge 1:1 year using "`dir_tmp_aux'/unemployment_rate_VENETO"
	keep if _merge == 3
	
*now do the steps to merge employment weighted average of lambdas and psi by year
	preserve
	
****save the psi
	use "`dir_tmp_aux'/psi_lambda_effects_BEAUDRY_DI_NARDO",clear
	keep psi firmidnew // firmidnew will be the proper link with the original estimation sample
	save "`dir_tmp_aux'/psi_BEAUDRY_DI_NARDO",replace  // saving a file with just contemporaneous firm effects
	
****save the lambdas	
	use "`dir_tmp_aux'/psi_lambda_effects_BEAUDRY_DI_NARDO",clear
	keep lambda firmidnew // firmidnew will be the proper link with the original estimation sample
	rename firmidnew firmid_lag
	save "`dir_tmp_aux'/lambda_BEAUDRY_DI_NARDO",replace  // saving a file with just lagged firm effects
	
****now merge the lambdas (note that they are already normalized)
	use "`dir_tmp_aux'/RAKIM_start",replace
	xtset id year
	bys id: gen firmid_lag = firmidnew[_n-1] 
	merge m:1 firmid_lag using "`dir_tmp_aux'/lambda_BEAUDRY_DI_NARDO", gen(merge_lambda)
	
****collapse
	drop if  firmid == 0 | firmid == -1 // only proper firm lagged effects
	collapse lambda, by(year)
	save "`dir_tmp_aux'/lambda_employers_BEAUDRY_DI_NARDO",replace			
	restore

*merge back
	merge 1:1 year using "`dir_tmp_aux'/lambda_employers_BEAUDRY_DI_NARDO", nogen
	describe
		
*do the line graph
	twoway (connected lambda_N year, yaxis(1) msymbol(square) lcolor(ltblue) mcolor(ltblue) legend(label(1 "{&lambda}{subscript:N}"))) (connected lambda_U year, yaxis(1) lcolor(gold) msymbol(square) mcolor(gold) legend(label(2 "{&lambda}{subscript:U}"))) (connected U_RATE year, yaxis(2) lcolor(orange) msymbol(triangle) mcolor(orange) legend(label(3 "Unemployment Rate in Veneto"))), ytitle("Firm Effects", axis(1)) ytitle("Unemployment Rate", axis(2)) xtitle("Year") legend(order(1 2 3)) legend(ring(2) position(6) rows(1))  plotregion(margin(zero))	xlabel(1984(2)2001)
	graph export "figures/BEAUDRY_DI_NARDO.pdf",replace
	graph save  "figures/BEAUDRY_DI_NARDO",replace
		
*do the scatter plot
	twoway (scatter lambda_N U_RATE, yaxis(1) msymbol(square) mcolor(ltblue) legend(label(1 "{&lambda}{subscript:N}"))) (scatter lambda_U U_RATE, yaxis(1) msymbol(square) mcolor(gold) legend(label(2 "{&lambda}{subscript:U}"))) (lfit lambda_N U_RATE, yaxis(1) lcolor(ltblue) lpattern(solid) legend(label(3 "{&lambda}{subscript:N}"))) (lfit lambda_U U_RATE, yaxis(1) lpattern(solid) lcolor(gold) legend(label(4 "{&lambda}{subscript:U}"))), ytitle("Firm Effects", axis(1)) xtitle("Unemployment Rate") legend(order(1 2)) legend(ring(2) position(6) rows(1)) note("Note {&lambda} normalized w.r.t. {&lambda}{subscript:N} in 1984" "Red Line correspond to average {&lambda} excluding N and U lagged employers")	yline(`mybar', lcolor(red))
	graph export "figures/BEAUDRY_DI_NARDO_scatter.pdf",replace
	graph save  "figures/BEAUDRY_DI_NARDO_scatter",replace
	
*do the scatter plot with lambda employers
	twoway (scatter lambda_N U_RATE, yaxis(1) msymbol(square) mcolor(ltblue) legend(label(1 "{&lambda}{subscript:N}"))) (scatter lambda_U U_RATE, yaxis(1) msymbol(square) mcolor(gold) legend(label(2 "{&lambda}{subscript:U}"))) (lfit lambda_N U_RATE, yaxis(1) lcolor(ltblue) lpattern(solid) legend(label(3 "{&lambda}{subscript:N}"))) (lfit lambda_U U_RATE, yaxis(1) lpattern(solid) lcolor(gold) legend(label(4 "{&lambda}{subscript:U} "))) (scatter lambda U_RATE, yaxis(1) mcolor(orange) msymbol(triangle) legend(label(5 "E[{&lambda}{subscript:l(j,t)} | year = t, l(j,t)=(1, ..., J)]"))) (lfit lambda U_RATE, yaxis(1) lpattern(solid) lcolor(orange) legend(label(6 "{&lambda}{subscript:J}"))), ytitle("Firm Effects", axis(1)) xtitle("Unemployment Rate") legend(order(1 2 5)) legend(ring(2) position(6) rows(2)) note("Note: {&lambda} normalized w.r.t. {&lambda}{subscript:N} in 1984")
	graph export "figures/BEAUDRY_DI_NARDO_scatter_aug.pdf",replace
	graph save  "figures/BEAUDRY_DI_NARDO_scatter_aug",replace
	


}
capture log close