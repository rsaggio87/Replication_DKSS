
*read the global
	local dir_tmp_aux = "/scratch/public/leave_out"	
	local file_type = "$file_type"
	local type_DEC	= 2

if 1 == 1{
*take the psi for male, keep only identified
	use "`dir_tmp_aux'/PSI_PVR_MALE_SAMPLE_group_1_2_3_ESTIMATES_WITH_CROSS_WALK",replace
	keep if identified_psi_male == 1
	merge 1:1 firmid using "`dir_tmp_aux'/firm_size_full_file",gen(merge_firm_size) 
	keep if merge_firm_size==3
	keep psi_male firmid log_firm_size
	save "`dir_tmp_aux'/PSI_PVR_MALE_SAMPLE_group_1_2_3_ESTIMATES_WITH_CROSS_WALK_NORMALIZED_IDENTIFIED",replace 
	
*take the psi for females, keep only identified	
	use "`dir_tmp_aux'/PSI_PVR_FEMALE_SAMPLE_group_1_2_3_ESTIMATES_WITH_CROSS_WALK",replace
	keep if identified_psi_female == 1
	merge 1:1 firmid using "`dir_tmp_aux'/firm_size_full_file",gen(merge_firm_size) 
	keep if merge_firm_size==3
	keep psi_female firmid log_firm_size
	save "`dir_tmp_aux'/PSI_PVR_FEMALE_SAMPLE_group_1_2_3_ESTIMATES_WITH_CROSS_WALK_NORMALIZED_IDENTIFIED",replace
	
*let's drop firms with firm size =1 (they are in the sample becayse we did not impose Pii<1 restrictions)
*keep if log_firm_size>0

*now merge psi(male), psi(females), keep only dual connected firms
	merge 1:1 firmid using "`dir_tmp_aux'/PSI_PVR_MALE_SAMPLE_group_1_2_3_ESTIMATES_WITH_CROSS_WALK_NORMALIZED_IDENTIFIED", gen(dual_merge)
	keep if dual_merge==3
	drop dual_merge
	
*impose the normalization now (lowest vingtile in largest connected set)
	xtile vinG = log_firm_size, nquantile(20)
	gen 	firm_size = exp(log_firm_size)
	sum 	psi_female [aw=firm_size] if vinG == 1
	scalar  psi_bar = r(mean)
	replace psi_female= psi_female - psi_bar
	sum 	psi_male  [aw=firm_size] if vinG == 1
	scalar  psi_bar = r(mean)
	replace psi_male= psi_male - psi_bar
	save "`dir_tmp_aux'/DUAL_CONNECTED_PSI_group_1_2_3_ESTIMATES_WITH_CROSS_WALK_NORMALIZED_IDENTIFIED",replace 
	
*let's now repeat for lambdas. Start with lambda for male
	use "`dir_tmp_aux'/LAMBDA_PVR_MALE_SAMPLE_group_1_2_3_ESTIMATES_WITH_CROSS_WALK",replace
	keep if identified_lambda == 1
	sum lambda if firmid == -1
	replace lambda = lambda - r(mean)
	keep lambda firmid_lag
	describe
	save "`dir_tmp_aux'/LAMBDA_PVR_MALE_SAMPLE_group_1_2_3_ESTIMATES_WITH_CROSS_WALK_NORMALIZED_IDENTIFIED",replace 
			
*Now lambda for females
	use "`dir_tmp_aux'/LAMBDA_PVR_FEMALE_SAMPLE_group_1_2_3_ESTIMATES_WITH_CROSS_WALK",replace
	keep if identified_lambda == 1
	sum lambda if firmid_lag == -1
	replace lambda = lambda - r(mean)
	keep lambda firmid_lag
	describe
	save "`dir_tmp_aux'/LAMBDA_PVR_FEMALE_SAMPLE_group_1_2_3_ESTIMATES_WITH_CROSS_WALK_NORMALIZED_IDENTIFIED",replace 	
	
*now merge lambda(male), lambda(females), keep only dual connected firms
	merge 1:1 firmid_lag using "`dir_tmp_aux'/LAMBDA_PVR_MALE_SAMPLE_group_1_2_3_ESTIMATES_WITH_CROSS_WALK_NORMALIZED_IDENTIFIED", gen(dual_merge)
	keep if dual_merge==3
	drop dual_merge
	save "`dir_tmp_aux'/DUAL_CONNECTED_LAMBDA_group_1_2_3_ESTIMATES_WITH_CROSS_WALK_NORMALIZED_IDENTIFIED",replace 
	
*now let's bring this info back into the auxiliary t=4 file
 	use id year firmid female age female log_dailywages group_type n_job using "`dir_tmp_aux'/OAXACA_POOLED",replace
	
*do the auxiliary steps to be able to age adjust the wages (age and year effects are gender specific now)
	preserve
	sum year
	local min_year = r(min)
	local max_year = r(max)
	import delimited "tables/controls_effects_MALE_SAMPLE_group_1_2_3.csv", encoding(ISO-8859-1) clear // this includes the year effects and the last two numbers are the coefficients in the age polynomial
	list 
	
	gen    obs_count = _n
	sum    v1
	scalar NNN=r(N)
	sum    v1 					if obs_count==NNN
	scalar age3_coeff_male=r(mean)
	sum    v1 					if obs_count==NNN-1
	scalar age2_coeff_male=r(mean)
	
	gen     year = `min_year'   if obs_count == 1
	replace year = year[_n-1]+1 if obs_count>1 
	drop if year > `max_year'
	
	rename v1 year_effect
	gen	   female = 0 
	list
	save "`dir_tmp_aux'/controls_effects_male_only",replace
	
	sum year
	local min_year = r(min)
	local max_year = r(max)
	import delimited "tables/controls_effects_FEMALE_SAMPLE_group_1_2_3.csv", encoding(ISO-8859-1) clear // this includes the year effects and the last two numbers are the coefficients in the age polynomial
	list 
	
	gen    obs_count = _n
	sum    v1
	scalar NNN=r(N)
	sum    v1 					if obs_count==NNN
	scalar age3_coeff_female=r(mean)
	sum    v1 					if obs_count==NNN-1
	scalar age2_coeff_female=r(mean)
	
	gen     year = `min_year'   if obs_count == 1
	replace year = year[_n-1]+1 if obs_count>1 
	drop if year > `max_year'
	
	rename v1 year_effect
	gen	   female = 1 
	list
	append using "`dir_tmp_aux'/controls_effects_male_only"
	list
	save "`dir_tmp_aux'/controls_effects_OAXACA",replace
	restore
	
*ready to perform the adjustment
	merge m:1 year female using "`dir_tmp_aux'/controls_effects_OAXACA", nogen
	gen age2 = ((age-40)/40)^2
	gen age3 = ((age-40)/40)^3
	bys female: sum log_dailywages
	replace log_dailywages = log_dailywages - age2*age3_coeff_male - age3*age2_coeff_male - year_effect     if female == 0
	replace log_dailywages = log_dailywages - age2*age3_coeff_female - age3*age2_coeff_female - year_effect if female == 1

*bring in the psi firm effects
	merge m:1 firmid using "`dir_tmp_aux'/DUAL_CONNECTED_PSI_group_1_2_3_ESTIMATES_WITH_CROSS_WALK_NORMALIZED_IDENTIFIED", gen(merge_psi)
		
*bring in the lambdas
	xtset id year
	bys id: gen double firmid_lag = firmid[_n-1]
	merge m:1 firmid_lag using "`dir_tmp_aux'/DUAL_CONNECTED_LAMBDA_group_1_2_3_ESTIMATES_WITH_CROSS_WALK_NORMALIZED_IDENTIFIED", gen(merge_lambda)	

*keep transitions for which we can identify both lambda and psi
	keep if merge_psi == 3 & merge_lambda == 3
	distinct id
	distinct firmid

*save the corresponding file
	save "`dir_tmp_aux'/PVR_AGE_GENDER_PROFILE_MEE_FULL_OAXACA",replace	
}	
	
*use corresponding file
	use "`dir_tmp_aux'/PVR_AGE_GENDER_PROFILE_MEE_FULL_OAXACA",replace
		
*now create the counterfactuals (DECOMP 1)
if `type_DEC' == 1 {
	sum 		lambda* 										if n_job 	== 1 // must be zero
	gen 		female_wage = log_dailywages 					if female 	== 1 
	gen 		male_wage 	= log_dailywages 					if female   == 0
	
	gen			female_psi	= psi_female						if female 	== 1
	gen 		male_psi	= psi_male							if female   == 0
	gen 		psi_male_fem= psi_female						if female 	== 0	
				 					 		
	gen			female_lam	= lambda_female						if female 	== 1
	gen 		male_lam	= lambda_male						if female   == 0
	gen 		lam_male_fem= lambda_female						if female 	== 0
	
	gen		    barga_psi	= psi_male-psi_female				if female   == 0
	gen			barga_lam	= lambda_male-lambda_female			if female 	== 0
}

if `type_DEC' == 2 {
	sum 		lambda* 										if n_job 	== 1 // must be zero
	gen 		female_wage = log_dailywages 					if female 	== 1 
	gen 		male_wage 	= log_dailywages 					if female   == 0
	
	gen			female_psi	= psi_female						if female 	== 1
	gen 		male_psi	= psi_male							if female   == 0
	gen 		psi_male_fem= psi_male							if female 	== 1	
				 					 		
	gen			female_lam	= lambda_female						if female 	== 1
	gen 		male_lam	= lambda_male						if female   == 0
	gen 		lam_male_fem= lambda_male						if female 	== 1
	
	gen		    barga_psi	= psi_male-psi_female				if female   == 1
	gen			barga_lam	= lambda_male-lambda_female			if female 	== 1
}
	
	sum   		female_psi male_psi psi_male_fem barga_psi

*CCK Time-out
	preserve
	reg 		psi_female psi_male
	scalar      beta_male = _b[psi_male]
	xtile 	    qqq=psi_male,nquantiles(100)
	collapse    psi_female psi_male, by(qqq)
	
	twoway (scatter psi_female psi_male,     msize(medlarge) msymbol(square)   mfcolor(ltblue)   legend(label(1 "scatter"))) (lfit psi_female psi_male,   lcolor(ltblue) legend(label(2 "Lfit"))), ///
	legend(off) xtitle("Psi Male") ytitle("Psi Females") graphregion(color(white))   bgcolor(white)  
	graph export "figures/scatter_CCK`file_type'_GROUPTYPE_`gg_coll'.pdf",replace
	graph save  "figures/scatter_CCK`file_type'_GROUPTYPE_`gg_coll'",replace
	restore
	
*Ready to collapse
		forval gg_coll = 1(1)4 {
		preserve

		if `gg_coll' == 3{
			keep if group_type == 3
			local title_frame = "NJJ"
		}

		if `gg_coll' == 2{
			keep if group_type == 2
			local title_frame = "NJUJJ"
		}

		if `gg_coll' == 1{
			keep if group_type == 1
			local title_frame = "NJJJ"
		}

		if `gg_coll' == 4{
			keep if id > 0
			local title_frame = ""
			
		}
		
		collapse female_wage male_wage female_psi male_psi psi_male_fem female_lam male_lam lam_male_fem barga_psi barga_lam, by(n_job)
		drop if n_job == 0
		save "`dir_tmp_aux'/full_oaxaca_sample",replace
		
*Create the gaps
		gen		   wage_gap	 	= 	male_wage-female_wage
		gen		   psi_gap		= 	male_psi-female_psi
		
		if 	`type_DEC' == 1 {
			gen 	   psi_sort	=   psi_male_fem-female_psi
			gen 	   lam_sort	= 	lam_male_fem-female_lam
		}
		
		if 	`type_DEC' == 2 {
			gen 	   psi_sort	=   male_psi-psi_male_fem
			gen 	   lam_sort	= 	male_lam-lam_male_fem
		}
		
		
		
*Do the graph
if `type_DEC' == 1 {
		twoway (connected wage_gap n_job,    lcolor(ltblue)  msize(medlarge) msymbol(square)   mfcolor(ltblue)   legend(label(1 "Gender Wage Gap (Age adjusted)"))) ///
			   (connected psi_sort n_job,    lcolor(gold) 	 msize(medlarge) msymbol(triangle) mfcolor(gold)     legend(label(2 "Move women to male employers -- {&psi} "))) ///
			   (connected lam_sort n_job ,   lcolor(pink)    msize(medlarge) msymbol(triangle) mfcolor(pink)     legend(label(3 "Move women to male lagged employers -- {&lambda} "))) ///
			   (connected barga_psi n_job ,  lcolor(gold)    msize(medlarge) msymbol(diamond)  mfcolor(gold)     legend(label(4 "Given men female {&psi} "))) ///		
			   (connected barga_lam n_job ,  lcolor(pink)    msize(medlarge) msymbol(diamond)  mfcolor(pink)     legend(label(5 "Given men female {&lambda} "))), ///		 		
			   legend(on) legend(order(2 3 4 5 1)) ylabel(-0.1 0 0.1 0.2 0.3 0.4 0.5) xtitle("") ytitle("Log Daily Wages") graphregion(color(white))   bgcolor(white)  legend(ring(2) position(6) rows(3)) xlabel(1 "1st Job" 2 "2nd Job" 3 "3rd Job") title(`title_frame')
		graph export "figures/DECOM_type_`type_DEC'_FULLOaxaca_jobs`file_type'_GROUPTYPE_`gg_coll'.pdf",replace
		graph save  "figures/DECOM_type_`type_DEC'_FULLOaxaca_jobs`file_type'_GROUPTYPE_`gg_coll'",replace
		
		twoway (connected wage_gap n_job,    lcolor(ltblue)  msize(medlarge) msymbol(square)   mfcolor(ltblue)   legend(label(1 "Gender Wage Gap (Age adjusted)"))) ///
			   (connected psi_sort n_job,    lcolor(gold) 	 msize(medlarge) msymbol(triangle) mfcolor(gold)     legend(label(2 "Move women to male employers -- {&psi} "))) ///
			   (connected lam_sort n_job ,   lcolor(pink)    msize(medlarge) msymbol(triangle) mfcolor(pink)     legend(label(3 "Move women to male lagged employers -- {&lambda} "))), ///
			   legend(on) legend(order(2 3 1)) xtitle("") ylabel(-0.1 0 0.1 0.2 0.3 0.4 0.5) ytitle("Log Daily Wages") graphregion(color(white))   bgcolor(white)  legend(ring(2) position(6) rows(3)) xlabel(1 "1st Job" 2 "2nd Job" 3 "3rd Job") title(`title_frame')
		graph export "figures/DECOM_type_`type_DEC'_SORTING_FULLOaxaca_jobs`file_type'_GROUPTYPE_`gg_coll'.pdf",replace
		graph save  "figures/DECOM_type_`type_DEC'_SORTING_FULLOaxaca_jobs`file_type'_GROUPTYPE_`gg_coll'",replace

	
		restore	
}


if `type_DEC' == 2 {
		twoway (connected wage_gap n_job,    lcolor(ltblue)  msize(medlarge) msymbol(square)   mfcolor(ltblue)   legend(label(1 "Gender Wage Gap (Age adjusted)"))) ///
			   (connected psi_sort n_job,    lcolor(gold) 	 msize(medlarge) msymbol(triangle) mfcolor(gold)     legend(label(2 "Move men to female employers -- {&psi} "))) ///
			   (connected lam_sort n_job ,   lcolor(pink)    msize(medlarge) msymbol(triangle) mfcolor(pink)     legend(label(3 "Move men to female lagged employers -- {&lambda} "))) ///
			   (connected barga_psi n_job ,  lcolor(gold)    msize(medlarge) msymbol(diamond)  mfcolor(gold)     legend(label(4 "Given women male {&psi} "))) ///		
			   (connected barga_lam n_job ,  lcolor(pink)    msize(medlarge) msymbol(diamond)  mfcolor(pink)     legend(label(5 "Given women male {&lambda} "))), ///		 		
			   legend(on) legend(order(2 3 4 5 1)) ylabel(-0.1 0 0.1 0.2 0.3 0.4 0.5) xtitle("") ytitle("Log Daily Wages") graphregion(color(white))   bgcolor(white)  legend(ring(2) position(6) rows(3)) xlabel(1 "1st Job" 2 "2nd Job" 3 "3rd Job") title(`title_frame')
		graph export "figures/DECOM_type_`type_DEC'_FULLOaxaca_jobs`file_type'_GROUPTYPE_`gg_coll'.pdf",replace
		graph save  "figures/DECOM_type_`type_DEC'_FULLOaxaca_jobs`file_type'_GROUPTYPE_`gg_coll'",replace
		
		twoway (connected wage_gap n_job,    lcolor(ltblue)  msize(medlarge) msymbol(square)   mfcolor(ltblue)   legend(label(1 "Gender Wage Gap (Age adjusted)"))) ///
			   (connected psi_sort n_job,    lcolor(gold) 	 msize(medlarge) msymbol(triangle) mfcolor(gold)     legend(label(2 "Move women to male employers -- {&psi} "))) ///
			   (connected lam_sort n_job ,   lcolor(pink)    msize(medlarge) msymbol(triangle) mfcolor(pink)     legend(label(3 "Move women to male lagged employers -- {&lambda} "))), ///
			   legend(on) legend(order(2 3 1)) xtitle("") ylabel(-0.1 0 0.1 0.2 0.3 0.4 0.5) ytitle("Log Daily Wages") graphregion(color(white))   bgcolor(white)  legend(ring(2) position(6) rows(3)) xlabel(1 "1st Job" 2 "2nd Job" 3 "3rd Job") title(`title_frame')
		graph export "figures/DECOM_type_`type_DEC'_SORTING_FULLOaxaca_jobs`file_type'_GROUPTYPE_`gg_coll'.pdf",replace
		graph save  "figures/DECOM_type_`type_DEC'_SORTING_FULLOaxaca_jobs`file_type'_GROUPTYPE_`gg_coll'",replace

	
		restore	
}
}	

	
	