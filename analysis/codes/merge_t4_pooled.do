
*read the global
	local dir_tmp_aux = "/scratch/public/leave_out"	
	local file_type = "$file_type"

*load the auxiliary t=4 file
 	use id year firmid female age female log_dailywages group_type n_job using "`dir_tmp_aux'/OAXACA_POOLED",replace
	 
*set up information to residualize wage
	preserve
	sum year
	local min_year = r(min)
	local max_year = r(max)
	import delimited "tables/controls_effects_`file_type'.csv", encoding(ISO-8859-1) clear // this includes the year effects and the last two numbers are the coefficients in the age polynomial
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
	save "`dir_tmp_aux'/controls_effects_`file_type'",replace
	restore
	
*residualize wages
	merge m:1 year using "`dir_tmp_aux'/controls_effects_`file_type'", nogen
	gen age2 = ((age-40)/40)^2
	gen age3 = ((age-40)/40)^3
	bys female: sum log_dailywages
	replace log_dailywages = log_dailywages - age2*age2_coeff - age3*age3_coeff - year_effect
			 
*bring in the psi firm effects (pooled across gender)
	merge m:1 firmid using "`dir_tmp_aux'/PSI_PVR_`file_type'_ESTIMATES_WITH_CROSS_WALK", gen(merge_psi)
		
*bring in the lambdas (pooled across gender)
	xtset id year
	bys id: gen double firmid_lag = firmid[_n-1]
	merge m:1 firmid_lag using "`dir_tmp_aux'/LAMBDA_PVR_`file_type'_ESTIMATES_WITH_CROSS_WALK", gen(merge_lambda)
	
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