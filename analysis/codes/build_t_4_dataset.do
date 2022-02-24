
*read the global
	local dir_tmp_aux = "/scratch/public/leave_out"	

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
	save "`dir_tmp_aux'/OAXACA_POOLED",replace