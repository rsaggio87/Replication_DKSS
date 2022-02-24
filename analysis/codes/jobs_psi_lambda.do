
*read the global
	local dir_tmp_aux = "/scratch/public/leave_out"	
	local file_type   = "$file_type"
	
*load the pool datasets with corresponding psi, lambda
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
		
*split by gender
		gen lambda_female	= lambda 	if female == 1
		gen psi_female		= psi 	 	if female == 1
		gen lambda_male		= lambda 	if female == 0
		gen psi_male		= psi	 	if female == 0	
	
*Read to collapse
		collapse log_dailywages lambda* psi*, by(n_job)
		drop if n_job == 0
		save "`dir_tmp_aux'/simple_graph",replace
	
*Do the graph
		twoway (connected psi n_job,       lcolor(gold)    msize(medlarge)  msymbol(triangle) mfcolor(gold)   legend(label(1 "Current firm assignments -- {&psi} "))) ///
			   (connected lambda n_job ,   lcolor(pink)    msize(medlarge)  msymbol(triangle) mfcolor(pink)   legend(label(2 "Lagged  firm assignments -- {&lambda} "))), ///		 		
			   legend(on) legend(order(1 2)) xtitle("") ytitle("Mean") graphregion(color(white))   bgcolor(white)  legend(ring(2) position(6) rows(1)) xlabel(1 "1st Job" 2 "2nd Job" 3 "3rd Job") title(`title_frame')
		graph export "figures/simple_jobs`file_type'_GROUPTYPE_`gg_coll'.pdf",replace
		graph save  "figures/simple_jobs`file_type'_GROUPTYPE_`gg_coll'",replace

		twoway (scatter psi n_job,       lcolor(gold)   msize(medlarge)  msymbol(triangle) mfcolor(gold)   legend(label(1 "Current firm assignments -- {&psi} "))) ///
			   (scatter lambda n_job ,   lcolor(pink)   msize(medlarge) msymbol(triangle) mfcolor(pink)   legend(label(2 "Lagged  firm assignments -- {&lambda} "))), ///		 		
			   legend(on) legend(order(1 2)) xtitle("") ytitle("Mean") graphregion(color(white))   bgcolor(white)  legend(ring(2) position(6) rows(1)) xlabel(1 "1st Job" 2 "2nd Job" 3 "3rd Job") title(`title_frame')
		graph export "figures/type2simple_jobs`file_type'_GROUPTYPE_`gg_coll'.pdf",replace
		graph save  "figures/type2simple_jobs`file_type'_GROUPTYPE_`gg_coll'",replace
		
		twoway (connected psi_male n_job,       lcolor(ltblue)     msize(medlarge)  msymbol(triangle) mfcolor(ltblue)    legend(label(1 "Men Current firm assignments -- {&psi} "))) ///
			   (connected lambda_male n_job ,   lcolor(ltblue)     msize(medlarge)  msymbol(square)   mfcolor(ltblue)   legend(label(2 "Men Lagged  firm assignments -- {&lambda} "))) ///	
			   (connected psi_female n_job,       lcolor(pink) lpattern(solid)     	 msize(medlarge)  msymbol(triangle) mfcolor(pink)    legend(label(3 "Women Current firm assignments -- {&psi} "))) ///
			   (connected lambda_female n_job ,   lcolor(pink)  lpattern(solid)      msize(medlarge)  msymbol(square)   mfcolor(pink)   legend(label(4 "Women Lagged  firm assignments -- {&lambda} "))), ///		 		
			   legend(on) legend(order(1 3 2 4)) xtitle("") ytitle("Mean") graphregion(color(white))   bgcolor(white)  legend(ring(2) position(6) rows(2)) xlabel(1 "1st Job" 2 "2nd Job" 3 "3rd Job") title(`title_frame')
		graph export "figures/BY_GENDER_type2simple_jobs`file_type'_GROUPTYPE_`gg_coll'.pdf",replace
		graph save  "figures/BY_GENDER_type2simple_jobs`file_type'_GROUPTYPE_`gg_coll'",replace
		
	restore		
	}