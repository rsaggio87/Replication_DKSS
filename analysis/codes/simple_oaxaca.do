
*read the global
	local dir_tmp_aux = "/scratch/public/leave_out"	
	local file_type = "$file_type"

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
		graph export "figures/Oaxaca_jobs`file_type'_GROUPTYPE_`gg_coll'.pdf",replace
		graph save  "figures/Oaxaca_jobs`file_type'_GROUPTYPE_`gg_coll'",replace
		
		twoway (scatter wage_gap n_job,   lcolor(ltblue)  msize(medlarge) msymbol(square)   mfcolor(ltblue) legend(label(1 "Gender Wage Gap (Age adjusted)"))) ///
			   (scatter psi_gap n_job,    lcolor(gold) 	  msize(medlarge) msymbol(triangle) mfcolor(gold)   legend(label(2 "Gender difference in current firm assignments -- {&psi} "))) ///
			   (scatter lam_gap n_job ,   lcolor(pink)    msize(medlarge) msymbol(triangle) mfcolor(pink)   legend(label(3 "Gender difference in lagged  firm assignments -- {&lambda} "))) ///
			   (scatter pe_gap n_job ,    lcolor(orange)  msize(medlarge) msymbol(diamond) mfcolor(orange)   legend(label(4 "Gender difference in person effects -- {&alpha} "))), ///		 		
			   legend(on) legend(order(1 2 3 4)) xtitle("") ytitle("Mean") graphregion(color(white))   bgcolor(white)  legend(ring(2) position(6) rows(2)) xlabel(1 "1st Job" 2 "2nd Job" 3 "3rd Job") title(`title_frame')
		graph export "figures/type2_Oaxaca_jobs`file_type'_GROUPTYPE_`gg_coll'.pdf",replace
		graph save  "figures/type2_Oaxaca_jobs`file_type'_GROUPTYPE_`gg_coll'",replace
	
		restore	
	}