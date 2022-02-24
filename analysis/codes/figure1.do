	local dir_tmp_aux 		= "../build/src/RAKIM/"
	local nQUANTILi   		= 4
	local SAMPLE   	= 0
	local three_jobs_only	= 0
	local baseline			= 0
	local perms_full_tine	= 0
	local adj_wages			= 1
	set seed 1234
	
	if `baseline' == 1 { 
	local SAMPLE = "BASELINE"	
	}
	
	if `baseline' == 0 { 
	local SAMPLE = "ESTIMATION_SAMPLE"	
	}
	
	if `adj_wages' == 1 { 
	local ADJ_WAGES = "RESID_WAGES"	
	}
	
*********************************************************************************

*					SET UP THE DATA TO RUN THE ANALYSIS

*********************************************************************************
if 1 == 1 {

if `baseline' == 1 { // Before pruning
	use id year firmid lagfirmid log_dailywages gap real_age_at_entry female age qualifica tipo_rapporto tipo_contratto reason date  using "`dir_tmp_aux'/DAVE_BUILD_RAKIM_APPENDED100_INVIND_data_1990_2016_2005_2016_RS_build",replace	
	sum year
}

if `baseline' == 0 { // Estimation file
	import delimited "`dir_tmp_aux'/PRUNED_POOLED_SAMPLE__INVIND_data_1990_2016_2005_2016_RS_build.csv", encoding(ISO-8859-1) clear
	
	rename v3 firmid 
	rename v4 lagfirmid
	collapse v1, by(firmid lagfirmid)
	drop v1
	save "`dir_tmp_aux'/identified_p_j_obs",replace
	use id year firmid lagfirmid log_dailywages gap real_age_at_entry female age qualifica tipo_rapporto tipo_contratto reason date  using "`dir_tmp_aux'/DAVE_BUILD_RAKIM_APPENDED100_INVIND_data_1990_2016_2005_2016_RS_build",replace	
	sum year // count # of obs and make sure that is the right starting point
	merge m:1 firmid lagfirmid using "`dir_tmp_aux'/identified_p_j_obs", gen(identified)
	keep if identified == 3
}

	
*Residualize wages
if `adj_wages' == 1 {
	egen cellMIA = group(year real_age_at_entry age female)
	reghdfe log_dailywages, abs(cellMIA) residuals(aux)
	sum aux,d //everytime you re-run this line you get slightly different answers in the residuals!
	replace log_dailywages = aux
}	
	
*now calculate the leave out average (poaching wages only)
	gen uno=1
	bys firmid: egen TOTALE = total(log_dailywages)
	bys firmid: egen DUDES =  total(uno) //at most one year with the firm
	
	gen leave_out = (TOTALE-log_dailywages)/(DUDES-1)
	xtile qqq = leave_out, nquantiles(`nQUANTILi')
	tab qqq
	
*generate the wage difference
	xtset id date
	bys id (date): gen DELTA = log_dailywages - 	log_dailywages[_n-1]
	sum DELTA,d	

*Winsorize
	replace DELTA =r(p99) if DELTA>=r(p99)	& DELTA!=.
	replace DELTA =r(p1) if DELTA<=r(p1)	
	
*now be very careful with the quantiles 
	xtset id date
	bys id: gen contami=_n
	xtset 		id contami
	gen 		qqq_use1      = qqq
	gen 		qqq_use2	  = L1.qqq
	
*save this initial file
	keep id contami lagfirmid firmid date reason gap DELTA qqq_use* uno year
	save "`dir_tmp_aux'/diagnostics",replace
}	
*********************************************************************************

*					BUILD THE TREATMENT

*********************************************************************************	
if 0 == 1{
*load this initial file
	use "`dir_tmp_aux'/diagnostics",replace
		
	
*do you want to keep people with three jobs only
if 	`three_jobs_only' == 1 {
	bys	 id: gen TTT = _N
	keep if 	TTT == 3
	}
*tell me the firmid of your first job, name it consistenly so I can merge the file with DAKM results	
	describe
	rename firmid firmid_true
	xtset id date
	by id (date): gen firmid=firmid_true[_n-1] // tell me the firm in your job #1	
	
*define treatment and control
	gen already_entered =0
	replace already_entered =1 if contami == 1 & lagfirmid>0
	bys id: egen ALREADY = mean(already_entered)
	drop if ALREADY >0 //only people that has initial job obtained from N 	
	keep if contami <= 2 // only the first two jobs
	list id date firmid firmid_true in 1/10

	gen 	  				treatment 	  = 1 		if reason	  != 2 & contami == 2
	replace					treatment 	  = 0 		if reason 	  == 2 & contami == 2 | gap == 1  & reason == . & contami == 2	
		
		
	bys id:		egen TRE		  = mean(treatment)
	
	tab TRE
	drop treatment
	rename TRE treatment		
	
*auxiliary for collapsing
	gen 		DELTA_U		  = DELTA	if 	treatment == 1 
	gen 		DELTA_J		  = DELTA 	if 	treatment == 0
	
	keep if  qqq_use1 !=. & 	qqq_use2 !=.	
	
	sum DELTA_U if qqq_use1 == 4 & qqq_use2 == 1,d
	sum DELTA_J if qqq_use1 == 4 & qqq_use2 == 1,d

	
	gen uno_T = uno if treatment == 1
	gen uno_C = uno if treatment == 0
	

*bring in the lambdas, tell me the average (firm-size weighted) origin effect fot those that were fired.
	preserve
	local J_J_STRINGA = "1990_2016_2005_2016_RS_build"	
	global JJSTRINGA  = "`J_J_STRINGA'"
	global file_type  = "POOLED_SAMPLE__INVIND_data_`J_J_STRINGA'"
	do "codes/firm_effects_PVR"
	restore
	preserve	
	merge m:1 firmid using "`dir_tmp_aux'/fe_POOLED_SAMPLE__INVIND_data_`J_J_STRINGA'", keepusing(lambda DIP10) //tell me the lambda
	keep if treatment == 0
	sum lambda
	collapse DIP10 lambda, by(firmid)
	sum lambda [aw=DIP10]
	restore
	
	
*now collapse	
	collapse DELTA* (sum) uno*, by(qqq_use1 qqq_use2)
	save "figures/dave_dataset`SAMPLE'_`nQUANTILi'`ADJ_WAGES'",replace
	
}

*********************************************************************************

*					DO THE PLOTS

*********************************************************************************		
	use "figures/dave_dataset`SAMPLE'_`nQUANTILi'`ADJ_WAGES'",replace
	
	reg DELTA_U DELTA_J [aw=uno]
	local constante= _b[_cons]
	local CCC = round(`constante',0.01)
	local SSS = round(_b[DELTA_J],0.01)
	local R222= round(e(r2),0.01)	

*plot	
	twoway (scatter DELTA_U DELTA_J, yaxis(1) msymbol(square) mcolor(dknavy)) (lfit DELTA_U DELTA_U, yaxis(1) lcolor(black%30) lpattern(shortdash) lwidth(vthin)), legend(off)  ytitle("Wage change among displaced workers") xtitle("Wage change among poached workers") caption("Constant: `CCC'; Slope: `SSS'; R2: `R222'", pos(11) place(s)) plotregion(margin(zero)) ylabel(-0.6(0.2)0.6) xlabel(-0.6(0.2)0.6) plotregion(margin(zero))
	graph export "figures/reality_check_1`SAMPLE'_`nQUANTILi'`ADJ_WAGES'.pdf",replace
	graph save  "figures/reality_check_1",replace
	
*animate 1	
	twoway (lfit DELTA_U DELTA_U, yaxis(1) lcolor(white) lpattern(solid) lwidth(vthin)), legend(off)  ytitle("Wage Change among Fired") xtitle("Wage Change among Voluntary Quits") plotregion(margin(zero)) ylabel(-0.6(0.2)0.6) xlabel(-0.6(0.2)0.6) plotregion(margin(zero))
	graph export "figures/reality_check_anim_1`SAMPLE'_`nQUANTILi'`ADJ_WAGES'.pdf",replace
	graph save  "figures/reality_anim_1",replace
	
*animate 2	
	twoway (lfit DELTA_U DELTA_U, yaxis(1) lcolor(black%30) lpattern(solid) lwidth(thin)), legend(off)  ytitle("Wage Change among Fired") xtitle("Wage Change among Voluntary Quits") plotregion(margin(zero)) ylabel(-0.6(0.2)0.6) xlabel(-0.6(0.2)0.6) plotregion(margin(zero))
	graph export "figures/reality_check_anim_2`SAMPLE'_`nQUANTILi'`ADJ_WAGES'.pdf",replace
	graph save  "figures/reality_anim_2",replace	
	
*animate 3	
	gen DELTA_shift=DELTA_U-0.05
	twoway (lfit DELTA_U DELTA_U, yaxis(1) lcolor(black%30) lpattern(solid) lwidth(thin)) (lfit DELTA_shift DELTA_U, yaxis(1) lpattern(shortdash) lwidth(thin) lcolor(red)), legend(off)  ytitle("Wage Change among Fired") xtitle("Wage Change among Voluntary Quits") plotregion(margin(zero)) ylabel(-0.6(0.2)0.6) xlabel(-0.6(0.2)0.6) plotregion(margin(zero))
	graph export "figures/reality_check_anim_3`SAMPLE'_`nQUANTILi'`ADJ_WAGES'.pdf",replace
	graph save  "figures/reality_anim_3",replace	
		
	
	
			
