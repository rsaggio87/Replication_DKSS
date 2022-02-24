	local nQUANTILi   		= 100
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
	local SAMPLE_TYPE = "POOLED"
	local J_J_STRINGA = "1990_2016_2005_2016_RS_build"	
	local file_type ="`SAMPLE_TYPE'_SAMPLE__INVIND_data_`J_J_STRINGA'"
	local dir_tmp_aux = "/scratch/public/leave_out"
*********************************************************************************

*					THIS IS TO DO THE 3-D PLOT

*********************************************************************************	
if 0 == 1{
*load this initial file
	use "`dir_tmp_aux'/diagnostics",replace
		
	
*define treatment and control
	gen 	  				treatment 	  = 1 		if reason	  != 2 & contami >= 1
	replace					treatment 	  = 0 		if reason 	  == 2 & contami >= 1 | gap == 1  & reason == . & contami >= 1	
		
		
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
	
*now collapse	
	collapse DELTA* (sum) uno*, by(qqq_use1 qqq_use2)
	save "figures/symmetry_plot",replace
}
	
*********************************************************************************

*					THIS IS TO DO THE NON-PARAMETRIC TEST

*********************************************************************************	
if 1 == 1{
*load this initial file
	use "`dir_tmp_aux'/diagnostics",replace
		
*let's just keep people that were fired from their first job
	gen already_entered =0
	replace already_entered =1 if contami == 1 & lagfirmid>0
	bys id: egen ALREADY = mean(already_entered)
	drop if ALREADY >0 //only people that has initial job obtained from N 	
	keep if contami <= 2 // only the first two jobs
	gen					treatment	  = 0
	replace 	  		treatment 	  = 1 if lagfirmid == 0
	bys id: 			egen TRE		  = mean(treatment)	
	keep if TRE>0
	list id date firmid  in 1/10
	
*reshape the data
	xtset id contami
	gen firmid_d = firmid
	bys id: gen firmid_o=firmid[_n-1]
	 
*keep only stuff that matters
	keep if contami == 2
	keep id firmid_o firmid_d DELTA
	list in 1/10
	
*let's drop crazy changes in Wages
	gen modulo = abs(DELTA)
	drop if modulo >=1	
	
*bring information on firm size
	gen firmid = firmid_o
	merge m:1 firmid using "`dir_tmp_aux'/fe_`file_type'", nogen keepusing(DIP10)
	drop firmid
	rename DIP10 DIP10_o
	gen firmid = firmid_d
	merge m:1 firmid using "`dir_tmp_aux'/fe_`file_type'", nogen keepusing(DIP10)
	drop firmid
	rename DIP10 DIP10_d

*calculate change in firm size
	gen change_size = (DIP10_d-DIP10_o)/DIP10_o
	sum change_size,d
	drop change_size
	sum DIP10_d,d
	replace DIP10_d=log(DIP10_d)
	replace DIP10_o=log(DIP10_o)
	gen change_size = DIP10_d-DIP10_o
	sum change_size,d	
		
*now order things
	gen direction = 0
	replace direction = 1 if change_size >0 
	preserve
	keep if direction == 1
	gen aux_d = firmid_d
	gen aux_o = firmid_o
	drop firmid*
	gen firmid_o = aux_d
	gen firmid_d = aux_o //swapping things
	rename DELTA DELTA_positive
	gen uno = 1
	collapse DELTA_positive (sum) uno, by(firmid_o firmid_d)
	rename uno teste_1
	save "`dir_tmp_aux'/direction",replace
	restore	
	keep if direction== 0
	gen uno = 1
	collapse DELTA change_size (sum) uno, by(firmid_o firmid_d) 
	rename DELTA DELTA_negative
	rename uno teste
	merge 1:1 firmid_o firmid_d using "`dir_tmp_aux'/direction", gen(identified)
	keep if identified ==3
	gen peso=teste+teste_1
	replace change_size = abs(change_size)
	
*save this dataset, this completes the build
	xtile qqq = change_size, nquantiles(`nQUANTILi')
	save "`dir_tmp_aux'/hotelling",replace
}		
*reshape the data for hotelling
	use "`dir_tmp_aux'/hotelling",replace
	gen pairs=1	
	gen differenza = 	DELTA_positive+DELTA_negative
	sum differenza,d
		
	
*binscatter	
	collapse (sum) pairs (mean) differenza DELTA_positive DELTA_negative change_size (semean) se_diff=differenza, by(qqq)
	gen T = differenza/se_diff
	sum T,d
	replace T = T^2
	sum T
	scalar Ttest = r(mean)*r(N)
	disp Ttest
	reg DELTA_negative DELTA_positive 
	local constante= _b[_cons]
	local CCC = round(`constante',0.01)
	local SSS = round(_b[DELTA_positive],0.01)
	gen meno = -DELTA_negative

*plot
	twoway (scatter DELTA_positive DELTA_negative, yaxis(1) msymbol(square) mcolor(cranberry)) (lfit meno DELTA_negative, yaxis(1) lcolor(cranberry) lpattern(solid)), legend(off) ytitle("Wage change (Positive Change in Firm Size)") xtitle("Wage change (Negative Change in Firm Size)") plotregion(margin(zero))   
	graph export "figures/symmetry_nP`nQUANTILi'`STRINGA'.pdf",replace
	graph save  "figures/symmetry_nP`nQUANTILi'`STRINGA'",replace
	
		