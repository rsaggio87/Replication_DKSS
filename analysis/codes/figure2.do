	local dir_tmp_aux 	= "../build/src/RAKIM/"
	local nQUANTILi   = 4
	local type_plot	  = 0 
	local baseline    = 0
	set scheme plotplain
*********************************************************************************

*			BUILD DATASET

*********************************************************************************
if 1 == 1 {

if `baseline' == 1 { // Before pruning
	use id year firmid lagfirmid log_dailywages gap real_age_at_entry female age qualifica tipo_rapporto tipo_contratto reason date  using "`dir_tmp_aux'/SUBMITTED_RAKIM_APPENDED_100_INVIND_data_1990_2016_2005_2016",replace	
}


if `baseline' == 0 { // Estimation file
	import delimited "`dir_tmp_aux'/PRUNED_POOLED_SAMPLE__INVIND_data_1990_2016_2005_2016_RS_build.csv", encoding(ISO-8859-1) clear
	rename v3 firmid 
	rename v4 lagfirmid
	collapse v1, by(firmid lagfirmid)
	drop v1
	save "`dir_tmp_aux'/identified_p_j_obs",replace
	use id year firmid lagfirmid log_dailywages gap real_age_at_entry female age qualifica tipo_rapporto tipo_contratto reason date  using "`dir_tmp_aux'/DAVE_BUILD_RAKIM_APPENDED100_INVIND_data_1990_2016_2005_2016_RS_build",replace	
	merge m:1 firmid lagfirmid using "`dir_tmp_aux'/identified_p_j_obs", gen(identified)
	keep if identified == 3
}
	
if `type_plot' == 2 {
	keep if female == 1
	local STRINGA = "FEMALES"
	local TITOLO "WOMEN - ONLY"
	}
	
if `type_plot' == 1 {
	keep if female == 0
	local STRINGA = "MALES"
	local TITOLO "MEN - ONLY"
		}
			
*Residualize wages
	egen cellMIA = group(year real_age_at_entry age female)
	reghdfe log_dailywages, abs(cellMIA) residuals(aux)
	replace log_dailywages = aux
	
*let's keep simple and just focus with people with at least three transitions
	xtset id date 
	bys id: gen TTT = _N
	keep if TTT >=3 
*now count
	xtset id date 
	bys id: gen contami = _n
	xtset id contami

*generate the wage difference
	bys id: gen DELTA = log_dailywages - 	log_dailywages[_n-1]

*Winsorize
	sum DELTA,d
	replace DELTA =r(p99) if DELTA>=r(p99)	& DELTA!=.
	replace DELTA =r(p1) if DELTA<=r(p1)	
	
*now be very careful with the firmids 
	gen 		firmid_graph1 = firmid
	gen 		firmid_graph2 = L1.firmid
	gen 		firmid_graph3 = L2.firmid	

*now calculate the leave out average (poaching wages only)
	gen uno=1
	bys firmid: egen TOTALE = total(log_dailywages)
	bys firmid: egen DUDES =  total(uno) //at most one year with the firm
	
	gen leave_out = (TOTALE-log_dailywages)/(DUDES-1)
	xtile qqq = leave_out, nquantiles(`nQUANTILi')
	
	
*now be very careful with the quantiles 
	xtset 		id contami
	gen 		qqq_use1      = qqq
	gen 		qqq_use2	  = L1.qqq
	gen 		qqq_use3	  = L2.qqq
	
*save this file
	save "`dir_tmp_aux'/diagnostics2",replace	
}	
*********************************************************************************

*			BUILD GROUPS

*********************************************************************************
if 1 == 1 {	
*load this file
	use "`dir_tmp_aux'/diagnostics2",replace			
	
*look for the guys that were unemployed between job#1 and job#2
	gen 		treatment 	  = 0 
	replace 	treatment 	  = 1 if contami == 2 & reason != 2 | contami == 3 & reason != 2
	replace 	treatment	  = 0 if contami == 2 & reason == . & gap == 1 | contami == 3 & reason == . & gap == 1 
	bys	id: egen TRE = total(treatment)
	replace TRE =1 if TRE == 2
	drop treatment
	rename TRE treatment

*now tell me guys that had two separate quits
	gen 		two_quits  	  = 0
	replace		two_quits	  = 1	if contami == 2 & reason == 2 | contami == 3 & reason == 2 
	bys id:		egen Quit	  = total(two_quits)
	replace 	Quit		  = 1 if Quit == 2
	
*now individual quit job#1, fired from job#2
	gen 		mixed	  	  = 0
	replace		mixed		  = 1	if contami == 2 & reason == 2 | contami == 3 & reason !=2  
	replace 	mixed	  	  = 0   if contami == 3 & reason == . & gap == 1	
	bys id:		egen Mixed	  = total(mixed)
	replace 	Mixed		  = 1 if Mixed == 2	
	
*now individual fired job#1, quit from job#2
	gen 		mixed_alt	  	  = 0
	replace		mixed_alt		  = 1	if contami == 3 & reason == 2 | contami == 2 & reason !=2  
	replace 	mixed_alt	  	  = 0   if contami == 2 & reason == . & gap == 1
	bys id:		egen Mixed_alt	  = total(mixed_alt)
	replace 	Mixed_alt		  = 1 if Mixed_alt == 2		
	
	
*auxiliary for collapsing
*	gen 		DELTA1 		  			  = DELTA	if 	qqq_use3 == 1 			& treatment == 1  			|   qqq_use3 == 2 			 	& treatment == 1	
*	gen 		DELTA`nQUANTILi' 		  = DELTA	if 	qqq_use3 == `nQUANTILi' & treatment == 1  			|   qqq_use3 == `nQUANTILi'-1 	& treatment == 1			
*	gen 		DELTA_quits	 			  = DELTA	if  qqq_use3 == `nQUANTILi' & Quit      == 1 	 	    |   qqq_use3 == `nQUANTILi'-1   & Quit      == 1	
	
	xtset id contami
	gen			leave_out_start			  = L2.leave_out
	xtile 		terziale=leave_out_start,nquantiles(3)
	gen 		DELTA1 		  			  = DELTA	if 	terziale == 1 			& treatment == 1  			
	gen 		DELTA`nQUANTILi' 		  = DELTA	if 	terziale == 3 			& treatment == 1  				
	gen 		DELTA_quits	 			  = DELTA	if  terziale == 3 			& Quit      == 1 
	gen 		DELTA_Mixed	 			  = DELTA	if  terziale == 3 			& Mixed     == 1 
	gen 		DELTA_Mixed_alt	 		  = DELTA	if  terziale == 3 			& Mixed_alt == 1 	 	 	

*focus on what you care about
	keep if contami == 3 & qqq_use1!=. & qqq_use2!=. & qqq_use3!=.	

	
*now collapse	
	collapse DELTA* (sum) uno, by(qqq_use1 qqq_use2)
	save "figures/thibaut_dataset",replace
}	

*********************************************************************************

*			PLOT

*********************************************************************************
	use "figures/thibaut_dataset",replace
	
	reg DELTA`nQUANTILi' DELTA1  [aw=uno],r
	local constante= round(_b[_cons],0.001)
	local SSS = round(_b[DELTA1],0.001)
	local R222= round(e(r2),0.001)
	
	reg DELTA_quits DELTA1  [aw=uno],r
	local constante_q= round(_b[_cons],0.001)
	gen cacca=_b[DELTA1]
	sum cacca
	local SSS_q = round(r(mean),0.001) // round gives a problem here fix it manually
	
	local R222_q= round(e(r2),0.001)
	
	reg DELTA_Mixed DELTA1  [aw=uno],r
	local constante_m= round(_b[_cons],0.001)
	local SSS_m = round(_b[DELTA1],0.001)
	local R222_m= round(e(r2),0.001)
	
	reg DELTA_Mixed_alt DELTA1  [aw=uno],r
	local constante_malt= round(_b[_cons],0.001)
	local SSS_malt = round(_b[DELTA1],0.001)
	local R222_malt= round(e(r2),0.001)
	
	hotelling DELTA`nQUANTILi' DELTA1 [aw=uno]
	scalar tsq = r(T2)
	scalar Nhot =r(N)
	scalar dfhot = r(df)
	scalar aiuto = Nhot-dfhot
	scalar Fstat = (dfhot/(Nhot-1))*0.5*tsq
	disp Fstat
	scalar  ppp	 = round(Ftail(aiuto, dfhot, Fstat),0.01)
	disp ppp
	local ppp = ppp
	
*plot
	twoway (scatter DELTA`nQUANTILi' DELTA1, yaxis(1) msymbol(square) mcolor(cranberry)) (lfit DELTA`nQUANTILi' DELTA1, yaxis(1) lcolor(cranberry) lpattern(solid)) (lfit DELTA1 DELTA1, yaxis(1) lcolor(black%30) lpattern(shortdash) lwidth(vthin)), legend(off) ytitle("Wage change among workers first hired by high-wage employer", size(small)) xtitle("Wage change among workers first hired by low-wage employer", size(small)) caption(" `TITOLO' Constant: `constante'; Slope: `SSS'; R2: `R222'", ring(0) place(e)) plotregion(margin(zero))   
	graph export "figures/past_matters`nQUANTILi'`STRINGA'.pdf",replace
	graph save  "figures/past_matters`nQUANTILi'`STRINGA'",replace

*plot Aug
	twoway (scatter DELTA_quits DELTA1, yaxis(1) msymbol(square) mcolor(dknavy) legend(label(1 "Voluntarily Separated from both Job#1 and Job#2 "))) (scatter DELTA`nQUANTILi' DELTA1, yaxis(1) msymbol(square) mcolor(cranberry) legend(label(2 "Involuntarily Separated from both Job#1 and Job#2 "))) (lfit DELTA1 DELTA1, yaxis(1) lcolor(black%90) lpattern(shortdash) lwidth(vthin)),  ytitle("Wage change among workers first hired by high-wage employer", size(small)) xtitle("Wage change among workers first hired by low-wage employer", size(small))  legend(ring(0) position(11) rows(2)) legend(size(vsmall)) plotregion(margin(zero)) legend(order(1 2)) text(-0.3 0.25 " Constant: `constante'; Slope: `SSS'", place(e) size(vsmall) color(cranberry)) text(-0.25 0.25 " Constant: `constante_q'; Slope: .942", place(e) size(vsmall) color(dknavy)) 
*twoway (scatter DELTA_quits DELTA1, yaxis(1) msymbol(square) mcolor(dknavy) legend(label(1 "Quit from both Job#1 and Job#2 "))) (scatter DELTA`nQUANTILi' DELTA1, yaxis(1) msymbol(square) mcolor(cranberry) legend(label(2 "Fired from both Job#1 and Job#2 "))) (lfit DELTA1 DELTA1, yaxis(1) lcolor(black%90) lpattern(shortdash) lwidth(vthin)),  ytitle("Wage Change b/w Job#2 and Job#3 --- Avg. Coworkers Wage in Last Tercile", size(vsmall)) xtitle("Wage Change b/w Job#2 and Job#3 --- Avg. Coworkers Wage Job#1 in First Tercile", size(small))  legend(ring(0) position(11) rows(2)) legend(size(vsmall)) plotregion(margin(zero)) legend(order(1 2)) text(-0.3 0.25 " Constant: `constante'; Slope: `SSS'", place(e) size(vsmall) color(cranberry)) text(-0.25 0.25 " Constant: `constante_q'; Slope: `SSS_q'", place(e) size(vsmall) color(dknavy)) 
*twoway (scatter DELTA_quits DELTA1, yaxis(1) msymbol(square) mcolor(dknavy) legend(label(1 "Quit from both Job#1 and Job#2 "))) (scatter DELTA`nQUANTILi' DELTA1, yaxis(1) msymbol(square) mcolor(cranberry) legend(label(2 "Fired from both Job#1 and Job#2 "))) (scatter DELTA_Mixed DELTA1, yaxis(1) msymbol(square) mcolor(gold) legend(label(3 "Quit from both Job#1, Fired from Job#2 "))) (lfit DELTA1 DELTA1, yaxis(1) lcolor(black%90) lpattern(shortdash) lwidth(vthin)),  ytitle("Wage Change b/w Job#2 and Job#3 --- Avg. Coworkers Wage in Last Tercile", size(vsmall)) xtitle("Wage Change b/w Job#2 and Job#3 --- Avg. Coworkers Wage Job#1 in First Tercile", size(small))  legend(ring(0) position(11) rows(2)) legend(size(vsmall)) plotregion(margin(zero)) legend(order(1 2 3)) text(-0.3 0.25 " Constant: `constante'; Slope: `SSS'", place(e) size(vsmall) color(cranberry)) text(-0.25 0.25 " Constant: `constante_q'; Slope: `SSS_q'", place(e) size(vsmall) color(dknavy)) text(-0.15 0.25 " Constant: `constante_m'; Slope: `SSS_m'", place(e) size(vsmall) color(gold)) 
*twoway (scatter DELTA_quits DELTA1, yaxis(1) msymbol(square) mcolor(dknavy) legend(label(1 "Quit from both Job#1 and Job#2 "))) (scatter DELTA`nQUANTILi' DELTA1, yaxis(1) msymbol(square) mcolor(cranberry) legend(label(2 "Fired from both Job#1 and Job#2 "))) (scatter DELTA_Mixed DELTA1, yaxis(1) msymbol(square) mcolor(gold) legend(label(3 "Quit from Job#1, Fired from Job#2 "))) (scatter DELTA_Mixed_alt DELTA1, yaxis(1) msymbol(square) mcolor(lime) legend(label(4 "Fired from Job#1, Quit from Job#2 "))) (lfit DELTA1 DELTA1, yaxis(1) lcolor(black%90) lpattern(shortdash) lwidth(vthin)),  ytitle("Wage Change b/w Job#2 and Job#3 --- Avg. Coworkers Wage in Last Tercile", size(vsmall)) xtitle("Wage Change b/w Job#2 and Job#3 --- Avg. Coworkers Wage Job#1 in First Tercile", size(small))  legend(ring(0) position(11) rows(2)) legend(size(vsmall)) plotregion(margin(zero)) legend(order(1 2 3 4)) text(-0.3 0.25 " Constant: `constante'; Slope: `SSS'", place(e) size(vsmall) color(cranberry)) text(-0.25 0.25 " Constant: `constante_q'; Slope: `SSS_q'", place(e) size(vsmall) color(dknavy)) text(-0.20 0.25 " Constant: `constante_malt'; Slope: `SSS_malt'", place(e) size(vsmall) color(lime))  text(-0.15 0.25 " Constant: `constante_m'; Slope: `SSS_m'", place(e) size(vsmall) color(gold)) 
*twoway (scatter DELTA_quits DELTA1, yaxis(1) msymbol(square) mcolor(dknavy) legend(label(1 "Quit from both Job#1 and Job#2 "))) (scatter DELTA`nQUANTILi' DELTA1, yaxis(1) msymbol(square) mcolor(cranberry) legend(label(2 "Fired from both Job#1 and Job#2 "))) (scatter DELTA_Mixed DELTA1, yaxis(1) msymbol(square) mcolor(gold) legend(label(3 "Quit from Job#1, Fired from Job#2 "))) (scatter DELTA_Mixed_alt DELTA1, yaxis(1) msymbol(square) mcolor(lime) legend(label(4 "Fired from Job#1, Quit from Job#2 "))) (lfit DELTA1 DELTA1, yaxis(1) lcolor(black%90) lpattern(shortdash) lwidth(vthin)),  ytitle("Wage Change b/w Job#2 and Job#3 --- Avg. Coworkers Wage in Last Tercile", size(vsmall)) xtitle("Wage Change b/w Job#2 and Job#3 --- Avg. Coworkers Wage Job#1 in First Quartile", size(small))  legend(ring(0) position(11) rows(2)) legend(size(vsmall)) plotregion(margin(zero)) legend(order(1 2 3 4)) text(-0.3 0.25 " Constant: `constante'; Slope: `SSS'", place(e) size(vsmall) color(cranberry)) text(-0.25 0.25 " Constant: `constante_q'; Slope: `SSS_q'", place(e) size(vsmall) color(dknavy)) text(-0.20 0.25 " Constant: `constante_malt'; Slope: `SSS_malt'", place(e) size(vsmall) color(lime))  text(-0.15 0.25 " Constant: `constante_m'; Slope: `SSS_m'", place(e) size(vsmall) color(gold)) 

	graph export "figures/power_past_matters`nQUANTILi'`STRINGA'.pdf",replace
	graph save  "figures/power_past_matters`nQUANTILi'`STRINGA'",replace
	
	twoway(scatter DELTA`nQUANTILi' DELTA1, yaxis(1) msymbol(square) mcolor(cranberry)) (lfit DELTA1 DELTA1, yaxis(1) lcolor(black%90) lpattern(shortdash) lwidth(vthin)),  ytitle("Wage change among workers first hired by high-wage employer", size(small)) xtitle("Wage change among workers first hired by low-wage employer", size(small))  legend(off)  plotregion(margin(zero)) legend(order(1 2)) text(-0.3 0.25 " Constant: `constante'; Slope: `SSS'", place(e) size(vsmall) color(cranberry)) 
*twoway (scatter DELTA_quits DELTA1, yaxis(1) msymbol(square) mcolor(dknavy) legend(label(1 "Quit from both Job#1 and Job#2 "))) (scatter DELTA`nQUANTILi' DELTA1, yaxis(1) msymbol(square) mcolor(cranberry) legend(label(2 "Fired from both Job#1 and Job#2 "))) (lfit DELTA1 DELTA1, yaxis(1) lcolor(black%90) lpattern(shortdash) lwidth(vthin)),  ytitle("Wage Change b/w Job#2 and Job#3 --- Avg. Coworkers Wage in Last Tercile", size(vsmall)) xtitle("Wage Change b/w Job#2 and Job#3 --- Avg. Coworkers Wage Job#1 in First Tercile", size(small))  legend(ring(0) position(11) rows(2)) legend(size(vsmall)) plotregion(margin(zero)) legend(order(1 2)) text(-0.3 0.25 " Constant: `constante'; Slope: `SSS'", place(e) size(vsmall) color(cranberry)) text(-0.25 0.25 " Constant: `constante_q'; Slope: `SSS_q'", place(e) size(vsmall) color(dknavy)) 
*twoway (scatter DELTA_quits DELTA1, yaxis(1) msymbol(square) mcolor(dknavy) legend(label(1 "Quit from both Job#1 and Job#2 "))) (scatter DELTA`nQUANTILi' DELTA1, yaxis(1) msymbol(square) mcolor(cranberry) legend(label(2 "Fired from both Job#1 and Job#2 "))) (scatter DELTA_Mixed DELTA1, yaxis(1) msymbol(square) mcolor(gold) legend(label(3 "Quit from both Job#1, Fired from Job#2 "))) (lfit DELTA1 DELTA1, yaxis(1) lcolor(black%90) lpattern(shortdash) lwidth(vthin)),  ytitle("Wage Change b/w Job#2 and Job#3 --- Avg. Coworkers Wage in Last Tercile", size(vsmall)) xtitle("Wage Change b/w Job#2 and Job#3 --- Avg. Coworkers Wage Job#1 in First Tercile", size(small))  legend(ring(0) position(11) rows(2)) legend(size(vsmall)) plotregion(margin(zero)) legend(order(1 2 3)) text(-0.3 0.25 " Constant: `constante'; Slope: `SSS'", place(e) size(vsmall) color(cranberry)) text(-0.25 0.25 " Constant: `constante_q'; Slope: `SSS_q'", place(e) size(vsmall) color(dknavy)) text(-0.15 0.25 " Constant: `constante_m'; Slope: `SSS_m'", place(e) size(vsmall) color(gold)) 
*twoway (scatter DELTA_quits DELTA1, yaxis(1) msymbol(square) mcolor(dknavy) legend(label(1 "Quit from both Job#1 and Job#2 "))) (scatter DELTA`nQUANTILi' DELTA1, yaxis(1) msymbol(square) mcolor(cranberry) legend(label(2 "Fired from both Job#1 and Job#2 "))) (scatter DELTA_Mixed DELTA1, yaxis(1) msymbol(square) mcolor(gold) legend(label(3 "Quit from Job#1, Fired from Job#2 "))) (scatter DELTA_Mixed_alt DELTA1, yaxis(1) msymbol(square) mcolor(lime) legend(label(4 "Fired from Job#1, Quit from Job#2 "))) (lfit DELTA1 DELTA1, yaxis(1) lcolor(black%90) lpattern(shortdash) lwidth(vthin)),  ytitle("Wage Change b/w Job#2 and Job#3 --- Avg. Coworkers Wage in Last Tercile", size(vsmall)) xtitle("Wage Change b/w Job#2 and Job#3 --- Avg. Coworkers Wage Job#1 in First Tercile", size(small))  legend(ring(0) position(11) rows(2)) legend(size(vsmall)) plotregion(margin(zero)) legend(order(1 2 3 4)) text(-0.3 0.25 " Constant: `constante'; Slope: `SSS'", place(e) size(vsmall) color(cranberry)) text(-0.25 0.25 " Constant: `constante_q'; Slope: `SSS_q'", place(e) size(vsmall) color(dknavy)) text(-0.20 0.25 " Constant: `constante_malt'; Slope: `SSS_malt'", place(e) size(vsmall) color(lime))  text(-0.15 0.25 " Constant: `constante_m'; Slope: `SSS_m'", place(e) size(vsmall) color(gold)) 
*twoway (scatter DELTA_quits DELTA1, yaxis(1) msymbol(square) mcolor(dknavy) legend(label(1 "Quit from both Job#1 and Job#2 "))) (scatter DELTA`nQUANTILi' DELTA1, yaxis(1) msymbol(square) mcolor(cranberry) legend(label(2 "Fired from both Job#1 and Job#2 "))) (scatter DELTA_Mixed DELTA1, yaxis(1) msymbol(square) mcolor(gold) legend(label(3 "Quit from Job#1, Fired from Job#2 "))) (scatter DELTA_Mixed_alt DELTA1, yaxis(1) msymbol(square) mcolor(lime) legend(label(4 "Fired from Job#1, Quit from Job#2 "))) (lfit DELTA1 DELTA1, yaxis(1) lcolor(black%90) lpattern(shortdash) lwidth(vthin)),  ytitle("Wage Change b/w Job#2 and Job#3 --- Avg. Coworkers Wage in Last Tercile", size(vsmall)) xtitle("Wage Change b/w Job#2 and Job#3 --- Avg. Coworkers Wage Job#1 in First Quartile", size(small))  legend(ring(0) position(11) rows(2)) legend(size(vsmall)) plotregion(margin(zero)) legend(order(1 2 3 4)) text(-0.3 0.25 " Constant: `constante'; Slope: `SSS'", place(e) size(vsmall) color(cranberry)) text(-0.25 0.25 " Constant: `constante_q'; Slope: `SSS_q'", place(e) size(vsmall) color(dknavy)) text(-0.20 0.25 " Constant: `constante_malt'; Slope: `SSS_malt'", place(e) size(vsmall) color(lime))  text(-0.15 0.25 " Constant: `constante_m'; Slope: `SSS_m'", place(e) size(vsmall) color(gold)) 

	graph export "figures/SIMPLE_power_past_matters`nQUANTILi'`STRINGA'.pdf",replace
	graph save  "figures/SIMPLE_power_past_matters`nQUANTILi'`STRINGA'",replace
	


*animate 1
	twoway (lfit DELTA1 DELTA1, yaxis(1) lcolor(white%30) lpattern(shortdash) lwidth(vthin)), legend(off) note("Sample: Workers fired from Job #1 and Job #2.", size(vsmall)) ytitle("Wage Change b/w Job#2 and Job#3 --- Avg. Coworkers Wage Job#1 > Median", size(vsmall)) xtitle("Wage Change b/w Job#2 and Job#3 --- Avg. Coworkers Wage Job#1 < Median", size(small)) plotregion(margin(zero))	
	graph export "figures/past_matters_a1.pdf",replace
	
*animate 2
	twoway (lfit DELTA1 DELTA1, yaxis(1) lcolor(black%30) lpattern(shortdash) lwidth(vthin)), legend(off) note("Sample: Workers fired from Job #1 and Job #2.", size(vsmall)) ytitle("Wage Change b/w Job#2 and Job#3 --- Avg. Coworkers Wage Job#1 > Median", size(vsmall)) xtitle("Wage Change b/w Job#2 and Job#3 --- Avg. Coworkers Wage Job#1 < Median", size(small)) plotregion(margin(zero))	
	graph export "figures/past_matters_a2.pdf",replace
	
*animate 3
	twoway (scatter DELTA`nQUANTILi' DELTA1, yaxis(1) msymbol(square) mcolor(cranberry)) (lfit DELTA`nQUANTILi' DELTA1, yaxis(1) lcolor(cranberry) lpattern(solid)) (lfit DELTA1 DELTA1, yaxis(1) lcolor(black%30) lpattern(shortdash) lwidth(vthin)), legend(off) note("Sample: Workers fired from Job #1 and Job #2.", size(vsmall)) ytitle("Wage Change b/w Job#2 and Job#3 --- Avg. Coworkers Wage Job#1 > Median", size(vsmall)) xtitle("Wage Change b/w Job#2 and Job#3 --- Avg. Coworkers Wage Job#1 < Median", size(small)) caption(" `TITOLO' Constant: `constante'; Slope: `SSS'; R2: `R222'", ring(0) place(e)) plotregion(margin(zero))  
	graph export "figures/past_matters_a3.pdf",replace
	graph save  "figures/past_matters_a3",replace


}	
