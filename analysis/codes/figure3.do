*read the global
	local dir_tmp_aux 	= "../build/src/RAKIM/"	
	local gg = 1
*********************************************************************************

*			Months of Unemployment b/w jobs --- STARTING SAMPLE

*********************************************************************************
*Load the full file
	use year gap id date firmid lagfirmid reason age real_age_at_entry using "`dir_tmp_aux'/DAVE_BUILD_RAKIM_APPENDED100_INVIND_data_1990_2016_2005_2016_RS_build",clear
	sum year
	
*tell me reason for gap = 1 and involuntary separated
	preserve
	keep if gap == 1 & lagfirmid == 0
	tab reason
	restore
	
	preserve
	keep if gap > 1 & lagfirmid == 0
	tab reason
	restore	
	
*Define group
	if `gg' == 1 {
	local group "startin_sample"
	local note_graph "Starting sample (prior to pruning)"
	}  	

	if `gg' == 2 {
	local group "estimation_sample"
	local note_graph "Estimation sample (after pruning)"
	preserve
	import delimited "`dir_tmp_aux'/PRUNED_POOLED_SAMPLE__INVIND_data_1990_2016_2005_2016_RS_build.csv", encoding(ISO-8859-1) clear
	rename v3 firmid 
	rename v4 lagfirmid
	collapse v1, by(firmid lagfirmid)
	drop v1
	save "`dir_tmp_aux'/identified_p_j_obs",replace
	restore
	merge m:1 firmid lagfirmid using "`dir_tmp_aux'/identified_p_j_obs", gen(identified)
	keep if identified == 3
	
	}  	

	bys year: tab reason, missing
	
	
*Define treatment and control
	gen  		treatment = 0 
	replace		treatment = 1 		if lagfirmid>0
	sum 		gap 				if treatment == 0,d
	local 		gap_C	  = round(r(p50), 0.01)
	local 		gap_C25	  = round(r(p25), 0.01)
	local 		gap_C75	  = round(r(p75), 0.01)
	sum 		gap 				if treatment == 1,d
	local 		gap_T	  = r(p50)
	local 		gap_T25	  = r(p25)
	local 		gap_T75	  = r(p75)
	
	tab 		gap if treatment == 1
	
*plot
	replace gap = 60 if gap>=60	
	twoway (histogram gap if treatment==1, start(1) width(1) bcolor(dknavy%20) lcolor(dknavy)) (histogram gap if treatment==0, start(1) width(1) bcolor(cranberry%20) lcolor(cranberry)), xtitle("Months of Non-employment b/w Jobs") legend(order(1 "Quit" 2 "Displaced" )) text(0.4 40 "Quit --- Median: `gap_T'; p25: `gap_T25'; p75: `gap_T75'", size(small) color(dknavy)) text(0.35 40 "Displaced --- Median: `gap_C'; p25: `gap_C25'; p75: `gap_C75'", size(small) color(cranberry)) ylabel(0(0.10)0.8) legend(ring(2) position(6) rows(1))
	graph export "figures/histogram_month_of_unemployment`group'.pdf", replace
	
