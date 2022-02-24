local dir_tmp_aux 	= "../build/src/RAKIM/"
local corrected   = 1 

if `corrected' == 1{
	local usami = "(Corrected)"
	local name_file = "SORKIN_GRAPH_KSS"
	
}
if `corrected' == 0{
	local usami = ""
	local name_file = "SORKIN_GRAPH_PI"
	
}


*********************************************************************************

*	LOAD .CSV WITH VARIANCE COMPONENT FOR EACH SECTOR CREATED BY MATLAB

*********************************************************************************
	import delimited "`dir_tmp_aux'/`name_file'.csv", encoding(ISO-8859-1) clear
	rename v1 var_psi
	rename v2 var_lambda
	rename v3 cov_psi_lambda
	gen bound_beta  = 0.5 + (var_psi-var_lambda)/(2*(var_psi+var_lambda+2*cov_psi_lambda))
	gen bound_rhos  = sqrt(var_psi/(var_psi+var_lambda+2*cov_psi_lambda))*(1-0.3*sqrt(var_lambda/(var_psi+var_lambda+2*cov_psi_lambda)))
	
	gen settore = _n-1
	*replace bound_beta = round(bound_beta,0.01)
	*label define industrie 1 "Retail" 2 "Construction" 3 "Restaurants/Hotels" 4 "Hairdressing/Care Centers" 5 "Law Firms" 6 "Manufacturing" 7 "Transportation" 8 "Cleaning/Security" 9 "Temp Agencies" 10 "Management/Consulting"
	*label values settore industrie


*Merge variance of log daily wages within each sector
if 1 == 1 {
	preserve
	import delimited "`dir_tmp_aux'/POOLED_SAMPLE__INVIND_data_1990_2016_2005_2016_RS_build.csv", encoding(ISO-8859-1) clear
	rename v1 log_dailywages
	rename v3 firmid
	merge m:1 firmid using "`dir_tmp_aux'/FIRM_SIZE",gen(sampling)
	keep if sampling == 3
	collapse (sd) log_dailywages [aw=DIP10], by(settore)
	save "`dir_tmp_aux'/varianze",replace
	restore	
}	
	merge 1:1 settore using "`dir_tmp_aux'/varianze"
	sort settore
	order settore log_dailywages var_psi  var_lambda cov_psi_lambda
	replace  var_psi = var_psi^0.5
	replace  var_lambda = var_lambda^0.5
	replace cov_psi_lambda = cov_psi_lambda/(var_psi*var_lambda)
	*drop if settore == 0
	list
	save "figures/`name_file'",replace
	gen aux = 0
	replace aux = 0.3 if settore == 1
	
*Relabel the sectors
	forval pp=1(1)10{
	sum bound_beta if settore == `pp'
	local meanizzami_`pp'=round(r(mean),0.001)
	
	}
	local meanizzami_1=0.83
	local meanizzami_3=0.94
	local meanizzami_7=0.94
	*label define industrie 1 "Retail ({&beta} >= `meanizzami_1')" 2 "Construction ({&beta} >= `meanizzami_2')" 3 "Restaurants/Hotels ({&beta} >= .94)" 4 "Hairdressing/Care Centers ({&beta} >= `meanizzami_4')" 5 "Law Firms ({&beta} >= `meanizzami_5')" 6 " " 7 "Transportation ({&beta} >= .94)" 8 "Cleaning/Security  ({&beta} >= `meanizzami_8')" 9 "Temp Agencies ({&beta} >= `meanizzami_9')" 10 "Management/Consulting/Tech({&beta} >= `meanizzami_10')"
	label define industrie 0 "Other" 1 "Retail" 2 "Construction" 3 "Restaurants/Hotels" 4 " " 5 "Law Firms" 6 "Manufacturing" 7 " " 8 "Transportation/Cleaning/Security" 9 "Temp Agencies" 10 " " 11 "Banking/Finance" 12 "Education/Health",replace
	label values settore industrie

	
*Plot	
	*twoway (scatter var_psi var_lambda, yaxis(1) msymbol(square) mcolor(dknavy) mlabel(settore)) (lfit aux aux, lcolor(dark) lcolor(black%30) lpattern(shortdash) lwidth(vthin)), legend(off)  ytitle("Std. Destination Effects `usami' ") xtitle("Std. Origin Effects `usami'") text(0.1583 0.0175 "Manufacturing ({&beta} >= .84)", size(small)) ylabel(0(0.10)0.3) xlabel(0(0.10)0.3) plotregion(margin(zero))
	twoway (scatter var_psi var_lambda, yaxis(1) msymbol(square) mcolor(dknavy) mlabel(settore)) (lfit aux aux, lcolor(dark) lcolor(black%30) lpattern(shortdash) lwidth(vthin)), legend(off)  ytitle("Std. Destination Effects `usami' ") xtitle("Std. Origin Effects `usami'") text(0.23 0.0375 "Hairdressing" "Care Centers", size(small) ) ylabel(0(0.10)0.3) text(0.275 0.052 "Managment" "Consulting/Tech", size(small)) ylabel(0(0.10)0.3) xlabel(0(0.10)0.3) plotregion(margin(zero))
	
	graph export "figures/`name_file'.pdf",replace
	graph save  "figures/`name_file'",replace
	save "figures/sorkin_graph",replace

	
