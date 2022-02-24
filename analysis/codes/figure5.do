*read the global
	local file_type 		= "$file_type"
	local dir_tmp_aux 		= "../build/src/RAKIM/"	

*Build firm effects via PVR function
if 1 == 1 {
	do "codes/firm_effects_PVR"
	
*Load the file created by firm_effects_PVR
	use "`dir_tmp_aux'/fe_`file_type'",replace		    
   
*Estimate the slope
	gen log_firm_size = log(DIP10)
	reg psi log_firm_size
	local slope_psi = round(_b[log_firm_size],0.001)
	
	reg lambda log_firm_size
	local slope_lambda = round(_b[log_firm_size],0.001)
	
*visualize the relationship
	xtile qqq=log_firm_size, nquantiles(100)
	
*collapse
	collapse psi lambda log_firm_size, by(qqq)	
			
*plot
	twoway (scatter psi log_firm_size, yaxis(1) msymbol(square) mcolor(ltblue) legend(label(1 "{&psi} -- Destination Firm Effects"))) (lfit psi log_firm_size, yaxis(1) lcolor(ltblue) lpattern(solid) legend(label(2 "Regression Line"))) (scatter lambda log_firm_size, yaxis(1) msymbol(square) mcolor(gold) legend(label(3 "{&lambda} -- Origin Firm Effects"))) (lfit lambda log_firm_size, yaxis(1) lcolor(gold) lpattern(solid) legend(label(4 "Regression Line"))), ytitle("Firm Effects") xtitle("Log Firm Size") note("Regression slope for {&psi}: `slope_psi' " "Regression slope for {&lambda}: `slope_lambda'",  pos(11) place(s)) plotregion(margin(zero)) legend(order(1 3)) legend(ring(2) position(6) rows(1))
	graph export "figures/figure_size`file_type'.pdf",replace
	graph save  "figures/figure_size`file_type'",replace
	
*plot
	twoway (scatter psi log_firm_size, yaxis(1) msymbol(square) mcolor(ltblue) legend(label(1 "{&psi} -- Destination Firm Effects"))) (lfit psi log_firm_size, yaxis(1) lcolor(ltblue) lpattern(solid) legend(label(2 "Regression Line"))) (scatter lambda log_firm_size, yaxis(1) msymbol(square) mcolor(gold) legend(label(3 "{&lambda} -- Origin Firm Effects"))) (lfit lambda log_firm_size, yaxis(1) lcolor(gold) lpattern(solid) legend(label(4 "Regression Line"))), ytitle("Firm Effects") xtitle("Log Firm Size") note("Regression slope for {&psi}: `slope_psi' " "Regression slope for {&lambda}: `slope_lambda'",  pos(11) place(s)) plotregion(margin(zero)) legend(order(1 3)) legend(ring(2) position(6) rows(1))
	graph export "figures/figure_size`file_type'.pdf",replace
	graph save  "figures/figure_size`file_type'",replace
	
	twoway (scatter psi log_firm_size, yaxis(1) msymbol(square) mcolor(ltblue) legend(label(1 "{&psi} -- Destination Firm Effects"))) (lfit psi log_firm_size, yaxis(1) lcolor(ltblue) lpattern(solid) legend(label(2 "Regression Line"))) (scatter lambda log_firm_size, yaxis(1) msymbol(square) mcolor(gold) legend(label(3 "{&lambda} -- Origin Firm Effects"))) (lfit lambda log_firm_size, yaxis(1) lcolor(gold) lpattern(solid) legend(label(4 "Regression Line"))), ytitle("Firm Effects") xtitle("Log Firm Size") plotregion(margin(zero)) legend(order(1 3)) legend(ring(2) position(6) rows(1)) text(0 4.5 "Slope: `slope_lambda'", size(small) color(gold)) text(0.21 4.3 "Slope: `slope_psi'", size(small) color(ltblue))
	graph export "figures/firm_size`file_type'.pdf",replace
	graph save  "figures/firm_size`file_type'",replace
	
	twoway (scatter psi log_firm_size, yaxis(1) msymbol(square) mcolor(ltblue) legend(label(1 "{&psi} -- Destination Firm Effects")))  (scatter lambda log_firm_size, yaxis(1) msymbol(square) mcolor(gold) legend(label(2 "{&lambda} -- Origin Firm Effects"))), ytitle("Firm Effects") xtitle("Log Firm Size") plotregion(margin(zero)) legend(order(1 2)) legend(ring(2) position(6) rows(1))
	graph export "figures/sizEE`file_type'.pdf",replace	
}	
if 1 == 1{	
*Load the file created by firm_effects_PVR
	use "`dir_tmp_aux'/fe_`file_type'",replace	

		
*visualize the relationship
	xtile qqq=log_VA_L [aw=DIP10] if psi!=., nquantiles(100)
	sum log_VA_L [aw=DIP10],d
	local mediana=r(p50)
	
*Fit the MF	
	gen 	log_VA_L_below=0	
	replace log_VA_L_below=log_VA_L if log_VA_L<=`mediana'

	gen 	log_VA_L_above=0	
	replace log_VA_L_above=log_VA_L if log_VA_L>`mediana'
	
	//standard errors obtained after running the lincom procedure for a DWL model, see matlab file "kss_rakim_levels"
	reg psi log_VA_L_below log_VA_L_above [aw=DIP10], r	
	predict fit_psi, xb
	local costante_psi = round(_b[_cons],0.001)
	local slope_psi_below = round(_b[log_VA_L_below],0.001)
	local SE_psi_below=round(0.00060295,0.001) 
	local slope_psi_above = round(_b[log_VA_L_above],0.001)
	local SE_psi_above=round(0.00045875,0.0001) 
	
	reg lambda log_VA_L_below log_VA_L_above [aw=DIP10], r	
	predict fit_lambda, xb
	local costante_lambda = round(_b[_cons],0.0001)
	local slope_lambda_below = round(_b[log_VA_L_below],0.001)
	local SE_lambda_below=round(0.0011682,0.001) 
	local slope_lambda_above = round(_b[log_VA_L_above],0.001)
	local SE_lambda_above=round(0.00088695,0.0001) 

	
*collapse
	collapse psi lambda log_VA_L fit* [aw=DIP10], by(qqq)
*trim
	save "figures/figure_VA`file_type",replace
	
*plot
	twoway (scatter psi log_VA_L, yaxis(1) msymbol(triangle) mcolor(dknavy) legend(label(1 "{&psi} -- Destination Firm Effects")))  (scatter lambda log_VA_L, yaxis(1) msymbol(square) mcolor(cranberry) legend(label(2 "{&lambda} -- Origin Firm Effects"))) (lfit fit_psi log_VA_L, lcolor(dknavy)) (lfit fit_lambda log_VA_L, lcolor(cranberry)), ytitle("Firm Effects") xtitle("Log Value Added per Worker") plotregion(margin(zero)) legend(order(1 2)) legend(ring(2) position(6) rows(1)) text(-0.13 4.1 "Regression slope for {&psi} (below median): `slope_psi_below' (`SE_psi_below') ", size(small) color(dknavy)) text(-0.1 4.1 "Regression slope for {&lambda} (below median): `slope_lambda_below' (`SE_lambda_below')", size(small) color(cranberry)) text(0.5 4.1 "Regression slope for {&psi} (above median): `slope_psi_above' (`SE_psi_above') ", size(small) color(dknavy)) text(0.53 4.1 "Regression slope for {&lambda} (above median): `slope_lambda_above' (`SE_lambda_above')", size(small) color(cranberry))
	graph export "figures/value_added`file_type'_WITH_SE.pdf",replace
}	
*Load the file created by firm_effects_PVR
if 1 == 1 {
	use "`dir_tmp_aux'/fe_`file_type'",replace	
		
*visualize the relationship
	xtile qqq=log_VA_L [aw=DIP10] if psi!=., nquantiles(100)
	
*create the sum of the effects
	gen sum_effects = psi+lambda  		
		
*collapse
	collapse psi lambda sum_effects  [aw=DIP10], by(qqq)
	
*Estimate the slope
	reg psi sum_effects if qqq>=50
	local slope_psi = round(_b[sum_effects],0.001)
	
	reg lambda sum_effects if qqq>=50
	local slope_lambda = round(_b[sum_effects],0.001)	
	
	reg psi sum_effects if qqq<50
	local slope_psi_below = round(_b[sum_effects],0.001)
	
	reg lambda sum_effects if qqq<50
	local slope_lambda_below = round(_b[sum_effects],0.001)
	
*plot
	twoway (scatter psi sum_effects, yaxis(1) msymbol(triangle) mcolor(dknavy) legend(label(1 "{&psi} -- Destination Firm Effects")))  (scatter lambda sum_effects, yaxis(1) msymbol(square) mcolor(cranberry) legend(label(2 "{&lambda} -- Origin Firm Effects"))) (lfit psi sum_effects if qqq<50, yaxis(1) lpattern(dash) lcolor(dknavy)) (lfit psi sum_effects if qqq>=50, yaxis(1) lpattern(dash) lcolor(dknavy)) (lfit lambda sum_effects if qqq<50, yaxis(1) lpattern(dash) lcolor(cranberry)) (lfit lambda sum_effects if qqq>=50, yaxis(1) lpattern(dash) lcolor(cranberry)), ytitle("Firm Effects") xtitle("Sum of Destination and Origin Effects") plotregion(margin(zero)) legend(order(1 2)) legend(ring(2) position(6) rows(1)) text(-0.02 -3.5 "Regression slope for {&psi} (below median): `slope_psi_below' ", size(small) color(dknavy)) text(-0.07 -3.5 "Regression slope for {&lambda} (below median): `slope_lambda_below'", size(small) color(cranberry)) text(0.5 -1.5 "Regression slope for {&psi} (above median): `slope_psi' ", size(small) color(dknavy)) text(0.45 -1.5 "Regression slope for {&lambda} (above median): `slope_lambda'", size(small) color(cranberry))
	graph export "figures/sum_of_effects_`file_type'.pdf",replace	
}
