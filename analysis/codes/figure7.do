*read the global
	local file_type 		= "$file_type"
	local dir_tmp_aux 		= "../build/src/RAKIM/"	

*load file from oaxaca build
	use "`dir_tmp_aux'/firm_level_combined_fe_COMBINED_SAMPLE_INVIND_data_1990_2016_2005_2016_RS_build",replace
		    
	
*visualize the relationship
	xtile qqq=log_VA_L [aw=DIP10] if psi_female!=., nquantiles(100)
	replace qqq=-1 if firmid==0
		
*collapse
	collapse psi* lambda* log_VA_L, by(qqq)
	save "figures/oaxaca_VA`file_type",replace
	
	
*slope
	reg psi_female psi_male
	local slope_psi = round(_b[psi_male],0.001)		
	local con_psi = round(_b[_cons],0.001)		
	
	reg lambda_female lambda_male
	local slope_lambda = round(_b[lambda_male],0.001)		
	local con_lambda = round(_b[_cons],0.001)

	
*plot
	twoway (scatter psi_female log_VA_L, msymbol(square) lcolor(pink) mcolor(pink) legend(label(1 "{&psi} -- Destination Firm Effects, Women"))) (scatter psi_male log_VA_L, msymbol(square) lcolor(dknavy) mcolor(dknavy) legend(label(2 "{&psi} -- Destination Firm Effects, Men"))), ytitle("Firm Effects") xtitle("Log Value Added per Worker")  legend(order(1 2)) legend(ring(2) position(6) rows(1))
	graph export "figures/cck_psi`file_type'.pdf",replace
	
	twoway (scatter lambda_female log_VA_L, msymbol(square) lcolor(pink) mcolor(pink) legend(label(1 "{&lambda} -- Origin Firm Effects, Women"))) (scatter lambda_male log_VA_L, msymbol(square) lcolor(dknavy) mcolor(dknavy) legend(label(2 "{&lambda} -- Origin Firm Effects, Men"))), ytitle("Firm Effects") xtitle("Log Value Added per Worker")  legend(order(1 2)) legend(ring(2) position(6) rows(1)) ylabel(0(0.05)0.2)
	graph export "figures/cck_lambda`file_type'.pdf",replace
	
	twoway (scatter psi_female psi_male, msymbol(square) lcolor(cranberry) mcolor(cranberry)) (lfit psi_female psi_male, lcolor(black%30) lpattern(shortdash) lwidth(vthin)), ytitle("Destination Firm Effects, Female (Normalized)") xtitle("Destination Firm Effects, Male  (Normalized)") text(0 0.35 "Constant:`con_psi'." "Slope: `slope_psi'.", size(small) color(cranberry)) legend(off)
	graph export "figures/slope_trimmed_cck_psi`file_type'.pdf",replace
	
	twoway (scatter lambda_female lambda_male if qqq==-1, msymbol(circle) lcolor(orange) mcolor(orange)) (scatter lambda_female lambda_male if qqq>0, msymbol(square) lcolor(dknavy) mcolor(dknavy)) (lfit lambda_female lambda_male if qqq>0, lcolor(black%30) lpattern(shortdash) lwidth(vthin)), ytitle("Origin Firm Effects, Female  (Normalized)") xtitle("Origin Firm Effects, Male (Normalized)") text(0.06 0.11 "Constant:`con_lambda'." "Slope: `slope_lambda'.", size(small) color(dknavy)) legend(off) xlabel(0.015(0.015)0.12) ylabel(0.015(0.015)0.12) text(0.022 0.0165 "{&lambda}{subscript:U}", size(small) color(orange)) 
	graph export "figures/slope_trimmed_cck_lambda`file_type'.pdf",replace

	
*plot trimmed
	drop if qqq<=3
	twoway (scatter psi_female log_VA_L, msymbol(square) lcolor(pink) mcolor(pink) legend(label(1 "{&psi} -- Destination Firm Effects, Women"))) (scatter psi_male log_VA_L, msymbol(square) lcolor(dknavy) mcolor(dknavy) legend(label(2 "{&psi} -- Destination Firm Effects, Men"))), ytitle("Firm Effects") xtitle("Log Value Added per Worker")  legend(order(1 2)) legend(ring(2) position(6) rows(1))
	graph export "figures/trimmed_cck_psi`file_type'.pdf",replace
	
	twoway (scatter lambda_female log_VA_L, msymbol(square) lcolor(pink) mcolor(pink) legend(label(1 "{&lambda} -- Origin Firm Effects, Women"))) (scatter lambda_male log_VA_L, msymbol(square) lcolor(dknavy) mcolor(dknavy) legend(label(2 "{&lambda} -- Origin Firm Effects, Men"))), ytitle("Firm Effects") xtitle("Log Value Added per Worker")  legend(order(1 2)) legend(ring(2) position(6) rows(1)) ylabel(0(0.05)0.2)
	graph export "figures/trimmed_cck_lambda`file_type'.pdf",replace
	

	
