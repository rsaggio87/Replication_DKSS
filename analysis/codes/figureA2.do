*read the global
	local file_type = "$file_type"
	local dir_tmp_aux = "../build/src/RAKIM/"	

*********************************************************************************

*							DAKM VS AKM

*********************************************************************************
*Build firm effects via PVR function
	do "codes/firm_effects_PVR"
	
*Load the file created by firm_effects_PVR
	use "`dir_tmp_aux'/fe_`file_type'",replace	

*evaluate
	reg psi psi_AKM	
	local slopEE = round(_b[psi_AKM],0.001)
	 	
*find the quantiles
	xtile qqq=psi_AKM, nquantiles(100)

*collapse and plot
	collapse psi psi_AKM, by(qqq)

*plot
	twoway (scatter psi psi_AKM, yaxis(1) msymbol(square) mcolor(cranberry) legend(label(1 "{&psi} -- Current Firm Effects"))) (lfit psi psi_AKM, yaxis(1) lcolor(cranberry) lpattern(solid) legend(label(2 "Regression Line"))), legend(off)  ytitle("Firm Effects - DWL") xtitle("Firm Effects - AKM") note("Regression slope: `slopEE'",  pos(11) place(s)) xlabel(-1.5(0.25)1.5) ylabel(-1.5(0.25)1.5) plotregion(margin(zero))
	graph export "figures/figure_AKM_DAKM_`file_type'.pdf",replace
	graph save  "figures/figure_AKM_DAKM_`file_type'",replace
