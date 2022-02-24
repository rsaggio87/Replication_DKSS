
local dir_tmp_aux = "../build/src/RAKIM/"	
local J_J_STRINGA = "1990_2016_2005_2016_RS_build"	
local file_type   = "POOLED_SAMPLE__INVIND_data_`J_J_STRINGA'"
*load
	use "`dir_tmp_aux'/Third_eye",replace

*hiring wage only for matched observations
	rename  log_dailywages log_dailywages_orig
	gen 	log_dailywages=log_dailywages_orig if merge_3==3
	
*wage growth (just for stayers)
	gen 	wage_growth = log_dailywages_3-log_dailywages
	sum wage_growth,d
	
*average at the firm level
	collapse wage_growth, by(firmid)
	
*attache the firm effects
	merge 1:1 firmid using "`dir_tmp_aux'/fe_`file_type'", keepusing(psi lambda DIP10 log_VA_L)
	xtile qqq_aux=log_VA_L [aw=DIP10] if psi!=., nquantiles(100)	
	gen productivity = psi+lambda
	distinct firmid if productivity !=. & log_VA_L!=. & wage_growth!=.
	collapse productivity wage_growth [aw=DIP10], by(qqq)
	merge 1:1 qqq_aux using "`dir_tmp_aux'/separation_bins",nogen
	drop if qqq==.
	save "figures/ucl",replace
}
if 1 == 1{
use "figures/ucl",replace	

*plot
	twoway (scatter wage_growth productivity,    lcolor(dknavy)  msymbol(dknavy)   mfcolor(dknavy)   legend(label(1 "Wage Growth"))) /// 	
		  (scatter separated productivity,    yaxis(2) lcolor(gold)   msymbol(triangle)   mfcolor(gold)   legend(label(2 "Fraction of Workers Separated over First Three Years"))), /// 		
		  legend(on) xtitle("Productivity (Sum of Destination + Origin Effect)") ytitle("Log wage in third year on job - log hiring wage") ytitle("Separation Rate", axis(2)) graphregion(color(white))   bgcolor(white)  legend(ring(0) position(5) rows(3))  ylabel(0(0.05)0.15, axis(1)) ylabel(0(0.2)1, axis(2)) 
		  graph export "figures/ucl.pdf",replace
}
		
