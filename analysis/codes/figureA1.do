local ss				= $it
local dir_tmp_aux 		= "../build/src/RAKIM"	
**************************************************************************

*	SECTION 1: Call the Sigma^2_i
	
**************************************************************************	
if 1 == 1 {
*call the sigma_i first
	import delimited "`dir_tmp_aux'/sigmaIPOOLED_SAMPLE__INVIND_data_1990_2016_2005_2016_RS_build.csv", encoding(ISO-8859-1) clear
	rename v1 sigma_i
	gen rows = _n
	save "`dir_tmp_aux'/sigmas",replace
	
*to do the merge
	import delimited "`dir_tmp_aux'/POOLED_SAMPLE__INVIND_data_1990_2016_2005_2016_RS_build.csv", encoding(ISO-8859-1) clear
	rename v3 firmid 
	rename v4 lagfirmid
	rename v6 fe
	rename v7 lag_fe
	rename v8 controls
	rename v19 age2 
	rename v20 age3
	gen year = 2004
	local aux = 1
	forval pp = 9(1)18 {
	replace year = 2004 + `aux' if v`pp' == 1
	local aux = `aux' +1
	}
	replace year = 2015 if year == 2004
	sum year

*merge sigma_i
	gen rows= _n
	merge 1:1 rows using "`dir_tmp_aux'/sigmas", nogen
	collapse sigma_i fe lag_fe controls, by(firmid lagfirmid age2 age3 year)
	save "`dir_tmp_aux'/identified_firms",replace

*merge in effects	
	use "`dir_tmp_aux'/DAVE_BUILD_RAKIM_APPENDED100_INVIND_data_1990_2016_2005_2016_RS_build",replace	
	gen age2= ((age-40)/40)^2
	gen age3= ((age-40)/40)^3
	merge m:1 firmid lagfirmid age2 age3 year using "`dir_tmp_aux'/identified_firms", gen(identified)
	keep if identified ==3
	
*merge in VA/L
	merge m:1 firmid using "`dir_tmp_aux'/va_file", gen(VA) keepusing(log_VA_L)
	keep if VA == 3	
	
*compute quantiles
	drop qqq
	xtile qqq = log_VA_L if fe!=., nquantile(20)	
	
*compute quantiles
	replace DIP10=log(DIP10)
	xtile qqq_size = DIP10 if fe!=., nquantile(20)		
	
*save
	save "`dir_tmp_aux'/hetero",replace
	 	
}
**************************************************************************

*	SECTION 2: NOW PLOT
	
**************************************************************************		
*load
	use sigma_i log_VA_L qqq using "`dir_tmp_aux'/hetero",replace
	
*collapse 
	collapse sigma_i log_VA_L, by(qqq)

*plotsigma_i
	twoway (scatter sigma_i log_VA_L,    lcolor(dknavy)  msymbol(square)   mfcolor(dknavy) legend(label(1 "Log Hiring Wage"))), /// 		
		legend(off) xtitle("Log Value Added per Worker") ytitle("Average Error Variance") graphregion(color(white))   bgcolor(white)  legend(ring(0) position(5) rows(3))  legend(order(1 2 4)) xlabel(1.50(1)5.5)
		graph export "figures/hetero_VA.pdf",replace
		
*load
	use sigma_i DIP10 qqq_size using "`dir_tmp_aux'/hetero",replace
	
*collapse 
	collapse sigma_i DIP10, by(qqq_size)


*plotsigma_i
	twoway (scatter sigma_i DIP10,    lcolor(dknavy)  msymbol(square)   mfcolor(dknavy) legend(label(1 "Log Hiring Wage"))), /// 		
		legend(off) xtitle("Log Firm Size") ytitle("Average Error Variance") graphregion(color(white))   bgcolor(white)  legend(ring(0) position(5) rows(3))  legend(order(1 2 4)) 
		graph export "figures/hetero_SIZE.pdf",replace
		
*load
	use sigma_i age female using "`dir_tmp_aux'/hetero",replace
	
*collapse 
	collapse sigma_i, by(age female)

*plotsigma_i
	twoway (scatter sigma_i age if female == 1,    lcolor(pink)  msymbol(square)   mfcolor(pink) legend(label(1 "Women"))) /// 		
		(scatter sigma_i age if female == 0,    lcolor(dknavy)  msymbol(square)   mfcolor(dknavy) legend(label(2 "Men"))), ///
		legend(on) xtitle("Age at Hiring") ytitle("Average Error Variance") graphregion(color(white))   bgcolor(white)  legend(ring(0) position(5) rows(3))  legend(order(1 2))
		graph export "figures/Age_female.pdf",replace				
	

}	
