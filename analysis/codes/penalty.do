
*start with pooled
	local SAMPLE_TYPE = "POOLED" //change this to have same calculations for W/M
	local J_J_STRINGA = "1990_2016_2005_2016_RS_build"	
	local file_type ="`SAMPLE_TYPE'_SAMPLE__INVIND_data_`J_J_STRINGA'"
	local dir_tmp_aux 	= "../build/src/RAKIM/"

*import the file with the estimates from the DAKM model
	import delimited "`dir_tmp_aux'/`file_type'.csv", encoding(ISO-8859-1) clear
	rename v2 id
	rename v3 firmid
	rename v4 lagfirmid
	rename v7 lambda
	
*gen the three states
	gen STATE = ""
	replace STATE = "U" if lagfirmid == 0
	replace STATE = "N" if lagfirmid == -1
	replace STATE = "J" if lagfirmid >0
	
*normalize
	sum lambda if STATE == "N"
	replace lambda = lambda-r(mean)
	
*show me the ATT
if 1 == 1 {
	preserve
	drop lambda
	gen rows = _n
	bys id (rows): gen contami = _n
	xtset id contami
	bys id: gen firmid_aux = firmid[_n-1]
	drop firmid
	gen firmid = firmid_aux
	merge m:1 firmid using "`dir_tmp_aux'/fe_`file_type'", nogen keepusing(lambda DIP10)
	keep if STATE == "U"
	sum lambda
	restore	
}

*show me the ATU
if 1 == 1 {
	preserve
	drop lambda
	gen rows = _n
	bys id (rows): gen contami = _n
	xtset id contami
	bys id: gen firmid_aux = firmid[_n-1]
	drop firmid
	gen firmid = firmid_aux
	merge m:1 firmid using "`dir_tmp_aux'/fe_`file_type'", nogen keepusing(lambda DIP10)
	keep if STATE == "J"
	sum lambda
	restore	
}	
	
*show me plot of quit rates by destination/origin effects
	preserve
	gen quit = 0
	replace quit = 1 if STATE == "J"
	drop lambda
	gen rows = _n
	bys id (rows): gen contami = _n
	xtset id contami
	bys id: gen firmid_aux = firmid[_n-1]
	drop firmid
	gen firmid = firmid_aux
	collapse quit,by(firmid)
	merge 1:1 firmid using "`dir_tmp_aux'/fe_`file_type'"
	xtile qqq=log_VA_L [aw=DIP10] if psi!=., nquantiles(100)	
	collapse quit psi log_VA_L lambda [aw=DIP10], by(qqq)
	reg quit log_VA_L
	twoway (scatter quit psi, yaxis(1) msymbol(triangle) mcolor(square) legend(label(1 "{&psi} -- Destination Firm Effects"))), ytitle("Quit Rates") xtitle("Destination Effects") plotregion(margin(zero)) legend(order(1 2)) legend(ring(2) position(6) rows(1))
	graph export "figures/quit_destin`file_type'.pdf",replace
	graph save  "figures/quit_destin`file_type'",replace
	twoway (scatter quit log_VA_L, yaxis(1) msymbol(square) mcolor(ltblue) legend(label(1 "{&psi} -- Destination Firm Effects"))), ytitle("Quit Rates") xtitle("Log Value Added per Worker") plotregion(margin(zero)) legend(order(1 2)) legend(ring(2) position(6) rows(1)) legend(off)
	graph export "figures/quit_VA`file_type'.pdf",replace
	graph save  "figures/quit_VA`file_type'",replace
	restore		
		
*report the shares
	gen 		pi = .
	gen 		pi_U = 0
	replace		pi_U = 1 if STATE == "U"
	sum 		pi_U 
	replace		pi = r(mean) if STATE == "U"
	gen 		pi_J = 0
	replace		pi_J = 1 if STATE == "J"
	sum 		pi_J 
	replace		pi = r(mean) if STATE == "J"	

*collapse
	collapse lambda pi, by(STATE)
	gen TYPE = "`SAMPLE_TYPE'"
	save "`dir_tmp_aux'/penalty",replace

*now men
	local SAMPLE_TYPE = "MALE"
	local J_J_STRINGA = "1990_2016_2005_2016_RS_build"	
	local file_type ="`SAMPLE_TYPE'_SAMPLE__INVIND_data_`J_J_STRINGA'"
	local dir_tmp_aux = "/scratch/public/leave_out"	

*import the file with the estimates from the DAKM model
	import delimited "`dir_tmp_aux'/`file_type'.csv", encoding(ISO-8859-1) clear
	rename v2 id
	rename v3 firmid
	rename v4 lagfirmid
	rename v7 lambda
	
*gen the three states
	gen STATE = ""
	replace STATE = "U" if lagfirmid == 0
	replace STATE = "N" if lagfirmid == -1
	replace STATE = "J" if lagfirmid >0
	
*normalize
	sum lambda if STATE == "N"
	replace lambda = lambda-r(mean)
	
*show me the ATT
if 1 == 1 {
	preserve
	drop lambda
	gen rows = _n
	bys id (rows): gen contami = _n
	xtset id contami
	bys id: gen firmid_aux = firmid[_n-1]
	drop firmid
	gen firmid = firmid_aux
	merge m:1 firmid using "`dir_tmp_aux'/fe_`file_type'", nogen keepusing(lambda DIP10)
	keep if STATE == "U"
	sum lambda
	restore	
}

*show me the ATU
if 1 == 1 {
	preserve
	drop lambda
	gen rows = _n
	bys id (rows): gen contami = _n
	xtset id contami
	bys id: gen firmid_aux = firmid[_n-1]
	drop firmid
	gen firmid = firmid_aux
	merge m:1 firmid using "`dir_tmp_aux'/fe_`file_type'", nogen keepusing(lambda DIP10)
	keep if STATE == "J"
	sum lambda
	restore	
}

*report the shares
	gen 		pi = .
	gen 		pi_U = 0
	replace		pi_U = 1 if STATE == "U"
	sum 		pi_U 
	replace		pi = r(mean) if STATE == "U"
	gen 		pi_J = 0
	replace		pi_J = 1 if STATE == "J"
	sum 		pi_J 
	replace		pi = r(mean) if STATE == "J"	

*collapse
	collapse lambda pi, by(STATE)
	gen TYPE = "`SAMPLE_TYPE'"
	append using "`dir_tmp_aux'/penalty"
	save "`dir_tmp_aux'/penalty",replace
	
*now women
	local SAMPLE_TYPE = "FEMALE"
	local J_J_STRINGA = "1990_2016_2005_2016_RS_build"	
	local file_type ="`SAMPLE_TYPE'_SAMPLE__INVIND_data_`J_J_STRINGA'"
	local dir_tmp_aux = "/scratch/public/leave_out"	

*import the file with the estimates from the DAKM model
	import delimited "`dir_tmp_aux'/`file_type'.csv", encoding(ISO-8859-1) clear
	rename v2 id
	rename v3 firmid
	rename v4 lagfirmid
	rename v7 lambda
	
*gen the three states
	gen STATE = ""
	replace STATE = "U" if lagfirmid == 0
	replace STATE = "N" if lagfirmid == -1
	replace STATE = "J" if lagfirmid >0
	
*normalize
	sum lambda if STATE == "N"
	replace lambda = lambda-r(mean)
	
*show me the ATT
if 1 == 1 {
	preserve
	drop lambda
	gen rows = _n
	bys id (rows): gen contami = _n
	xtset id contami
	bys id: gen firmid_aux = firmid[_n-1]
	drop firmid
	gen firmid = firmid_aux
	merge m:1 firmid using "`dir_tmp_aux'/fe_`file_type'", nogen keepusing(lambda DIP10)
	keep if STATE == "U"
	sum lambda
	restore	
}

*show me the ATU
if 1 == 1 {
	preserve
	drop lambda
	gen rows = _n
	bys id (rows): gen contami = _n
	xtset id contami
	bys id: gen firmid_aux = firmid[_n-1]
	drop firmid
	gen firmid = firmid_aux
	merge m:1 firmid using "`dir_tmp_aux'/fe_`file_type'", nogen keepusing(lambda DIP10)
	keep if STATE == "J"
	sum lambda
	restore	
}

*report the shares
	gen 		pi = .
	gen 		pi_U = 0
	replace		pi_U = 1 if STATE == "U"
	sum 		pi_U 
	replace		pi = r(mean) if STATE == "U"
	gen 		pi_J = 0
	replace		pi_J = 1 if STATE == "J"
	sum 		pi_J 
	replace		pi = r(mean) if STATE == "J"	

*collapse
	collapse lambda pi, by(STATE)
	gen TYPE = "`SAMPLE_TYPE'"
	append using "`dir_tmp_aux'/penalty"
	save "`dir_tmp_aux'/penalty",replace
	list
	
*compute b/w component
	use "`dir_tmp_aux'/penalty",replace
	drop if STATE == "N"
	gen stato = 0
	replace stato = 1 if STATE == "J"
	reshape wide lambda pi STATE, i(TYPE) j(stato)
	list
	gen between = pi0*(1-pi0)*lambda0^2+pi1*(1-pi1)*lambda1^2-2*pi1*pi0*lambda0*lambda1
	gen totale  = .
	replace totale= 0.0019232 if TYPE == "POOLED"
	replace totale= 0.0020615 if TYPE == "MALE"
	replace totale= 0.0018603 if TYPE == "FEMALE"
	gen within  = (totale-between)/pi1
	replace within = within^0.5
	replace between = between^0.5
	replace totale = totale^0.5
	list

		
