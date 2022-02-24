********************************************************************************
* TELL ME WHAT SPECIFICATION YOU ARE RUNNING AND WOULD LIKE TO GET TABLES OUT OF
********************************************************************************
local dir_tmp_aux "/scratch/public/leave_out"	
local jj_flag = $jjSPEC

if `jj_flag' == 3{
	local J_J_STRINGA = "group_1_2_3"
}

if `jj_flag' == 2{
	local J_J_STRINGA = "group_1_2"	
}

if `jj_flag' == 1{
	local J_J_STRINGA = "group_1"	
}

if `jj_flag' == 4{
	local J_J_STRINGA = "group_1_2_3_no_diff"	
}

if `jj_flag' == 5{
	local J_J_STRINGA = "group_1_2_3_BEAUDRY_DI_NARDO"	
}
********************************************************************************
* CALCULATE THE MEANS
********************************************************************************
*load the gender specific samples	
 	forval pp = 1(1)3 {
	 
	 if `pp' == 3{
	 	local file_type = "FEMALE_SAMPLE"
	 }

	 if `pp' == 2{
		local file_type = "MALE_SAMPLE"
	 }

	 if `pp' == 1{
		local file_type = "POOLED_SAMPLE"	
	 }
*load
	use "`dir_tmp_aux'/PVR_`file_type'_ESTIMATES_WITH_CROSS_WALK",replace 	 
	
*keep identified only
	keep if identified_psi == 1 & identified_lambda ==1 
	
*normalize lambda
	sum lambda if firmid == -1
	replace lambda = lambda - r(mean)
	sum lambda 
	local MEDIA_lambda_`pp' = r(mean)

*normalize psi
	merge 1:1 firmid using "`dir_tmp_aux'/firm_size_full_file", nogen
	xtile vinG = log_firm_size, nquantile(20)
	sum psi 	if vinG ==1
	replace psi=psi-r(mean)
	sum psi 
	local MEDIA_psi_`pp' = r(mean)
	
}
stop
********************************************************************************
* PIN DOWN THE TABLE
********************************************************************************
 forval tt = 1(1)4 {
	
	if `tt' ==  1{
	 	local tableToLoad = "tableSum"
	}

	if `tt' ==  2{
	 	local tableToLoad = "tableR2"
	}

	if `tt' ==  3{
	 	local tableToLoad = "tableDecomp"
	}
	 
	if `tt' ==  4{
	 	local tableToLoad = "tableFirmMoments"
	}
 

*import the template that Mikkel created
	import delimited "tables/`tableToLoad'_template.csv", delimiter(";") varnames(1) encoding(ISO-8859-1) clear
	gen rows = _n +1 

	save "tables/template_Mikkel",replace

*now read the file from matlab
	 forval pp = 1(1)3 {
	 
	 if `pp' == 3{
	 	local typeSAMPLE = "FEMALE_SAMPLE"
		local myName	 = "female"
	 }

	 if `pp' == 2{
		local typeSAMPLE = "MALE_SAMPLE"
		local myName	 = "male"
	 }

	 if `pp' == 1{
		local typeSAMPLE = "POOLED_SAMPLE"	
		local myName	 = "all"
	 }
	 
	 
    import delimited "tables/TABLE_`tt'`typeSAMPLE'_`J_J_STRINGA'.csv", delimiter(";") encoding(ISO-8859-1) clear
	replace v1 = . if v1 == -9999999
	tostring v1, generate(v2) format(%12.4fc) force
	tostring v1, generate(v3) format(%9.4gc) force
	gen `myName' = v2 if v1<100
	replace `myName' = v3 if v1>100
	replace `myName' = "" if `myName' == "."
	drop v1 v2 v3
	if `tt' == 1 {
	gen rows=_n
	}
	
	if `tt' > 1 {
	gen rows=_n+1
	}

	if `tt' == 2 & `pp' > 1 {
	replace `myName' = "" if rows== 12 | rows == 7 // don't show R2 interacted model in female and male samples only
	}
	
	
	save "tables/aux_`pp'",replace	 
	}
	
	
*concatenate the template and the matlab tables
	use "tables/template_Mikkel",replace
	forval pp = 1(1)3 {
	merge 1:1 rows using "tables/aux_`pp'", nogen
	}

	
*minor fixes	
	drop rows
	drop if type == ""
	rename followtitle followTitle
	
	
*ciao
	export delimited using "tables/`tableToLoad'_`J_J_STRINGA'_RAFFA", datafmt replace  delimiter(";") 

}

