*read the global
	local file_type = "$file_type"
	local dir_tmp_aux = "/scratch/public/leave_out"	
	local gg = 2
*********************************************************************************

*The idea here is to take the baseline leave out connected sample, attach
*the manager dummy and export a .CSV which is going to be feed into the build
*routine in Matlab

*********************************************************************************
*tell me manager infomation
if 0 == 1 {
	use "`dir_tmp_aux'/DAVE_BUILD_RAKIM_APPENDED100_INVIND_data_1990_2016_2005_2016_RS_build",replace
	replace log_dailywages=round(log_dailywages,0.00001)
	gen			manager   = 0 
	replace		manager   = 1 		if qualifica == "Q" | qualifica == "3"
	gen 		age2	  = ((age-40)/40)^2
	gen 		age3	  = ((age-40)/40)^3
	bys log_dailywages firmid lagfirmid age2 age3 year: gen TTT=_N
	keep if TTT==1
	collapse manager, by(log_dailywages firmid lagfirmid age2 age3 year)
	save "`dir_tmp_aux'/managers",replace
}
	
*Load pruned sample	
	import delimited "`dir_tmp_aux'/POOLED_SAMPLE__INVIND_data_1990_2016_2005_2016_RS_build.csv", encoding(ISO-8859-1) clear
	gen sortami=_n
	list v1 v2 v3 v4 in 1/50
	gen log_dailywages=v1
	gen id=v2
	gen firmid=v3
	gen lagfirmid=v4
	gen age2=v19
	gen age3=v20
	gen year = 2004
	local aux = 1
	forval pp = 9(1)18 {
	replace year = 2004 + `aux' if v`pp' == 1
	local aux = `aux' +1
	}
	replace year = 2015 if year == 2004
	sum year
	replace log_dailywages=round(log_dailywages,0.00001)
	merge m:1 log_dailywages firmid lagfirmid age2 age3 year using "`dir_tmp_aux'/managers", gen(identified) keepusing(manager)
	keep if identified == 1 | identified == 3
	replace manager =0 if manager==.
	drop if sortami == .
*Export
	sort sortami
	sum manager
	export delimited manager  using "`dir_tmp_aux'/MANAGER", replace novarnames nolabel
