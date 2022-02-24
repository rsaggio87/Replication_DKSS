local dir_tmp_aux = "/scratch/public/leave_out" 
*********************************************************************************

*	EXPORT A CSV FILE THAT FOR EACH DESTINATION FIRM IT PROVIDES ME WITH SIZE AND SECTOR 
*********************************************************************************
if 1 == 1{
	import delimited "`dir_tmp_aux'/POOLED_SAMPLE__INVIND_data_1990_2016_2005_2016_RS_build.csv", encoding(ISO-8859-1) clear
	rename v1 log_dailywages
	rename v3 firmid
	collapse log_dailywages, by(firmid)
	drop  log_dailywages
	merge 1:m firmid using "`dir_tmp_aux'/DAVE_BUILD_RAKIM_APPENDED100_INVIND_data_1990_2016_2005_2016_RS_build", keepusing(id_impresa) gen(cross_walk)
	merge m:1 id_impresa using "`dir_tmp_aux'/firms_universe_upd", gen(my_sample) keepusing(settore)
	keep if my_sample == 3
	collapse my_sample settore, by(firmid)
	drop my_sample
	merge m:1 firmid using "`dir_tmp_aux'/va_file", gen(my_sample) keepusing(DIP10)
	keep if my_sample==3
	replace DIP10 = ceil(DIP10)
	sum DIP10,d
	export delimited firmid DIP10 settore using "`dir_tmp_aux'/FIRM_SIZE.csv", replace novarnames nolabel
	save "`dir_tmp_aux'/FIRM_SIZE",replace
}
	