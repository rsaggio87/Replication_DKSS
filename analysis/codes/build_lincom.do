local dir_tmp_aux = "/scratch/public/leave_out" 
*********************************************************************************

*	EXPORT A CSV FILE THAT FOR EACH DESTINATION FIRM IT PROVIDES ME WITH SIZE AND LOG VA/L 
*********************************************************************************
if 1 == 1{
	use "`dir_tmp_aux'/va_file",replace
	keep if log_VA_L!=.
	export delimited firmid DIP10 log_VA_L using "`dir_tmp_aux'/LINCOM", replace novarnames nolabel
}
	
