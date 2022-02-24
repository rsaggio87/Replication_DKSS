
*read the global
	local dir_tmp_aux = "/scratch/public/leave_out"	
	local file_type   = "$file_type"
	local female	  = $FEMALE
	
	if `female' == 0 {
	local type_string = "_male"
	}
	
	if `female' == 1 {
	local type_string = "_female"
	}

*save auxiliary files
	use "`dir_tmp_aux'/PVR_`file_type'_ESTIMATES_WITH_CROSS_WALK",replace 
	
	preserve
	keep psi* identified_psi* firmid
	rename psi psi`type_string'
	rename identified_psi identified_psi`type_string'
	save "`dir_tmp_aux'/PSI_PVR_`file_type'_ESTIMATES_WITH_CROSS_WALK",replace 
	restore
	
	preserve
	keep lambda* identified_lambda* firmid
	rename firmid firmid_lag
	rename lambda lambda`type_string'
	rename identified_lambda identified_lambda`type_string'
	save "`dir_tmp_aux'/LAMBDA_PVR_`file_type'_ESTIMATES_WITH_CROSS_WALK",replace 
	restore		