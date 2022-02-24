*read the global
	local file_type = "$file_type"
	local dir_tmp_aux = "/scratch/public/leave_out"	
	
*********************************************************************************

*	BUILD THE SAMPLE, NEED TO BRING INFO WAGE ON PREV JOB

*********************************************************************************
if 0 == 1 {
*LOAD THE MONTHLY PANEL
		use "/scratch/public/board/RAKIM/CLEANED_SPELL_GIANT_INVIND_data_1990_2016", replace 
		describe
			
*CREATE THE TENURE VARIABLE	
		bys id id_impresa year: gen n_spells_year=_N 
		bys id id_impresa year: gen count_spells_year=_n
		bys id id_impresa: egen min_year =min(year)	
		

		gen 			temp_contract = 0
		replace 		temp_contract = 1 if tipo_contratto == "D" | tipo_contratto == "I" & qualifica == "5" | tipo_contratto == "I" & qualifica == "4"
		bys id id_impresa year: egen on_a_temp_contract=max(temp_contract)
		keep if  year==min_year	
		save "`dir_tmp_aux'/hiring_year", replace 
		}



*LOAD
		use "`dir_tmp_aux'/hiring_year", replace 
		sum n_spells_year if count_spells_year==1, d
		sum n_spells_year if count_spells_year==1 & on_a_temp_contract==0, d
		
		capture drop more_than_one
		gen 	more_than_one=0
		replace more_than_one=1 if n_spells_year>1
		
		sum more_than_one if count_spells_year==1
		sum more_than_one if on_a_temp_contract==0 & count_spells_year==1
		
*UNDERSTAND WHICH SPELL CAME FIRST
		gen start_month = 12
		forval t=11(-1)1{
		replace start_month = `t' if decl`t'==1
		}

*NOW SORT
		sort id year id_impresa start_month
		keep if n_spells_year>=2
		gen log_daily_wages=log(imponibile/giorni_retribuiti)
		bys id year id_impresa (start_month): gen conteggio=_n
		gen ssrtar_wage=log_daily_wages if conteggio==1
		bys id year id_impresa: egen START_WAGE = mean(ssrtar_wage)

*NOW TELL ME WAGE INCREASES
		gen 	 wage_increase = .
		replace  wage_increase = 0 if log_daily_wages<START_WAGE & conteggio>1
	    replace  wage_increase = 1 if log_daily_wages==START_WAGE & conteggio>1
		replace	 wage_increase = 2 if log_daily_wages>START_WAGE & conteggio>1
		tab wage_increase
		tab wage_increase  if on_a_temp_contract==0
		
}		

