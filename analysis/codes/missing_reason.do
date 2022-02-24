*read the global
	local file_type = "$file_type"
	local dir_tmp_aux = "/scratch/public/leave_out"	
	local nsamples = 1
	global var_regression i.year south female age maxTenure temp_contract firmsize log_dailywages closed_firm ended_Dec
	
*********************************************************************************

*	START BY LOADING THE MONTHLY FILES, KEEP LAST MONTH BEFORE SEP (PARALLEL)

*********************************************************************************
if 0 == 1 {
forval ss=1(1)`nsamples'{
*LOAD THE MONTHLY PANEL
		use "/scratch/public/board/RAKIM/CLEANED_MONTHLY_PANEL_GIANT_`ss'_INVIND_data_1990_2016", replace 
			
*CREATE THE TENURE VARIABLE	
		xtset id date
		bys id id_impresa working (date): gen tenure =_n
		replace tenure = 0 if working == 0			

*KEEP LAST MONTH< JOB MUST LAST AT LEAST THREE MONTHS, SEPARATION MUST OCCUR IN 2002
		bys id id_impresa: egen maxTenure = max(tenure)
		keep if tenure == maxTenure & maxTenure>=3 & year>=2005 // notice that I am not doing much about cycles of employment: started at firm A then moved to firm B then back at Firm A...
		xtset id date							
		  
*NOW UNDERSTAND WHETHER IT WAS AN HIRE FROM UNEMPLOYMENT OR AN HIRE FROM EMPLOYMENT
		xtset id date
		bys id: gen gap = date[_n] -  date[_n-1]
		replace gap = . if tenure !=1
		sum gap,d
			    
*CLEAN THE REASON VARIABLE 
		gen 	reason = .
		replace reason = 1 if motivo_cessazione == "1A" // objective reason
		replace reason = 2 if motivo_cessazione == "1B" // resign
		replace reason = 3 if motivo_cessazione == "1C" // End of temp contract
		replace reason = 4 if motivo_cessazione == "1D" // Subjective reason
		replace reason = 5 if motivo_cessazione != "" & motivo_cessazione != "1A"  & motivo_cessazione != "1B" & motivo_cessazione != "1C" & motivo_cessazione != "1D" 
		drop motivo_cessazione
		
				 	
*FINAL TOUCHES
		gen 	   female 			 = 0
		replace    female 			 = 1 		 	 if sesso == "F"
		gen		   firmsize			 = log(DIP10)
		
		tostring id_azienda, replace
		replace id_impresa = id_azienda if id_impresa ==""

		bys id:  gen conteggio = _n
		xtset id date

*DROP CASES DUE TO RIGHT TRUNCATION (sometimes there is recall but can't fully observe that in 2015 because the data ends there)
		sum year
		drop if year ==r(max)
		
*GENERATE VARIABLE OF INTEREST
		gen missing_reason = .
		replace missing_reason = 0 if reason !=. 
		replace missing_reason = 1 if reason ==. 		
		
*SAVE
		if `ss' == 1  {				
		save "`dir_tmp_aux'/missing_analysis", replace 
		}
		
		if `ss'>1{
		append using "`dir_tmp_aux'/missing_analysis"
		save "`dir_tmp_aux'/missing_analysis", replace 
		}		
}

*VERIFY TIME-SERIES PATTERN
		preserve
		collapse missing_reason, by(year)
		twoway (connected missing_reason year, msymbol(square) mcolor(dknavy) msize(medlarge)), ytitle("Missing Reason for Separation") xtitle("Year") ylabel(0.15(0.01)0.20) xlabel(2005(1)2014)
		graph export "figures/missing_reason_balanced_n_of_sample_`nsamples'.pdf",replace
		restore	

}
**************************************************************************

*	 PREPARE THE REGRESSION
	
**************************************************************************
if 1 == 1 {	
*LOAD 
	use "`dir_tmp_aux'/missing_analysis", replace 
	keep if missing_reason!=.
	
*Dummy for whether job ended in December
	replace month = month(dofm(date)) 		
	gen ended_Dec=0
	replace ended_Dec=1 if month==12
	
	
*Closing firm
if 0 == 1 {
	preserve
	use CF_PERSG ID_PERSF D_CESSAZ year using "/scratch/public/leave_out/matr_bitalia_1990_2013NEW_INVIND", replace	
	sum year
	local MASSIMO=r(max)
	gen id_impresa = CF_PERSG
	tostring ID_PERSF,replace
	replace id_impresa= ID_PERSF if CF_PERSG==""
	collapse (max) firm_closure=year, by(id_impresa)
	drop if firm_closure == `MASSIMO'
	save "/scratch/public/leave_out/closing_firms",replace
	restore
}	
	merge m:1 id_impresa using "/scratch/public/leave_out/closing_firms", gen(firm_close)
	drop if firm_close==2
	gen closed_firm =0
	replace closed_firm=1 if year==firm_closure 
	
*Auxiliaries
	gen temp_contract = 0
	replace temp_contract = 1 if tipo_contratto == "D" | tipo_contratto == "I" & qualifica == "5" | tipo_contratto == "I" & qualifica == "4"
	rename provincia_lavoro provincia_recovered
	merge m:1 provincia_recovered using "/scratch/public/leave_out/PROVINCE",nogen
	gen south =0
	replace south =1 if descrizione_regione == "SICILIA" | descrizione_regione == "CAMPANIA" | descrizione_regione == "SARDEGNA" | descrizione_regione == "PUGLIA" | descrizione_regione == "CALABRIA" | descrizione_regione == "BASILICATA" | descrizione_regione == "MOLISE" | descrizione_regione == "ABRUZZO"


*regression table
	reg missing_reason $var_regression, cluster(id_impresa)
	eststo rr_1, title("All Jobs in INPS-INVIND")
	sum missing_reason if e(sample) == 1
	estadd local MEDIA=round(r(mean),0.01),replace
	estadd local adjR2=round(e(r2_a),0.01),replace
	
*run the simplest IV regression now
	reg missing_reason $var_regression if temp_contract==1, cluster(id_impresa)
	eststo rr_2, title("All Temp Jobs in INPS-INVIND")
	sum missing_reason if e(sample) == 1
	estadd local MEDIA=round(r(mean),0.01),replace
	estadd local adjR2=round(e(r2_a),0.01),replace
	
*start adding firm level controls
	reg missing_reason $var_regression if temp_contract==0, cluster(id_impresa)
	eststo rr_3, title("All Non-Temp Jobs in INPS-INVIND")
	sum missing_reason if e(sample) == 1
	estadd local MEDIA=round(r(mean),0.01),replace
	estadd local adjR2=round(e(r2_a),0.01),replace
	
*start adding firm level controls: polynomial
	reg missing_reason $var_regression if ended_Dec==0, cluster(id_impresa)
	eststo rr_4, title("All Jobs not terminated in December")
	sum missing_reason if e(sample) == 1
	estadd local MEDIA=round(r(mean),0.01),replace
	estadd local adjR2=round(e(r2_a),0.01),replace
				
*export the IV regression
	esttab rr_* using "tables/missing_reason.csv", label wide star(* 0.10 ** 0.05 *** 0.01) cells(b(star fmt(%9.4f)) se(par fmt(%9.4f))) stats(N MEDIA adjR2, fmt(%9.0f %9.4f  %9.4f) labels("Observations" "Mean Outcome" "Adj. R2")) replace
}		