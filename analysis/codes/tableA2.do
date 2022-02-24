*read the global
	local dir_tmp_aux = "../build/src/RAKIM/"	
	global var_cuts log_dailywages
	local nsamples=100
	local tenure_restriction=12
	local prev_wage_to_use = 1
	local margini_y=2
	local margin_x=1
	local margini_y_=`margini_y'-0.05
	
*read locals	
	if `tenure_restriction'>3 {
	local tenure_restriction_NAME "tenure_restric_`tenure_restriction'"
	}	
*********************************************************************************

*	BUILD THE SAMPLE, NEED TO BRING INFO WAGE ON PREV JOB

*********************************************************************************
if 1 == 1 {
forval ss=1(1)`nsamples'{
*LOAD THE MONTHLY PANEL
		use  "../build/src/RAKIM/CLEANED_MONTHLY_PANEL_GIANT_`ss'_INVIND_data_1990_2016", replace 	

*CONVERT TO NOMINAL
		bys year: sum log_dailywages
		gen cpi = 0 
		replace cpi = 55.8390000000000000 if year == 1990
		replace cpi = 59.3289000000000000 if year == 1991
		replace cpi = 62.4559000000000000 if year == 1992
		replace cpi = 65.3456000000000000 if year == 1993
		replace cpi = 67.9933000000000000 if year == 1994
		replace cpi = 71.5530000000000000 if year == 1995
		replace cpi = 74.4201000000000000 if year == 1996
		replace cpi = 75.9406000000000000 if year == 1997
		replace cpi = 77.4253000000000000 if year == 1998
		replace cpi = 78.7133000000000000 if year == 1999
		replace cpi = 80.7107000000000000 if year == 2000
		replace cpi = 82.9587000000000000 if year == 2001
		replace cpi = 85.0039000000000000 if year == 2002
		replace cpi = 87.2757000000000000 if year == 2003
		replace cpi = 89.2016000000000000 if year == 2004
		replace cpi = 90.9725000000000000 if year == 2005
		replace cpi = 92.8746000000000000 if year == 2006
		replace cpi = 94.5740000000000000 if year == 2007
		replace cpi = 97.7402000000000000 if year == 2008
		replace cpi = 98.4974000000000000 if year == 2009
		replace cpi = 100 if year == 2010
		replace cpi = 102.7810000000000000 if year == 2011
		replace cpi = 105.9070000000000000 if year == 2012
		replace cpi = 107.1990000000000000 if year == 2013
		replace cpi = 107.4570000000000000 if year == 2014
		replace cpi = 107.4990000000000000 if year == 2015
		replace cpi = 107.3980000000000000 if year == 2016
		replace log_dailywages=exp(log_dailywages)
		replace log_dailywages=(log_dailywages*cpi)/100
		replace log_dailywages=log(log_dailywages)
		drop cpi
		bys year: sum log_dailywages
			
*CREATE THE TENURE VARIABLE	
		xtset id date
		bys id id_impresa working (date): gen tenure =_n
		replace tenure = 0 if working == 0	
		bys id id_impresa: egen maxTenure = max(tenure)
		
*CALCULATE WAGE IN MATCH
		by id id_impresa: egen avg_wage =mean(log_dailywages)
		gen last_wage=log_dailywages if tenure==maxTenure
		gen first_wage=log_dailywages if tenure==1
		gen job12_wage=log_dailywages if tenure==maxTenure-12
		bys id id_impresa: egen last_WAGE=mean(last_wage)
		bys id id_impresa: egen first_WAGE=mean(first_wage)
		bys id id_impresa: egen job12_WAGE=mean(job12_wage)
		
*FOR RAKIM I ONLY NEED START AND END OF A TENURE
		keep if tenure == 1 | tenure == maxTenure & tenure> 0 // notice that I am not doing much about cycles of employment: started at firm A then moved to firm B then back at Firm A...
		
*KEEP JOBS
		keep if maxTenure>=3											
		  
*IT'S IMPORTANT TO UNDERSTAND THE TENURE VARIABLE
		sum maxTenure,d
			    
*NOW UNDERSTAND WHETHER IT WAS AN HIRE FROM UNEMPLOYMENT OR AN HIRE FROM EMPLOYMENT
		xtset id date
		bys id: gen gap = date[_n] -  date[_n-1]
		replace gap = . if tenure !=1
		sum gap,d
			    
*BRING FORWARD THE INFORMATION ON WHY THE JOB WAS DESTROYED (SO NOW IT SHOWS WHEN TENURE==1). Notice that for jobs that lasted for less than a year you will see the information on separation both for rows when tenure=1 and when tenure=MaxTenure, see here: https://www.dropbox.com/s/aqsw4wz1wou341d/screen_RAKIM.png?dl=1 
		xtset id date
		bys id: gen MOTIVO = motivo_cessazione[_n-1]
		bys id: gen past_wage=avg_wage[_n-1]
		bys id: gen past_last_wage=last_WAGE[_n-1]
		bys id: gen past_first_wage=first_WAGE[_n-1]
		bys id: gen past_12_wage=job12_WAGE[_n-1]
		
		replace gap = . if tenure !=1
		sum gap,d	
		 
*CLEAN THE REASON VARIABLE 
		gen 	reason = .
		replace reason = 1 if MOTIVO == "1A" // objective reason
		replace reason = 2 if MOTIVO == "1B" // resign
		replace reason = 3 if MOTIVO == "1C" // End of temp contract
		replace reason = 4 if MOTIVO == "1D" // Subjective reason
		replace reason = 5 if MOTIVO != "" & MOTIVO != "1A"  & MOTIVO != "1B" & MOTIVO != "1C" & MOTIVO != "1D" 
		drop MOTIVO motivo_cessazione
		
				 	
*FINAL TOUCHES
		gen 	   female 			 = 0
		replace    female 			 = 1 		 	 if sesso == "F"
		gen		   firmsize			 = log(DIP10)	
		tostring id_azienda, replace
		replace id_impresa = id_azienda if id_impresa ==""
		egen firmid 	= group(id_impresa)
		bys id: gen lagfirmid = firmid[_n-1]
		bys id:  gen conteggio = _n
		replace lagfirmid 	  =	-1 if conteggio==1		
		
*DROP JOBS FOR WHICH I HAVE A RIGHT CENSORING PROBLEM (CAN'T REALLY SAY WHETHER LAGFIRM = 0 or not)
		gen 		right_cen 		 = 0
		replace		right_cen		 = 1			 if conteggio == 1 & year == 1990
		drop if right_cen == 1 // notice that by imposing this condition, there is going to be some individuals (those who had a job in 1990) for which gap will always be non-missing. 		
		
*DROP DUDES OBSERVED WITH ONLY ONE JOB
		distinct id
		bys id: gen TTT = _N
		drop if TTT == 1
		distinct id
		drop TTT
				
*DEFINE LAGGED STATE
		xtset id date
		drop lagfirmid conteggio
		gen lagfirmid = 0
		bys id: gen pre_firm = firmid[_n-1]
		bys id:  gen conteggio = _n
		replace lagfirmid 	  =	-1 				if conteggio==1
		replace lagfirmid 	  =  pre_firm		if reason == 2 & conteggio > 1 | gap==1 & reason == . & conteggio>1	
		
*Restrict to sample of jobs that started in 2005 or forward
		keep if tenure == 1 & year>=2005
					
*SAVE
		if `ss' == 1  {				
		save "`dir_tmp_aux'/cuts_analysis", replace 
		}
		
		if `ss'>1{
		append using "`dir_tmp_aux'/cuts_analysis"
		save "`dir_tmp_aux'/cuts_analysis", replace 
		}		
}

		save "`dir_tmp_aux'/cuts_analysis", replace 
}

*********************************************************************************

*	MERGE VA

*********************************************************************************
if 1 == 1{
*LOAD
		use  "`dir_tmp_aux'/cuts_analysis", replace 
		capture drop codfis_mod merge_VA VA_L valagg

		
*BRING INFO ON VA
		destring id_impresa, gen(codfis_mod) force
		tostring codfis_mod, format(%016.0f) replace
		merge m:1 codfis_mod year using "`dir_tmp_aux'/cerved_smaller",gen(merge_VA)
		keep if merge_VA == 1 | merge_VA == 3
		replace DIP10 = 1 if DIP10 == 0
*NOMINAL
		gen cpi = 0 
		replace cpi = 55.8390000000000000 if year == 1990
		replace cpi = 59.3289000000000000 if year == 1991
		replace cpi = 62.4559000000000000 if year == 1992
		replace cpi = 65.3456000000000000 if year == 1993
		replace cpi = 67.9933000000000000 if year == 1994
		replace cpi = 71.5530000000000000 if year == 1995
		replace cpi = 74.4201000000000000 if year == 1996
		replace cpi = 75.9406000000000000 if year == 1997
		replace cpi = 77.4253000000000000 if year == 1998
		replace cpi = 78.7133000000000000 if year == 1999
		replace cpi = 80.7107000000000000 if year == 2000
		replace cpi = 82.9587000000000000 if year == 2001
		replace cpi = 85.0039000000000000 if year == 2002
		replace cpi = 87.2757000000000000 if year == 2003
		replace cpi = 89.2016000000000000 if year == 2004
		replace cpi = 90.9725000000000000 if year == 2005
		replace cpi = 92.8746000000000000 if year == 2006
		replace cpi = 94.5740000000000000 if year == 2007
		replace cpi = 97.7402000000000000 if year == 2008
		replace cpi = 98.4974000000000000 if year == 2009
		replace cpi = 100 if year == 2010
		replace cpi = 102.7810000000000000 if year == 2011
		replace cpi = 105.9070000000000000 if year == 2012
		replace cpi = 107.1990000000000000 if year == 2013
		replace cpi = 107.4570000000000000 if year == 2014
		replace cpi = 107.4990000000000000 if year == 2015
		replace cpi = 107.3980000000000000 if year == 2016
		replace valagg=(valagg*cpi)/100

		gen VA_L = valagg/DIP10	
		capture drop log_VA_L
		gen log_VA_L = log(VA_L)
		capture drop MEDIA
		bys id_impresa: egen MEDIA=mean(log_VA_L)
		replace   	log_VA_L=MEDIA
		sum 		log_VA_L,d
		replace 	log_VA_L =r(p5)  if log_VA_L<r(p5)
		replace 	log_VA_L =r(p95) if log_VA_L>r(p95) & log_VA_L!=. 
		xtset id date
		capture drop past_VA
		bys id: gen past_VA=log_VA_L[_n-1]
		capture drop codfis_mod merge_VA VA_L valagg
		save "`dir_tmp_aux'/cuts_analysis_use", replace 
}
*********************************************************************************

*	BASELINE STATS

*********************************************************************************
if 1 == 1 {
*load
			use gap log* past* firmid lagfirmid maxTenure tipo_contratto qualifica if maxTenure>=`tenure_restriction' using "`dir_tmp_aux'/cuts_analysis_use", replace
			forval prev_wage_to_use=1(1)3{
			preserve
			if `prev_wage_to_use'==1 {
				local prev_wage_use_NAME "avg_wage"
				local USAMI "past_wage"
				local labella "avg wage in prev. job"
			}
	
	
			if `prev_wage_to_use'==4 {
				local prev_wage_use_NAME "last_Wage"
				local USAMI "past_last_wage"
				local labella "last wage in prev. job"
			}
	
			if `prev_wage_to_use'==2 {
				local prev_wage_use_NAME "first_wage"
				local USAMI "past_first_wage"
				local labella "poaching wage in prev. job"
				local margini_y=1
				local margin_x=2
				local margini_y_=`margini_y'-0.05
			}
	
	
			if `prev_wage_to_use'==3 {
				local prev_wage_use_NAME "wage_12"
				local USAMI "past_12_wage"
				local labella "wage in prev. job, 12 months before separation"
	
			}


*temp contract
			gen 			temp_contract = 0
			replace 		temp_contract = 1 if tipo_contratto == "D" | tipo_contratto == "I" & qualifica == "5" | tipo_contratto == "I" & qualifica == "4"

*change in wages 		
			gen			change_wage = log_dailywages-`USAMI'
			gen 		change_VA=log_VA_L-past_VA
		
			gen 		share_negative_wage = .	
			replace 	share_negative_wage = 0 if change_wage>=0 & change_wage!=.
			replace 	share_negative_wage = 1 if change_wage<0
		
			gen 		share_negative_VA = .	
			replace 	share_negative_VA = 0 if change_VA>=0 & change_VA!=. | change_VA==. & log_VA_L!=. & past_VA==.
			replace 	share_negative_VA = 1 if change_VA<0 | change_VA==. & log_VA_L==. & past_VA!=.

			sum 		share_negative_VA  if gap == 1
			sum 		share_negative_VA  if gap >1 & gap!=.
		
			gen  		treatment = 0 
			replace		treatment = 1 		if lagfirmid>0	
		
		
			collapse  share_negative_wage share_negative_VA, by(treatment temp_contract)
			rename share_negative_wage share_negative_wage_`prev_wage_to_use'
			save "tables/baseline_stats`prev_wage_to_use'",replace
			restore
			}
			
			use "tables/baseline_stats1",replace
			merge 1:1 treatment temp_contract using "tables/baseline_stats2",nogen
			merge 1:1 treatment temp_contract using "tables/baseline_stats3",nogen
			rename treatment voluntary
			order share_negative_wage_2 share_negative_wage_1 share_negative_wage_3 share_negative_VA voluntary temp_contract
			sort temp_contract voluntary
			list
			save "tables/baseline_stats",replace
			
}
