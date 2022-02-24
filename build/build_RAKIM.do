********************************************************************************

						*DECLARE AND INSTALL* 	
					
********************************************************************************
set seed 1234
set more off
local small_sample		= 0 //<---------------------------Instead of working with the enormous file take only,say, 1000 workers from each file
local ss				= $it
local total_files_app	= 100
if `small_sample' == 1 {
	local sample_load = "small"
}
local type_data "INVIND_data_1990_2016`sample_load'"
********************************************************************************

*	SECTION 1: LOAD INITIAL FILE
					
********************************************************************************
if 1 == 1 {
*LOAD SRC FILE 
		use "src/giantNEW_INVIND`sample_load'.dta",replace
		rename mt_qdrb_anno year
		
*Drop variables we are not going to use for sure
		drop mt_qdrb_assicopalt mt_qdrb_assicopds mt_qdrb_assicopivs mt_qdrb_assicopivs  mt_qdrb_codcomune mt_qdrb_partcf_min  assco asscofg flag_anf trasfrl tfr_euro_tc trasfrl numero_anf classe_anf tabella_anf vet_insurance unemp_insurance other_insurance tax_exem miner

*BRING IN INFORMATION ON FISCAL CODE
		merge m:1 id_azienda year using "src/matr_bitalia_1990_2013NEW_INVIND.dta", keepusing(prov codice_contributivo CF_PERSG ID_PERSF DIP10) gen(merge_prov)		
		drop if merge_prov == 2
		
*ENLARGE LIST OF ID WORKER
		format %20.0g id

*DROP DUPLICATES
		duplicates drop id year id_azienda mt_qdrb_ggretrib mt_qdrb_settimane mt_qdrb_imponibile, force
		
**************************************************************************

*	SECTION 2: CROSS WALK INVIND DATA
	
**************************************************************************
		rename id id_soggetto
		rename C_CONTRATTO codice_contratto
		gen livello_inquadramento = .
		rename mt_qdrb_qualif2 tipo_rapporto 
		rename mt_qdrb_qualif1 qualifica 
		rename mt_qdrb_qualif3 tipo_contratto 
		rename codcomune comune_residenza
		rename mt_qdrb_imponibile imponibile
		rename mt_qdrb_ggretrib giorni_retribuiti
		rename mt_qdrb_settimane settimane_retribuite
		rename mt_qdrb_settutili settimane_utili
		rename mt_qdrb_tipcontri tipo_contribuzione
		rename provincia provincia_residenza
		drop provlav		
		gen id_impresa = CF_PERSG
		tostring ID_PERSF,replace
		replace id_impresa= ID_PERSF if CF_PERSG==""
		keep if merge_prov==3 | merge_prov ==1
		drop merge_prov
		rename prov provincia_lavoro
		rename codice_contributivo industry_code
		gen comune_lavoro = .
		rename mt_qdrb_assunzione data_assunzione
		rename mt_qdrb_tipoassunz motivo_assunzione
		rename mt_qdrb_cessazione data_cessazione 
		rename mt_qdrb_tipo_cessazione motivo_cessazione 
		rename stato_nascita cittadinanza //will have to do one more fix on this
		gen anno_pensionamento =.
		gen mese_pensionamento =.
		gen provincia_contribuzione=. 
		gen tipo_politica = . 
		egen id=group(id_soggetto)
		bys id: egen anno_inizio_lavoro=min(year)
**************************************************************************

*	SECTION 3: CLEAN SPELLS
	
**************************************************************************			
*CPI ADJUSTED EARNINGS
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
		replace imponibile=(imponibile*100)/cpi
		gen daily_wage=imponibile/ giorni_retribuiti
		gen log_dailywages=log(daily_wage)
		drop cpi
		
*WATERFALL OF RESTRICTIONS:
		capture file close myfile	
		file open myfile using "src/Observations_Cleaning_SPELL_monthly_`type_data'.txt", text write replace
		
*-------*BEGIN RESTRICTIONS	
		distinct id
		local ndist = r(ndistinct)
		local nn    = r(N)
		file write myfile %9s "Starting point: number of individuals: " %7.0f (`ndist') %9s ", number of spell-year-observations: " %7.0f (`nn') _n	
	
		*Dropping spells where the worker is either too young(<18) or too old (>60)
		gen age=year-anno_nascita
		keep if age>=18 & age<=60
		distinct id
		local ndist = r(ndistinct)
		local nn    = r(N)
		file write myfile %9s "After dropping spells with too young or too older workers: number of individuals: " %7.0f (`ndist') %9s ", number of spell-year-observations: " %7.0f (`nn') _n	
	
		*Drop spells with crazy amount of days worked. Flag these workers
		gen flag_G = 0 
		replace flag_G=1 if giorni_retribuiti >=500
		replace giorni_retribuiti=365 if giorni_retribuiti>=365
		bys id: egen flag_G0=mean(flag_G)
		drop if flag_G==1
		distinct id
		local ndist = r(ndistinct)
		local nn    = r(N)
		file write myfile %9s "After dropping spells with crazy amount of days worked: number of individuals: " %7.0f (`ndist') %9s ", number of spell-year-observations: " %7.0f (`nn') _n	
		drop flag_G
		
		*Drop spells with 0 Earnings, 0 spells. Flag these workers.
		gen flag_0 = 0 
		replace flag_0=1 if imponibile == 0 & giorni_retribuiti == 0 & settimane_retribuite == 0 | imponibile < 0 & giorni_retribuiti > 0 & settimane_retribuite > 0 & giorni_retribuiti !=. & settimane_retribuite != . | giorni_retribuiti == 0 & imponibile > 0
		sum flag_0, //how many spells report a zero.
		bys id: egen flag_w0=mean(flag_0)
		drop if flag_0==1
		distinct id
		local ndist = r(ndistinct)
		local nn    = r(N)
		file write myfile %9s "After dropping spells with 0 earnings or days worked or weeks worked: number of individuals: " %7.0f (`ndist') %9s ", number of spell-year-observations: " %7.0f (`nn') _n	
		drop flag_0
		
		*Drop spells with incredibly small pay (<2 real daily euros)
		drop if daily_wage<2
		distinct id
		local ndist = r(ndistinct)
		local nn    = r(N)
		file write myfile %9s "After dropping spells with daily wage less than 2 euros : number of individuals: " %7.0f (`ndist') %9s ", number of spell-year-observations: " %7.0f (`nn') _n	
	
		*Drop workers that had at some point a crazy amount of days worked
		drop if flag_G0>0
		distinct id
		local ndist = r(ndistinct)
		local nn    = r(N)
		file write myfile %9s "After dropping individuals with crazy amount of days worked: number of individuals: " %7.0f (`ndist') %9s ", number of spell-year-observations: " %7.0f (`nn') _n	
		drop flag_G0
		
		*Drop workers that at some point reported a crazy amount of jobs in a year.	
		bys id year: gen total_jobs_in_year=_N
		tab total_jobs_in_year 
		gen too_many_jobs=0
		replace too_many_jobs=1 if total_jobs_in_year>=10
		bys id: egen too_many=mean(too_many_jobs)
		drop if too_many>0
		drop too_many_jobs too_many	
		distinct id
		local ndist = r(ndistinct)
		local nn    = r(N)
		file write myfile %9s "After dropping workers with too many jobs: number of individuals: " %7.0f (`ndist') %9s ", number of person-year-observations: " %7.0f (`nn') _n	
		
		*Drop individuals that entered in the INPS data for the first time when they were very old or way too young.
		gen real_age_at_entry=anno_inizio_lavoro-anno_nascita
		sum real_age_at_entry
		drop if real_age_at_entry>=55
		drop if real_age_at_entry<=14
		distinct id
		local ndist = r(ndistinct)
		local nn    = r(N)
		drop anno_inizio_lavoro
		file write myfile %9s "After dropping individuals that entered way too soon or way too late: " %7.0f (`ndist') %9s ", number of person-year-observations: " %7.0f (`nn') _n	

		*Drop individuals that for some reason did not have their gender reported
		gen no_gender=0
		replace no_gender=1 if sesso!="F" & sesso!="M" 
		bys id: egen no_gender_aux=mean(no_gender)
		drop if no_gender_aux>0
		drop no_gender_aux no_gender
		distinct id
		local ndist = r(ndistinct)
		local nn    = r(N)
		file write myfile %9s "After dropping individuals with no gender reported: " %7.0f (`ndist') %9s ", number of person-year-observations: " %7.0f (`nn') _n	
		
		*Drop individuals that reported negative or zero earnings
		gen neg_earnings=0
		replace neg_earnings=1 if imponibile <= 0
		bys id: egen neg_earnings_aux=mean(neg_earnings)
		drop if neg_earnings_aux>0
		drop neg_earnings_aux neg_earnings
		distinct id
		local ndist = r(ndistinct)
		local nn    = r(N)
		file write myfile %9s "After dropping individuals with negative earnings: " %7.0f (`ndist') %9s ", number of person-year-observations: " %7.0f (`nn') _n	
*-------*END OF RESTRICTIONS


*SAVE THIS CLEANED SPELL DATA
		save "$dirVP/RAKIM/CLEANED_SPELL_GIANT_`type_data'", replace 


**************************************************************************

*	SECTION 4: SPLIT THE DATA TO RUN WITH PARALLEL
	
**************************************************************************	
*LOAD CLEANED SPELL DATA
		use id year using "$dirVP/RAKIM/CLEANED_SPELL_GIANT_`type_data'", replace

*QUANTILES OF THE DATA
		collapse year, by(id)
		drop year
		gen   conteggio = _n
		xtile qqq = conteggio, nquantiles(100)
		drop conteggio
		save "$dirVP/RAKIM/AUX",replace

*LOAD CLEANED SPELL DATA
		use "$dirVP/RAKIM/CLEANED_SPELL_GIANT_`type_data'", replace	
		merge m:1 id using 	"$dirVP/RAKIM/AUX", nogen
		
*SAVE THIS AUXILIARY FILE
		save "$dirVP/RAKIM/Q_CLEANED_SPELL_GIANT_`type_data'", replace
}		
**************************************************************************

*	SECTION 5: CREATE THE MONTHLY PANEL (PARALLEL CODING)
	
**************************************************************************	
if 1 == 1 {
*LOAD CLEANED SPELL DATA
		use if qqq == `ss' using "$dirVP/RAKIM/Q_CLEANED_SPELL_GIANT_`type_data'", replace 

*CLEAN MONTHLY DECLARATION VARIABLES
		forval tt=1(1)12{
		tab decl`tt', missing
		}
		
		forval tt=1(1)12{
		replace decl`tt'=0 if decl`tt'==3 | decl`tt'==8 | decl`tt'==9 | decl`tt'== 7 
		replace decl`tt'=1 if decl`tt'==2 | decl`tt'==5
		}
		
		forval tt=1(1)12{
		tab decl`tt', missing
		}

*SIMPLE TRICK
		gen double id_monthly = _n	
	    drop tipo_politica flag_w0 total_jobs_in_year anno_pensionamento mese_pensionamento livello_inquadramento ID_PERSF codcatn mese_morte anno_morte anno_nascita mt_qdrb_tipoc CF_PERSG flag_diversa_ settimane_ret settimane_utili
	    save "$dirVP/RAKIM/can_be_deleted`ss'", replace 
	    
*NOW ONLY IDS and DECLARATIONS
		keep id_monthly decl*
		reshape long decl, i(id_monthly) j(month)
		merge m:1 id_monthly using "$dirVP/RAKIM/can_be_deleted`ss'",  nogenerate
		drop decl1 decl2 decl3 decl4 decl5 decl6 decl7 decl8 decl9 decl10 decl11 decl12	
		capture rm 		"$dirVP/RAKIM/can_be_deleted`ss'"
		capture erase 	"$dirVP/RAKIM/can_be_deleted`ss'"

*CHECK SPELLS THAT APPEARED AS DUPLICATES
		bys id id_azienda year month decl: gen TTT=_N
		sum TTT,d 
		gsort id id_azienda year month -imponibile
		by id id_azienda year month: gen top_spell=_n
		keep if top_spell == 1
		drop top_spell TTT
			
*GEN INFO AT THE MONTHLY LEVEL		
		rename decl working
		bys id_monthly: egen months_worked = total(working)
		gen earnings_monthly = working*(imponibile/months_worked)
		gen days_monthly = working*(giorni_retribuiti/months_worked)
		
*Drop spells for those that never worked in a year. 
		sum months_worked,d
		gen      zero_months_worked_spell = 0 
		replace  zero_months_worked_spell = 1 if months_worked == 0 
		drop if zero_months_worked_spell  == 1
		drop zero_months_worked_spell
		
*TOTALS IN A GIVEN ID - YEAR - MONTH CELLS		
		bys id year month: egen total_days_month 	 = total(days_monthly)
		bys id year month: egen total_earnings_month = total(earnings_monthly)
		drop imponibile giorni_retribuiti
		
*COLLAPSE TO TOP JOB IN A GIVEN ID-YEAR-MONTH CELL
		gsort id year month -earnings_monthly
		by id year month: gen top_job=_n
		keep if top_job == 1
		drop top_job

**************************************************************************

*	SECTION 6: FILL GAPS IN THE ID-YEAR-MONTH PANEL USING TSFILL
	
**************************************************************************	
*BALANCE THE PANEL (TSFILL)
		gen date=ym(year,month)	
		format date %tm
		drop year month
		xtset id date
		tsfill
		list id date id_azienda in 1/100
		bysort id: carryforward sesso, replace
		bysort id: carryforward cittadinanza, replace
		bysort id: carryforward birth,replace
		bysort id: carryforward real_age_at_entry,replace
		drop age
		gen year  = yofd(dofm(date))
		gen month = mofd(dofm(date))
		gen age = year - birth

*THE IMPORTANT STUFF
		gen out_of_sample = 0
		replace out_of_sample        = 1 if working == .
		replace working              = 0 if out_of_sample == 1
		replace total_days_month     = 0 if out_of_sample == 1
		replace total_earnings_month = 0 if out_of_sample == 1
		replace earnings_monthly     = 0 if out_of_sample == 1
		replace days_monthly         = 0 if out_of_sample == 1
		replace id_impresa 			 = "unemployed" if working == 0 
		replace id_azienda 			 = -9 if working == 0 
		
*SAVE THE GIANT MONTHLY PANEL -- FUTURE REFERENCE FOR OTHER PROJECTS
		distinct id
		save "$dirVP/RAKIM/CLEANED_MONTHLY_PANEL_GIANT_`ss'_`type_data'", replace 
}		
**************************************************************************

*	SECTION 7: NOW IMPOSE THE RESTRICTIONS IN RAKIM (parallel coding)
	
**************************************************************************
if 1 == 1 {	
*LOAD THE MONTHLY PANEL
		use "$dirVP/RAKIM/CLEANED_MONTHLY_PANEL_GIANT_`ss'_`type_data'", replace 
			
*CREATE THE TENURE VARIABLE	
		xtset id date
		bys id id_impresa working (date): gen tenure =_n
		replace tenure = 0 if working == 0			

*FOR RAKIM I ONLY NEED START AND END OF A TENURE
		bys id id_impresa: egen maxTenure = max(tenure)
		keep if tenure == 1 | tenure == maxTenure & tenure> 0 // notice that I am not doing much about cycles of employment: started at firm A then moved to firm B then back at Firm A...
		xtset id date							
		  
*IT'S IMPORTANT TO UNDERSTAND THE TENURE VARIABLE
		sum maxTenure,d
		
********RESTRICTION: AT LEAST 3 MONTHS TO BE A PROPER JOB
		distinct id
		keep if maxTenure>=3
		distinct id
			    
*NOW UNDERSTAND WHETHER IT WAS AN HIRE FROM UNEMPLOYMENT OR AN HIRE FROM EMPLOYMENT
		xtset id date
		bys id: gen gap = date[_n] -  date[_n-1]
		replace gap = . if tenure !=1
		sum gap,d
			    
*BRING FORWARD THE INFORMATION ON WHY THE JOB WAS DESTROYED (SO NOW IT SHOWS WHEN TENURE==1). Notice that for jobs that lasted for less than a year you will see the information on separation both for rows when tenure=1 and when tenure=MaxTenure, see here: https://www.dropbox.com/s/aqsw4wz1wou341d/screen_RAKIM.png?dl=1 
		xtset id date
		bys id: gen MOTIVO = motivo_cessazione[_n-1] 
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
		
*KEEP ONLY POACHING WAGES (WE HAVE CALCULATED THE GAP)		
		keep if tenure == 1
				 	
*FINAL TOUCHES
		gen 	   female 			 = 0
		replace    female 			 = 1 		 	 if sesso == "F"
		gen		   firmsize			 = log(DIP10)

*SAVE				
		save "$dirVP/RAKIM/RAKIM_CLEANED_MONTHLY_PANEL_GIANT_`ss'_`type_data'", replace 
}		
**************************************************************************

*	SECTION 8: APPEND THE FILES
	
**************************************************************************
if 1 == 1 {	
*LOAD THE FIRST ONE
		use "$dirVP/RAKIM/RAKIM_CLEANED_MONTHLY_PANEL_GIANT_1_`type_data'", replace
		
		forval ss = 2(1)`total_files_app' {
		append using "$dirVP/RAKIM/RAKIM_CLEANED_MONTHLY_PANEL_GIANT_`ss'_`type_data'",
		}
		
*NOW GENERATE THE AUXILIARY LAGFIRMID
		capture drop firmid lagfirmid conteggio right_cen
		tostring id_azienda, replace
		replace id_impresa = id_azienda if id_impresa ==""
		egen firmid 	= group(id_impresa)
		bys id: gen lagfirmid = firmid[_n-1]
		bys id:  gen conteggio = _n
		replace lagfirmid 	  =	-1 if conteggio==1
		
*FLAG RIGHT CENSORED OBS
		gen 		right_cen 		 = 0
		replace		right_cen		 = 1			 if conteggio == 1 & year == 1990
		
*READY,SET, GO
		replace     lagfirmid 	 	 = 0 			 if gap >= 2 & gap!=. | reason == 1 | reason == 4 // hired from unemployment or following a layoff		


*SAVE THE DATASET
		save "$dirVP/RAKIM/RAKIM_APPENDED`total_files_app'_`type_data'", replace 
}	
**************************************************************************

*	DONE
	
**************************************************************************									 	    
