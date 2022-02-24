
local dir_tmp_aux 		= "src/RAKIM" 
*********************************************************************************

*			THIS DO-FILE BULDS THE TABLES USED IN TABLE 4

*********************************************************************************
*LOAD FILE SAVED BY MATLAB
if 1 == 1{
	import delimited "`dir_tmp_aux'/POOLED_SAMPLE__INVIND_data_1990_2016_2005_2016_RS_build.csv", encoding(ISO-8859-1) clear
	rename v1 log_dailywages
	rename v2 id
	rename v3 firmid 
	rename v4 lagfirmid
	rename v6 fe
	rename v7 lag_fe
	rename v8 controls
	rename v19 age2 
	rename v20 age3
	gen year = 2004
	local aux = 1
	forval pp = 9(1)18 {
	replace year = 2004 + `aux' if v`pp' == 1
	local aux = `aux' +1
	}
	replace year = 2015 if year == 2004
	keep log_dailywages id firmid lagfirmid age2 age3 year
	bys id: gen contami = _n
	bys id: gen TOTALE  =_N 
	xtset id contami
	save "`dir_tmp_aux'/estimation_sample_pooled",replace	
	gen uno=1
	collapse uno, by(firmid lagfirmid age2 age3 log_dailywages year)
	drop uno
	save "`dir_tmp_aux'/estimation_sample_pooled_coll",replace
}
if 1 == 1{	
*now load the original src file with hiring wages (prior to pruning)
	use log_dailywages id_impresa id_azienda id firmid lagfirmid year age female firmsize date provincia_lavoro using "`dir_tmp_aux'/DAVE_BUILD_RAKIM_APPENDED100_INVIND_data_1990_2016_2005_2016_RS_build",replace
	gen provincia_recovered = provincia_lavoro
	merge m:1 provincia_recovered using "tables/PROVINCE",nogen //merge province cross-walk
	bys id: gen contami = _n
	gen age2= ((age-40)/40)^2
	gen age3= ((age-40)/40)^3
	
*focus on the estimation sample	
	merge m:1 firmid lagfirmid age2 age3 year log_dailywages  using "`dir_tmp_aux'/estimation_sample_pooled_coll", gen(merge_back)
	keep if merge_back == 3
	sum log_dailywages,d
	distinct id //slight over-coverage.
	distinct firmid
	distinct id_impresa

*now collapse to identified and pruned matches	
	collapse log_dailywages firmid female lagfirmid firmsize (first) descrizione_regione, by(id id_impresa)
	save "`dir_tmp_aux'/identified_matches",replace
	
*find these identified matches in the Giant Spell-level data	
	merge 1:m id id_impresa using "/scratch/public/board/RAKIM/Q_CLEANED_SPELL_GIANT_INVIND_data_1990_2016", gen(merge_spells) keepusing(imponibile year age)
	keep if year>=2005
	keep if merge_spells == 3
	gsort id year -imponibile
	by id year: gen top_job=_n
	
*retain only dominant jobs within the year but keep hiring + incumbent wages (SAMPLE IN COLUMN 3)
	keep if top_job == 1 // only dominant jobs in the year
	drop lagfirmid
	gen lagfirmid =1
	xtset id year
	distinct id
	distinct firmid 
	export delimited id firmid lagfirmid year log_dailywages age female using "`dir_tmp_aux'/INVIND_RAKIM_PY_AKM_POOLED_SAMPLE__INVIND_data_1990_2016_2005_2016_RS_build",replace novarnames nolabel
	distinct firmid
	save "`dir_tmp_aux'/INVIND_RAKIM_PY_AKM_POOLED_SAMPLE__INVIND_data_1990_2016_2005_2016_RS_build",replace
	
*retain only dominant jobs within the year with hiring wage (SAMPLE IN COLUMN 2)	
	bys id id_impresa (year): gen tenure =_n
	keep if tenure == 1 //just hiring wage frequencies
	distinct id
	distinct firmid
	distinct id_impresa
	xtset id year 		
	export delimited id firmid lagfirmid year log_dailywages age female using "`dir_tmp_aux'/INVIND_RAKIM_DOMINANT_AKM_POOLED_SAMPLE__INVIND_data_1990_2016_2005_2016_RS_build",replace novarnames nolabel
	
	
*call the sample in column 3,retain the firms there 
	use "`dir_tmp_aux'/INVIND_RAKIM_PY_AKM_POOLED_SAMPLE__INVIND_data_1990_2016_2005_2016_RS_build",clear
	
	collapse firmid, by(id_impresa)
	merge 1:m id_impresa using "/scratch/public/board/RAKIM/Q_CLEANED_SPELL_GIANT_INVIND_data_1990_2016", gen(merge_firms) keepusing(id log_dailywages imponibile year age sesso) //merge the original src spell data
	keep if year>=2005
	gsort id year -imponibile
	by id year: gen top_job=_n
	keep if top_job == 1 // keep only dominant jobs
	xtset id year
	bys id: gen lagfirmid = id_impresa[_n-1]
	bys id: gen contami   = _n

*find the stayers of those firms
	gen 		mover_aux = 0
	replace  	mover_aux = 1 if id_impresa !=lagfirmid & contami >1
	bys id:	    egen MOVER = mean(mover_aux)
	gen			always_firm_aux = 0
	replace		always_firm_aux = 1			if merge_firms == 3
	bys id: 	egen ALWAYS = mean(always) 
	keep 		if ALWAYS == 1 & MOVER == 0 // keep these firm stayers.
	
*now append
	drop lagfirmid
	gen lagfirmid =1
	gen female = 0
	replace female =1 if sesso == "F"
	keep id firmid lagfirmid year log_dailywages age female
	save "`dir_tmp_aux'/STAYERS",replace
}	
	use id firmid lagfirmid year log_dailywages age female using "`dir_tmp_aux'/INVIND_RAKIM_PY_AKM_POOLED_SAMPLE__INVIND_data_1990_2016_2005_2016_RS_build",clear
	distinct id
	distinct firmid
	sum log_dailywages,d
	xtset id year
	append using "`dir_tmp_aux'/STAYERS"
	distinct id
	distinct firmid
	sum log_dailywages,d
*xtset id year // if you were to run this line it would break because of the py collapse some dudes are now "stayers" yet they would live in the estimation sample. id=87 for instance went from firm 1866524 in 2014 to some other firm in 2014 but then back to 1866524 in 2015
	bys id year: gen contami = _n
	keep if contami == 1
	xtset id year 
	export delimited id firmid lagfirmid year log_dailywages age female using "`dir_tmp_aux'/INVIND_RAKIM_WITH_STAYERS_PY_AKM_POOLED_SAMPLE__INVIND_data_1990_2016_2005_2016_RS_build",replace novarnames nolabel

	
	
	
