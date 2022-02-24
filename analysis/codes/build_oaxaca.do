local J_J_STRINGA = "$JJSTRINGA"
local dir_tmp_aux 	= "../build/src/RAKIM/"	
*********************************************************************************

*			START BY FINDING AND ORGANIZING THE FIRM EFFECTS FOR FEMALES

*********************************************************************************
if 1 == 1 {
global file_type = "FEMALE_SAMPLE__INVIND_data_`J_J_STRINGA'"
do "codes/firm_effects_PVR"
*********************************************************************************

*			NOW DO THE SAME FOR MALES

*********************************************************************************
global file_type = "MALE_SAMPLE__INVIND_data_`J_J_STRINGA'"
do "codes/firm_effects_PVR"

*********************************************************************************

*			NOW COMBINE THE FILES JUST CREATED INTO ONE

*********************************************************************************
local file_type = "FEMALE_SAMPLE__INVIND_data_`J_J_STRINGA'"
use "`dir_tmp_aux'/fe_`file_type'",replace
rename psi psi_female
rename lambda lambda_female
rename psi_AKM psi_AKM_female
rename merge_dual merge_dual_female
rename merge_AKM merge_AKM_female
rename merge_size merge_size_female
drop DIP10
save "`dir_tmp_aux'/renamed_fe_`file_type'",replace
**
local file_type = "MALE_SAMPLE__INVIND_data_`J_J_STRINGA'"
use "`dir_tmp_aux'/fe_`file_type'",replace
rename psi psi_male
rename lambda lambda_male
rename psi_AKM psi_AKM_male
rename merge_dual merge_dual_male
rename merge_AKM merge_AKM_male
rename merge_size merge_size_male
**now merge
local file_type = "FEMALE_SAMPLE__INVIND_data_`J_J_STRINGA'"
merge 1:1 firmid using "`dir_tmp_aux'/renamed_fe_`file_type'",gen(dual_gender_merge)
**keep only dual connected firms 
keep if dual_gender_merge == 3
drop dual_gender_merge
*now save
local file_type = "COMBINED_SAMPLE_INVIND_data_`J_J_STRINGA'"
save "`dir_tmp_aux'/combined_fe_`file_type'",replace
save "`dir_tmp_aux'/firm_level_combined_fe_`file_type'",replace
*********************************************************************************

*			NOW CREATE TWO SEPARATE FILES: ONE FOR LAMBDA AND ONE FOR PSI

*********************************************************************************
*contemporary effects
use "`dir_tmp_aux'/combined_fe_`file_type'",replace
keep firmid psi*
keep if psi_male !=. & psi_female !=.
save "`dir_tmp_aux'/combined_psi_fe_`file_type'",replace
*lagged effects
use "`dir_tmp_aux'/combined_fe_`file_type'",replace
keep firmid lambda* 
rename firmid lagfirmid
keep if lambda_male !=. & lambda_female !=.
save "`dir_tmp_aux'/combined_lambda_fe_`file_type'",replace
*********************************************************************************

*			NOW MERGE THIS FILE INTO THE STARTING SAMPLE -- THIS COMPLETES THE BUILD FOR THE OAXACA ANALYSIS

*********************************************************************************
}
use "`dir_tmp_aux'/DAVE_BUILD_RAKIM_APPENDED100_INVIND_data_1990_2016_2005_2016_RS_build",replace
merge m:1 firmid 	using "`dir_tmp_aux'/combined_psi_fe_`file_type'", gen(dual_mee_psi)
merge m:1 lagfirmid using "`dir_tmp_aux'/combined_lambda_fe_`file_type'", gen(dual_mee_lambda)
keep if dual_mee_psi == 3 & dual_mee_lambda == 3 // identified person-job observations for which we can find a (psi,lambda) both in the femele and male sample
if 1 == 1 {
preserve
import delimited "`dir_tmp_aux'/FEMALE_SAMPLE__INVIND_data_`J_J_STRINGA'.csv", encoding(ISO-8859-1) clear
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
sum year
collapse controls, by(age2 age3 year)
rename controls controls_female
gen female = 1
save "`dir_tmp_aux'/aux_file_females",replace
restore
gen age2= ((age-40)/40)^2
gen age3= ((age-40)/40)^3
merge m:1 female age2 age3 year using "`dir_tmp_aux'/aux_file_females", gen(merge_female) 
preserve
import delimited "`dir_tmp_aux'/MALE_SAMPLE__INVIND_data_`J_J_STRINGA'.csv", encoding(ISO-8859-1) clear
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
sum year
collapse controls, by(age2 age3 year)
rename controls controls_male
gen female = 0
save "`dir_tmp_aux'/aux_file_males",replace
restore
merge m:1 female age2 age3 year using "`dir_tmp_aux'/aux_file_males", gen(merge_male)    keepusing(controls)
}
gen res		 	= log_dailywages-psi_male-lambda_male-controls_male       	if female == 0
replace res  	= log_dailywages-psi_female-lambda_female-controls_female 	if female == 1
bys id: egen pe = mean(res)
gen prediction = psi_male+lambda_male+controls_male 						if female == 0
replace prediction = psi_female+lambda_female+controls_female 				if female == 1
save "`dir_tmp_aux'/oaxaca_sample_`file_type'",replace
