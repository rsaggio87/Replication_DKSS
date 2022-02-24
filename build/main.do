*******************************************

*DECLARATIONS

*******************************************	
global SPEC				= `1' // reads local specified from bash-submission
global it				= `2'
local array_id 			= `3'
set seed 1234
scalar drop _all
matrix drop _all
scalar drop _all
capture log close
set more off
global dirVP "src/" // path to main directory
global TYPE_DATA "INVIND_data"
cd $dirVP
set scheme plotplain
*******************************************

*BUILD

*******************************************	
if `array_id' == 1 {
	capture log close
	log using "logs/build_RAKIM.log", replace
	do "codes/build_RAKIM.do"
	capture log close
}

if `array_id' == 2 {
	capture log close
	log using "logs/export_matlab.log", replace
	do "codes/export_matlab.do"
	capture log close
}

if `array_id' == 3 {
	capture log close
	log using "logs/summary_build.log", replace
	do "codes/summary_build.do"
	capture log close
}

if `array_id' == 4 {
	capture log close
	log using "logs/build_AKM.log", replace
	do "codes/build_AKM.do"
	capture log close
}
