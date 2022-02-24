global it				= `1' //read local from bash submission
********************************************************************************

*		Figures and Tables from Di Addario, Kline, Saggio and Sølvsten 

********************************************************************************
set more off
global DIR_TMP_AUX "../build/src/RAKIM/"	
set scheme plotplain
local J_J_STRINGA = "1990_2016_2005_2016_RS_build"	
********************************************************************************

*FIGURE 1							

********************************************************************************
if 1 == 1 {
	capture log close
	log using "logs/figure1.log", replace
	do "codes/figure1"
	capture log close
}
********************************************************************************

*FIGURE 2							

********************************************************************************
if 1 == 1 {
	capture log close
	log using "logs/figure2.log", replace
	do "codes/figure2"
	capture log close
}
********************************************************************************

*FIGURE 3 AND TABLE 7 							

********************************************************************************
if 1 == 1 {
	capture log close
	log using "logs/figure3.log", replace
	do "codes/figure3"
	capture log close
}
********************************************************************************

*FIGURE 4							

********************************************************************************
if 1 == 1 {
	capture log close
	log using "logs/figure4.log", replace
	do "codes/figure4"
	capture log close
}
********************************************************************************

*FIGURE 5							

********************************************************************************
if 1 == 1 {
	capture log close
	global file_type = "POOLED_SAMPLE__INVIND_data_`J_J_STRINGA'"
	global JJSTRINGA = "`J_J_STRINGA'"
	log using "logs/figure5.log", replace
	do "codes/figure5"
	capture log close
}
********************************************************************************

*BUILD FOR OUT-OF-SAMPLE FIT OF DWL MODEL DONE IN FIGURE 6

********************************************************************************
if 1 == 1 {
	capture log close
	log using "logs/build_figure6.log", replace
	do "codes/build_figure6"
	capture log close
}
********************************************************************************

*FIGURE 6							

********************************************************************************
if 1 == 1 {
	capture log close
	log using "logs/figure6.log", replace
	do "codes/figure6"
	capture log close
}
********************************************************************************

*BUILD FOR OAXACA DECOMPOSITION

********************************************************************************
if 1 == 1 {
	capture log close
	global JJSTRINGA = "`J_J_STRINGA'"
	log using "logs/build_oaxaca`J_J_STRINGA'.log", replace
	do "codes/build_oaxaca"
	capture log close
}	
********************************************************************************

*FIGURE 7 AND FIGURE A3				

********************************************************************************
if 1 == 1 {
	capture log close
	log using "logs/figure7.log", replace
	do "codes/figure7"
	capture log close
}
********************************************************************************

*FIGURE 8							

********************************************************************************
if 1 == 1 {
	capture log close
	log using "logs/figure8.log", replace
	do "codes/figure8"
	capture log close
}
********************************************************************************

*FIGURE A1							

********************************************************************************
if 1 == 1 {
	capture log close
	global JJSTRINGA = "`J_J_STRINGA'"
	log using "logs/figureA1.log", replace
	do "codes/figureA1"
	capture log close
}
********************************************************************************

*FIGURE A2						

********************************************************************************
if 1 == 1 {
	capture log close
	global file_type = "POOLED_SAMPLE__INVIND_data_`J_J_STRINGA'"
	global JJSTRINGA = "`J_J_STRINGA'"
	log using "logs/figureA2.log", replace
	do "codes/figureA2"
	capture log close
}
********************************************************************************

*FIGURE A4					

********************************************************************************
if 1 == 1 {
	capture log close
	global file_type = "POOLED_SAMPLE__INVIND_data_`J_J_STRINGA'"
	global JJSTRINGA = "`J_J_STRINGA'"
	log using "logs/figureA4.log", replace
	do "codes/figureA4"
	capture log close
}
********************************************************************************

*TABLE A1					

********************************************************************************
if 1 == 1 {
	capture log close
	log using "logs/tableA1.log", replace
	do "codes/tableA1"
	capture log close
}
********************************************************************************

*TABLE A2					

********************************************************************************
if 1 == 1 {
	capture log close
	global file_type = "POOLED_SAMPLE__INVIND_data_`J_J_STRINGA'"
	global JJSTRINGA = "`J_J_STRINGA'"
	log using "logs/tableA2.log", replace
	do "codes/tableA2"
	capture log close
}
********************************************************************************

*CALCULATIONS OF PENALTY FROM BEING HIRED FROM UNEMPLOYMENT (ADDENDUM TO TABLE 5)					

********************************************************************************
if 1 == 1 {
	capture log close
	log using "logs/penalty.log", replace
	do "codes/penalty"
	capture log close
}

