local J_J_STRINGA = "$JJSTRINGA"
local dir_tmp_aux = "/scratch/public/leave_out"	
local file_type   = "$file_type"
*********************************************************************************

*			ORGANIZE THE FIRM EFFECTS FOUND IN THE B-N DAKM SPECIFICATION

*********************************************************************************
if 0 == 1 {
*import the file with the estimates from the DAKM model
	import delimited "`dir_tmp_aux'/POOLED_SAMPLE__INVIND_data_`J_J_STRINGA'.csv", encoding(ISO-8859-1) clear
	rename v1 y
	rename v2 id
	rename v3 firmid
	rename v4 lagfirmid
	rename v5 pe
	rename v6 psi
	rename v7 lambda
	rename v8 gamma
	save "`dir_tmp_aux'/`file_type'",replace
	
*create the firm level file now
	preserve
	collapse psi, by(firmid)
	sum psi
	save "`dir_tmp_aux'/psi_fe_BN",replace	
	restore
	collapse lambda, by(lagfirmid)
	save "`dir_tmp_aux'/lambda_fe_BN",replace

*load the original starting file
	use "`dir_tmp_aux'/B_N_BUILD_RAKIM_APPENDED100_INVIND_data_`J_J_STRINGA'",clear
	
*attach psi
	merge m:1 firmid    using "`dir_tmp_aux'/psi_fe_BN", gen(merge_psi)

*attach lambda
	merge m:1 lagfirmid using "`dir_tmp_aux'/lambda_fe_BN", gen(merge_lambda)
	
*just identified person-job observations
	keep if merge_psi == 3 & merge_lambda == 3

*collapse and reshape now
	preserve
	keep if old_lagfirmid>0
	collapse lambda, by(year)
	save "`dir_tmp_aux'/LAMBDA_no_state`file_type'",replace	 	
	restore
	preserve
	keep if old_lagfirmid== 0
	collapse lambda, by(year)
	rename lambda lambda_U 
	save "`dir_tmp_aux'/LAMBDA_U`file_type'",replace	 	
	restore 
	keep if old_lagfirmid== -1
	collapse lambda, by(year)
	rename lambda lambda_N
	merge 1:1 year using "`dir_tmp_aux'/LAMBDA_no_state`file_type'",nogen
	merge 1:1 year using "`dir_tmp_aux'/LAMBDA_U`file_type'",nogen

*merge unemployment rate in Italy
	merge 1:1 year using "`dir_tmp_aux'/u_rate", nogen

*save this auxiliary dataset
	save "`dir_tmp_aux'/BN_SAMPLE`file_type'",replace
}	

*load this auxiliary dataset
	use "`dir_tmp_aux'/BN_SAMPLE`file_type'",replace
	
*impose normalization
	describe
	sum lambda_N if year == 2005
	replace lambda_N = lambda_N-r(mean)
	replace lambda_U = lambda_U-r(mean)
	replace lambda 	 = lambda - r(mean)
	
	reg lambda_N U_RATE
	local coeff_N = round(_b[U_RATE],0.001)
	predict l_N,xb
	
	reg lambda_U U_RATE
	local coeff_U = round(_b[U_RATE],0.001)
	predict l_U,xb
	
	reg lambda U_RATE
	local coeff_E = round(_b[U_RATE],0.001)
	predict l_J,xb
	
	sum lambda
	local mybar = r(mean)
	
		
*do the line graph
	twoway (connected lambda_N year, yaxis(1) msymbol(square) lcolor(ltblue) mcolor(ltblue) legend(label(1 "{&lambda}{subscript:N} - State: Not yet in the Labor Market"))) (connected lambda_U year, yaxis(1) lcolor(gold) msymbol(square) mcolor(gold) legend(label(2 "{&lambda}{subscript:U} - State: Unemployment"))) (connected U_RATE year, yaxis(2) lcolor(orange) msymbol(triangle) mcolor(orange) legend(label(3 "Unemployment Rate"))), ytitle("{&lambda} --- State Effects", axis(1)) ytitle("Unemployment Rate in Italy", axis(2)) xtitle("Year") legend(order(1 2 3)) legend(ring(0) position(11) rows(3))  plotregion(margin(zero)) ylabel(-0.02(0.02)0.10)
	graph export "figures/BEAUDRY_DI_NARDO.pdf",replace
	graph save  "figures/BEAUDRY_DI_NARDO",replace
	
*do the line graph
	twoway (connected lambda_N year, yaxis(1) msymbol(square) lcolor(white) mcolor(white) legend(label(1 "{&lambda}{subscript:N} - State: Not yet in the Labor Market"))) (connected lambda_U year, yaxis(1) lcolor(white) msymbol(square) mcolor(white) legend(label(2 "{&lambda}{subscript:U} - State: Unemployment"))) (connected U_RATE year, yaxis(2) lcolor(orange) msymbol(triangle) mcolor(orange) legend(label(3 "Unemployment Rate"))), ytitle("{&lambda} --- State Effects", axis(1)) ytitle("Unemployment Rate in Italy", axis(2)) xtitle("Year") legend(off)  plotregion(margin(zero)) ylabel(-0.02(0.02)0.10)
	graph export "figures/BEAUDRY_DI_NARDO_1.pdf",replace
	graph save  "figures/BEAUDRY_DI_NARDO_1",replace	
		
*do the scatter plot
	twoway (scatter lambda_N U_RATE, yaxis(1) msymbol(square) mcolor(ltblue) legend(label(1 "{&lambda}{subscript:N}"))) (scatter lambda_U U_RATE, yaxis(1) msymbol(square) mcolor(gold) legend(label(2 "{&lambda}{subscript:U}"))) (lfit lambda_N U_RATE, yaxis(1) lcolor(ltblue) lpattern(solid) legend(label(3 "{&lambda}{subscript:N}"))) (lfit lambda_U U_RATE, yaxis(1) lpattern(solid) lcolor(gold) legend(label(4 "{&lambda}{subscript:U}"))), ytitle("{&lambda} --- State Effects", axis(1)) xtitle("Unemployment Rate") legend(order(1 2)) legend(ring(2) position(6) rows(1)) note("Note {&lambda} normalized w.r.t. {&lambda}{subscript:N} in 2005" "Red Line correspond to average {&lambda} excluding N and U lagged employers")	yline(`mybar', lcolor(red))
	graph export "figures/BEAUDRY_DI_NARDO_scatter.pdf",replace
	graph save  "figures/BEAUDRY_DI_NARDO_scatter",replace
	
*do the scatter plot with lambda employers 
	twoway (scatter lambda_N U_RATE, yaxis(1) msymbol(square) mcolor(ltblue) legend(label(1 "{&lambda}{subscript:N} - State: Not yet in the Labor Market"))) (scatter lambda_U U_RATE, yaxis(1) msymbol(square) mcolor(gold) legend(label(2 "{&lambda}{subscript:U}- State: Non-Employed "))) (line l_N U_RATE, yaxis(1) lcolor(ltblue) lpattern(solid) legend(label(3 "{&lambda}{subscript:N}"))) (line l_U U_RATE, yaxis(1) lpattern(solid) lcolor(gold) legend(label(4 "{&lambda}{subscript:U} "))) (scatter lambda U_RATE, yaxis(1) mcolor(orange) msymbol(triangle) legend(label(5 "{&lambda}{subscript:J} - State: Job-Job Transition"))) (line l_J U_RATE, yaxis(1) lpattern(solid) lcolor(orange) legend(label(6 "{&lambda}{subscript:J}"))), ytitle("{&lambda} --- State Effects", axis(1)) xtitle("Unemployment Rate") xlabel(6(2)13) legend(order(1 2 5)) legend(ring(0) position(7) rows(3))  text(0.09 12 "Slope: `coeff_E'", size(small) color(orange)) text(-0.01 12.3 "Slope: `coeff_N'", size(small) color(ltblue)) text(0.0075 12.3 "Slope: `coeff_U'", size(small) color(gold))
	graph export "figures/BEAUDRY_DI_NARDO_scatter_aug.pdf",replace
	graph save  "figures/BEAUDRY_DI_NARDO_scatter_aug",replace
		