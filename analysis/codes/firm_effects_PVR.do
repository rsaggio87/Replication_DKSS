
*read the global
	local J_J_STRINGA = "$JJSTRINGA"
	local file_type = "$file_type"
	local dir_tmp_aux 	= "../build/src/RAKIM/"	

*import the file with the estimates from the DAKM model
	import delimited "`dir_tmp_aux'/`file_type'.csv", encoding(ISO-8859-1) clear
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
	save "`dir_tmp_aux'/fe_`file_type'",replace	
	restore
	collapse lambda, by(lagfirmid)
	rename lagfirmid firmid
	merge 1:1 firmid using "`dir_tmp_aux'/fe_`file_type'", gen(merge_dual)
	save "`dir_tmp_aux'/fe_`file_type'",replace

	
*attach firm effects from AKM 	
	import delimited "`dir_tmp_aux'/AKM`file_type'.csv", encoding(ISO-8859-1) clear
	rename v1 y
	rename v2 id
	rename v3 firmid
	rename v4 lagfirmid
	rename v5 psi_AKM
	
	collapse psi_AKM, by(firmid)
	merge 1:1 firmid using "`dir_tmp_aux'/fe_`file_type'", gen(merge_AKM) // merge_AKM == 2 are firms that have only an identified lambda, not a psi.

	save "`dir_tmp_aux'/fe_`file_type'",replace 
	
*bring information on firm size, VA, codice fiscale
	use DIP10 id_impresa firmid year using "`dir_tmp_aux'/SUBMITTED_RAKIM_APPENDED_100_INVIND_data_`J_J_STRINGA'",clear
	destring id_impresa, gen(codfis_mod) force
	tostring codfis_mod, format(%016.0f) replace
	collapse DIP10 (first) codfis_mod, by(firmid year) 
	merge m:1 codfis_mod year using "`dir_tmp_aux'/cerved_smaller",gen(merge_VA) 
	keep if merge_VA == 1 | merge_VA == 3
	
*clean value added per worker
	replace DIP10 = 1 if DIP10 == 0
	gen VA_L = valagg/DIP10	
	
*windsorize values at 5 and 95 in each year
	sum year 
	local bottom = r(min)
	local up = r(max)
	forval tt=`bottom'(1)`up'{
	sum VA_L if year == `tt',d
	replace 	VA_L = r(p5)  if VA_L<=r(p5)  & year == `tt'
	replace 	VA_L = r(p95) if VA_L>=r(p95) & VA_L!=. & year == `tt'
	}
	gen log_VA_L = log(VA_L)
	sum 	log_VA_L,d
	
*now collapse at the firm level
	collapse DIP10 log_VA_L, by(firmid) // this is going to average size and VA across years.
	save "`dir_tmp_aux'/va_file",replace
	merge 1:1  firmid using "`dir_tmp_aux'/fe_`file_type'", gen(merge_size) // firms that cannot be merged are firms that only appear as lagged firms and never as contemporaneous
	keep if merge_size == 3 | merge_size == 2
	
*normalize the lag firm effects now
	sum lambda if firmid == -1
	replace lambda = lambda - r(mean)
	
*normalize cont firm effects now (using firm size)
	xtile vinG = DIP10 if psi!=., nquantile(20)	
	sum psi [aw=DIP10] if vinG == 1
	scalar psi_bar = r(mean)
	gen psi_size	   = psi - psi_bar
	sum psi_AKM [aw=DIP10] if vinG == 1
	scalar psi_bar = r(mean)
	gen psi_AKM_size	 = psi_AKM - psi_bar
	drop vinG
	
*normalize cont firm effects now (using VA/L)
	xtile vinG = log_VA_L if psi!=., nquantile(20)	
	sum psi [aw=DIP10] if vinG == 1
	scalar psi_bar = r(mean)
	replace psi	   = psi - psi_bar
	sum psi_AKM [aw=DIP10] if vinG == 1
	scalar psi_bar = r(mean)
	replace psi_AKM	 = psi_AKM - psi_bar
	
*done
	save "`dir_tmp_aux'/fe_`file_type'",replace		
