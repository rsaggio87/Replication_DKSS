local dir_tmp_aux = "/scratch/public/leave_out"	

*call file created by matlab
import delimited "`dir_tmp_aux'/symmetry.csv", encoding(ISO-8859-1) clear
rename v1 y
rename v4 firmid_orig
rename v5 firmid_destin
rename v6 change_AKM
rename v7 symmetry
rename v8 unrestricted

*merge VA information
gen firmid = firmid_orig
merge m:1 firmid using "`dir_tmp_aux'/va_file",gen(merge_origin)
keep if merge_origin ==3 
rename DIP10 DIP10_orig
rename log_VA_L log_VA_L_orig
drop firmid
gen firmid = firmid_destin
merge m:1 firmid using "`dir_tmp_aux'/va_file",gen(merge_destin)
keep if merge_destin ==3 
rename DIP10 DIP10_destin
rename log_VA_L log_VA_L_destin
drop firmid
gen change_VA  = log_VA_L_destin - log_VA_L_orig
replace DIP10_destin=log(DIP10_destin)
replace DIP10_orig=log(DIP10_orig)
gen change_SIZE = DIP10_destin-DIP10_orig
save "`dir_tmp_aux'/symmetry_build",replace

*collapse
use "`dir_tmp_aux'/symmetry_build",replace
xtile qqq=change_VA, nquantiles(100)
collapse change_AKM unrestricted,by(qqq)
reg unrestricted change_AKM
local sloppami 	= round(_b[change_AKM],0.01)
local constante = round(_b[_cons],0.01)
twoway (scatter unrestricted change_AKM, msymbol(square) mcolor(cranberry)) (lfit change_AKM change_AKM, lcolor(cranberry)), xtitle("Change in AKM effects") ytitle("OD Cell (Directed)") note("45 Degree Line Shown." "Each Bin corresponds to a centile of the change in Log VA/L b/w destination and origin firm" "Slope from cell data: `sloppami'; Constant: `constante'") legend(off)
graph export "figures/symmetry_OD.pdf",replace

*collapse (using firm size)
use "`dir_tmp_aux'/symmetry_build",replace
xtile qqq=change_SIZE, nquantiles(100)
collapse change_AKM unrestricted,by(qqq)
reg unrestricted change_AKM
local sloppami 	= round(_b[change_AKM],0.01)
local constante = round(_b[_cons],0.01)

twoway (scatter unrestricted change_AKM, msymbol(square) mcolor(cranberry)) (lfit change_AKM change_AKM, lcolor(cranberry)), xtitle("Change in AKM effects") ytitle("OD Cell (Directed)") note("45 Degree Line Shown" "Each Bin corresponds to a centile of the change in Log Firm Size b/w destination and origin firm" "Slope from cell data: `sloppami'; Constant: `constante'") legend(off)
graph export "figures/symmetry_OD_size.pdf",replace
