function new_firmid = normalize_firm_effects(firmid,to_norm)
	%Renormalize the set of firm effects according to the indicator specified by key_firm
	
	new_firmid			= firmid;
	firmid_saved		= firmid;
	sel			 		= firmid==max(firmid);
	new_firmid(sel)	 	= to_norm;
	
	sel			 		= firmid_saved==to_norm;
	new_firmid(sel)  	= max(firmid);
	
end