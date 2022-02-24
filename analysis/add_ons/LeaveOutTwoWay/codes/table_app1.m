function TABELLA=table_app1(y,id,firmid,leave_out_level,year,controls,resid_controls,andrews_estimates,eigen_diagno,subsample_llr_fit,restrict_movers,do_SE, type_of_algorithm,epsilon,filename)
	controls_orig=controls;
	
	
	%Run KSS
	[one] = leave_out_COMPLETE(y,id,firmid,leave_out_level,year,controls,resid_controls,andrews_estimates,eigen_diagno,subsample_llr_fit,restrict_movers,do_SE, type_of_algorithm,epsilon,filename);
	TABELLA(1)=one;
	
	%Run leaving a match out 
	leave_out_level = 'matches'
	[one] = leave_out_COMPLETE(y,id,firmid,leave_out_level,year,controls,resid_controls,andrews_estimates,eigen_diagno,subsample_llr_fit,restrict_movers,do_SE, type_of_algorithm,epsilon,filename);
	TABELLA(2)=one;
	
	%Run leaving a worker out
	leave_out_level = 'workers'	   	
	one = leave_out_FD(y,id,firmid,leave_out_level,0,controls,type_of_algorithm,eigen_diagno,0,0,1,epsilon,1,filename)
	TABELLA(3)=one;
	
end


