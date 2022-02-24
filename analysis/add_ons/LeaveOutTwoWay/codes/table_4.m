function TABELLA=table_4(y,id,firmid,leave_out_level,leave_2_out,controls,type_of_algorithm,eigen_diagno,eigen_fast,do_montecarlo,restrict_movers,epsilon,filename);
	
	%Run KSS
	[one two] = leave_out_FD(y,id,firmid,leave_out_level,leave_2_out,controls,type_of_algorithm,eigen_diagno,eigen_fast,do_montecarlo,restrict_movers,epsilon,1,filename);
	
	%Run KSS
	[one two] = leave_out_FD(y,id,firmid,leave_out_level,leave_2_out,controls,type_of_algorithm,eigen_diagno,eigen_fast,do_montecarlo,restrict_movers,epsilon,2,filename);
		
	TABELLA=0;	
	
end


