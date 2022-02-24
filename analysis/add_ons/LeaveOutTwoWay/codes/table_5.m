function TABELLA=table_5(y,id,firmid,leave_out_level,leave_2_out,controls,type_of_algorithm,eigen_diagno,eigen_fast,do_montecarlo,restrict_movers,epsilon,filename);

	%Run KSS
	[one two] = leave_out_FD(y,id,firmid,leave_out_level,leave_2_out,controls,type_of_algorithm,eigen_diagno,eigen_fast,do_montecarlo,restrict_movers,epsilon,1,filename);
	
	%Load the auxiliary file for table
	s=[filename 'MCCC_FD_NEW_V3'];
	load(s)
	
	%Create the table
	TABELLA(1) 		= sigma2_true;
	TABELLA(2) 		= mean(theta);
	TABELLA(3) 		= std(theta);
	TABELLA(4) 		= mean(theta_AKM);
	TABELLA(5) 		= std(theta_AKM);
	TABELLA(6) 		= mean(theta_homo);
	TABELLA(7) 		= std(theta_homo);
	TABELLA(8) 		= mean(V_SIM);
	TABELLA(9) 		= mean(coverage);
	TABELLA(10) 	= mean(mean(coverage_AM));
	
end


