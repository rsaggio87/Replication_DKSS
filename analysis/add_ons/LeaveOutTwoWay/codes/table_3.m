function TABELLA_FINALE=table_3(y,id,firmid,leave_out_level,year,controls,resid_controls,andrews_estimates,eigen_diagno,subsample_llr_fit,restrict_movers,do_SE, type_of_algorithm,epsilon,filename)

	%Run KSS
	placeFile='results/';
	filename=[placeFile 'RESULTS_ALL_VENETO_LEAVE1_OUT'];
	%[one] = leave_out_COMPLETE(y,id,firmid,leave_out_level,year,[year==2001],resid_controls,andrews_estimates,eigen_diagno,subsample_llr_fit,restrict_movers,do_SE, type_of_algorithm,epsilon,filename);
						 
	%Load results associated with Leave One Out Connected Set
	results=importdata('results/RESULTS_ALL_VENETO_LEAVE1_OUT.csv');
	y=results(:,1);
	firmid=results(:,2);
	id=results(:,3);
	old_dummy=controls;
	clear results

	%Load the matrix of statistical leverages saved by leave_out_COMPLETE
	Lambda_P=load('results/RESULTS_ALL_VENETO_LEAVE1_OUT_Lambda_P');
	Lambda_P=Lambda_P.Lambda_P;

	%Now create the X associated with the two-way model (dropping as usual last firm).
	NT=size(y,1);
	D=sparse(1:NT,id',1); 
	F=sparse(1:NT,firmid',1);
	GRANDE=F'*F;
	GRANDE=diag(GRANDE);
	median_firm_size=quantile(GRANDE,0.5)
	median_firm_size=quantile(log(GRANDE),0.5)
	GRANDE=log(F*GRANDE);
	N=size(D,2);
	J=size(F,2);
	S=speye(J-1);
	S=[S;sparse(-zeros(1,J-1))];
	F=F*S;
	X=[D,F]; %Design Matrix.
	
	%Test on FIRM EFFECTS - MODEL 1
	Transform=[sparse(NT,N) F]; %Transform*beta will therefore give back the firm effects in the person-year space.
	Z=[old_dummy];
	labels={'Old Dummy'};
	
	%Run Lincom KSS
	%In the code v is always defined as v=(Z'*Z)^(-1)*Z'*Transform and the
	%function "lincom_KSS" always adds a constant to the user specified matrix
	%"Z".

	%Therefore in this case the second element of v'*beta returns the regression
	%coefficient that captures the difference in firm effects between region 2
	%and region 1. 

	%We want to conduct robust inference on this particular
	%regression coefficient. The function lincom_KSS outputs the t-test
	%of this particular linear combination.


	%Run
	[~, TABELLA_FINALE(:,1), TABELLA_FINALE(:,2), TABELLA_FINALE(:,3),  TABELLA_FINALE(:,4) stat, pvalue]=lincom_KSS(y,X,Z,Transform,[],Lambda_P,labels); %Note: the [] implies that we will provide heteroskedatic robust inference (i.e. no clustering) which is consistent with the prior step where we conducted leave out estimation by leaving a single observation out each time (leave_out_level='obs';)
	TABELLA_FINALE=TABELLA_FINALE';	

	
	
	%Test on FIRM EFFECTS - MODEL 2
	Transform=[sparse(NT,N) F]; %Transform*beta will therefore give back the firm effects in the person-year space.
	Z=[old_dummy GRANDE old_dummy.*GRANDE]; 
	labels={'Old Dummy','Firm Size','Interaction'};
	
	%Run Lincom KSS
	%In the code v is always defined as v=(Z'*Z)^(-1)*Z'*Transform and the
	%function "lincom_KSS" always adds a constant to the user specified matrix
	%"Z".

	%Therefore in this case the second element of v'*beta returns the regression
	%coefficient that captures the difference in firm effects between region 2
	%and region 1. 

	%We want to conduct robust inference on this particular
	%regression coefficient. The function lincom_KSS outputs the t-test
	%of this particular linear combination.


	%Run
	[~, TABELLA(:,1), TABELLA(:,2), TABELLA(:,3),  TABELLA(:,4) stat, pvalue]=lincom_KSS(y,X,Z,Transform,[],Lambda_P,labels); %Note: the [] implies that we will provide heteroskedatic robust inference (i.e. no clustering) which is consistent with the prior step where we conducted leave out estimation by leaving a single observation out each time (leave_out_level='obs';)
	TABELLA_FINALE=[TABELLA_FINALE -9.*ones(4,1) -9.*ones(4,1) TABELLA'];	

end


