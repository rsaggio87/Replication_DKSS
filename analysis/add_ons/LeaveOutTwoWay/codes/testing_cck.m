function [test V student]= testing_cck(key_firm,r,filename,results_0,results_1,leave_out_level,type_of_algorithm,epsilon)
	
	%Start by loading results for leave out connected set of workers in group 0
	results						= results_0;
	y_group0					= results(:,1);
	firmid_group0				= results(:,2);
	id_group0					= results(:,3);
	group0						= results(:,7);
	firmid_orig_group0			= results(:,4);
	id_orig_group0				= results(:,5);
	clear results
	
	%Start by loading results for leave out connected set of workers in group 1
	results                     = results_1;
	y_group1					= results(:,1);
	firmid_group1				= results(:,2);
	id_group1					= results(:,3);
	group1						= results(:,7);
	firmid_orig_group1			= results(:,4);
	id_orig_group1				= results(:,5);
	clear results results_0 results_1

	%Now stack entirely the two datasets
	id_stack_orig				= [id_orig_group0; id_orig_group1];
	id_stack					= [id_group0; max(id_group0)+id_group1];
	firmid_orig_stack			= [firmid_orig_group0; firmid_orig_group1];
	[~,~,firmid_stack]			= unique(firmid_orig_stack);
	y_stack						= [y_group0; y_group1];
	group_stack					= [group0; group1];
	
	%Normalize the firm effects
	to_norm						= mean(firmid_stack(firmid_orig_stack==key_firm))
	firmid_stack				= normalize_firm_effects(firmid_stack,to_norm);
	
	%Begin by running the restricted model in first differences, save leverages and predictors
	strinGA						= [filename '_RESTRICTED_MODEL']
	sigma2_psi  				= leave_out_FD(y_stack,id_stack,firmid_stack,leave_out_level,0,group_stack,type_of_algorithm,0,0,0,0,epsilon,1,strinGA);
	load([strinGA 'FD_completed'],'ydelta','firmid_delta','firmid_delta_f','id_movers','Lambda_P','ydelta_pred','controls_delta')
	y_pred_restrict				= ydelta_pred;
	Pii_restrict				= diag(Lambda_P);
	ydelta_stack				= ydelta;
	Ndelta						= size(ydelta,1);
	index_C						=(1:Ndelta)';	
	GROUP_delta					= controls_delta;
	firmid_delta				= firmid_delta;
	firmid_delta_f				= firmid_delta_f;
	
	%Run model on group 0, save leverages, predictors and KSS sigmas
	strinGA						= [filename '_UNRESTRICTED_MODEL_GROUP0'];
	sel							= group_stack==0;
	sigma2_psi  				= leave_out_FD(y_stack(sel),id_stack(sel),firmid_stack(sel),leave_out_level,0,[],type_of_algorithm,0,0,0,0,epsilon,1,strinGA);
	load([strinGA 'FD_completed'],'Lambda_P','ydelta_pred','eta_h','ydelta')
	y_pred_unrestrict			= ydelta_pred;
	Pii_unrestrict				= diag(Lambda_P);
	sigma_i						= ydelta.*eta_h;

	%Run model on group 1, save leverages, predictors and KSS sigmas
	strinGA						= [filename '_UNRESTRICTED_MODEL_GROUP1'];
	sel							= group_stack==1;
	sigma2_psi  				= leave_out_FD(y_stack(sel),id_stack(sel),firmid_stack(sel),leave_out_level,0,[],type_of_algorithm,0,0,0,0,epsilon,1,strinGA);
	load([strinGA 'FD_completed'],'Lambda_P','ydelta_pred','eta_h','ydelta')
	y_pred_unrestrict			= [y_pred_unrestrict; ydelta_pred];
	Pii_unrestrict				= [Pii_unrestrict;diag(Lambda_P)];
	sigma_i						= [sigma_i;ydelta.*eta_h];

	%Build the KSS quadratic form
	Bii							=  (Pii_unrestrict-Pii_restrict);
	test						=  sum((y_pred_unrestrict-y_pred_restrict).^2)-sum(Bii.*sigma_i);
	test						=  test/r
	test_unadjusted				=  (sum((y_pred_unrestrict-y_pred_restrict).^2))/r
	
	%Build the new vec object with this updated C. The function C_build_cck allows me to compute C*v for any generic vector v. 
	vec 						= C_build_cck(ydelta_stack,GROUP_delta,group_stack,firmid_stack,id_stack,firmid_delta_f,firmid_delta,Pii_unrestrict,Pii_restrict);
	disp('verify --- MUST REPORT THE SAME AS test ')
	(ydelta_stack'*vec)/r
	
	%obtain the split sample estimates for group 0
	sel=GROUP_delta==0;
	[sigma_1_group0 sigma_2_group0 D_group0 adjustment_group0] = split_sample_inference_cck(id_movers(sel),firmid_delta(sel),firmid_delta_f(sel),ydelta_stack(sel),vec(sel),index_C,Ndelta,GROUP_delta,group_stack,firmid_stack,id_stack,firmid_delta_f,firmid_delta,Pii_unrestrict,Pii_restrict);
	
	%obtain the split sample estimates for group 1
	sel=GROUP_delta==1;
	[sigma_1_group1 sigma_2_group1 D_group1 adjustment_group1] = split_sample_inference_cck(id_movers(sel),firmid_delta(sel),firmid_delta_f(sel),ydelta_stack(sel),vec(sel),index_C,Ndelta,GROUP_delta,group_stack,firmid_stack,id_stack,firmid_delta_f,firmid_delta,Pii_unrestrict,Pii_restrict);
	
	%now combine
	sigma_1						= sigma_1_group0+sigma_1_group1;
	sigma_2						= sigma_2_group0+sigma_2_group1;
	D							= D_group0+D_group1;
	
	%Now get the standard error of the KSS quadratic form
	first_piece 				= 4*sum((vec.^2).*sigma_2);
	
	%Now the most complicated term (via Hutchinson, as usual)
	NSIM						= 1000;
	aux_SIM						= zeros(NSIM,1);
	XSIMUL						= randn(Ndelta,NSIM);
	VCM 						= parallel.pool.Constant(diag(sigma_1));

    parfor s=1:NSIM    
        xsimul			= XSIMUL(:,s);       
	 	left			= C_build_cck(xsimul,GROUP_delta,group_stack,firmid_stack,id_stack,firmid_delta_f,firmid_delta,Pii_unrestrict,Pii_restrict); 	
		xsimul			= VCM.Value*xsimul;
		right			= C_build_cck(xsimul,GROUP_delta,group_stack,firmid_stack,id_stack,firmid_delta_f,firmid_delta,Pii_unrestrict,Pii_restrict);
		aux_SIM(s)		= left'*VCM.Value*right;
		
    end 
    
   second_piece			= 2*mean(aux_SIM);
   Dsym					= 2*D + 2*D';
   third_piece			= 2*sigma_1'*Dsym*sigma_1;
   fourth_piece			= adjustment_group0 + adjustment_group1; %already have a times two in front.
   V					= first_piece - second_piece - third_piece - fourth_piece;
   V					= sqrt(V/r^2);
   
   %Done!
   student				= test/V;

end	