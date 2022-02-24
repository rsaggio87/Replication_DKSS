function [key_firm,r]	= dual_connected_set(filename,results_0,results_1,leave_out_level,type_of_algorithm,epsilon);

	%Start by loading results for leave out connected set of workers in group 0
	results						= results_0;
	y_group0					= results(:,1);
	firmid_group0				= results(:,2);
	J_group0					= max(firmid_group0);
	id_group0					= results(:,3);
	group0						= results(:,7);
	firmid_orig_group0			= results(:,4);
	id_orig_group0				= results(:,5);
	size_group0					= accumarray(firmid_group0,1);
	size_group0					= size_group0(firmid_group0);
	clear results
	
	
	%Start by loading results for leave out connected set of workers in group 1
	results                     = results_1;
	y_group1					= results(:,1);
	firmid_group1				= results(:,2);
	J_group1					= max(firmid_group1);
	id_group1					= results(:,3);
	group1						= results(:,7);
	firmid_orig_group1			= results(:,4);
	id_orig_group1				= results(:,5);
	size_group1					= accumarray(firmid_group1,1);
	size_group1					= size_group1(firmid_group1);
	clear results results_0 results_1
	
	%Find overlapping firms  
	sel_group0    				= ismember(firmid_orig_group0,firmid_orig_group1);
	sel_group1    				= ismember(firmid_orig_group1,firmid_orig_group0);
	
	%Summarize Dual Connected set of firms
	id							= [id_orig_group0(sel_group0); id_orig_group1(sel_group1)];
	[~,~,id]					= unique(id);
	firmid_orig					= [firmid_orig_group0(sel_group0); firmid_orig_group1(sel_group1)];
	[~,~,firmid]				= unique(firmid_orig);
	r							= max(firmid);
	y							= [y_group0(sel_group0); y_group1(sel_group1)];
	group						= [group0(sel_group0); group1(sel_group1)];
	SIZE						= accumarray(firmid,1);
	SIZE						= SIZE(firmid);				
	
	%Auxiliaries
	gcs 						= [NaN; id(1:end-1)];
	gcs 						= id~=gcs;
	lagfirmid					= [NaN; firmid(1:end-1)];
	lagfirmid(gcs==1)			= NaN; %%first obs for each worker
	stayer						= (firmid==lagfirmid);
	stayer(gcs==1)				= 1;
	stayer						= accumarray(id,stayer);
	T							= accumarray(id,1);
	stayer						= T==stayer;
	movers						= stayer~=1;
	movers						= movers(id);
	id_movers					= id(movers);
	[~,~,n]						= unique(id_movers);
	Nmovers						= max(n);
	
	
	s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s);
	disp(s);
	s=['Summarizing Dual Connected Set'];
	disp(s);
	s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s);
	disp(s);
	s=['mean wage: ' num2str(mean(y))];
	disp(s)
	s=['variance of wage: ' num2str(var(y))];
	disp(s)
	s=['# of Movers: ' num2str(Nmovers)];
	disp(s);
	s=['# of Firms: ' num2str(max(firmid))];
	disp(s);
	s=['# of Person Year Observations: ' num2str(size(y,1))];
	disp(s)
	s=['Share in GROUP 1: ' num2str(mean(group))];
	disp(s);
	
	clear gcs lagfirmd stayer T movers id_movers
	
	%Find largest firm in dual connected set
	key_firm					= find(SIZE==max(SIZE));
	key_firm					= mean((firmid_orig(key_firm)))
	
	%Renormalize the set of firm effects for group 0 according to key_firm 
	to_norm						= mean(firmid_group0(firmid_orig_group0==key_firm));
	firmid_group0 				= normalize_firm_effects(firmid_group0,to_norm);
	
	%Renormalize the set of firm effects for group 1 according to key_firm
	to_norm						= mean(firmid_group1(firmid_orig_group1==key_firm));
	firmid_group1				= normalize_firm_effects(firmid_group1,to_norm);
	
	%Estimate the two sets of firm effects again, save a .csv file that will be read by Stata.
	save_AKM 					=['results/gruop0_estimated_renormalized'];
	[fe_group0] 				= quick_AKM(y_group0,id_group0,firmid_group0,id_orig_group0,firmid_orig_group0,save_AKM);
	save_AKM 					=['results/group1_estimated_renormalized'];
	[fe_group1] 				= quick_AKM(y_group1,id_group1,firmid_group1,id_orig_group1,firmid_orig_group1,save_AKM);
	
	%Some auxiliaries
	R_group1					= splitapply(@max,sel_group1,firmid_group1); 
	disp('share of firms in leave out CS for group 1 that are in dual CS')
	mean(R_group1)
	dual_group1					= R_group1(firmid_group1);

	R_group0					= splitapply(@max,sel_group0,firmid_group0);
	disp('share of firms in leave out CS for group 0 that are in dual CS')
	mean(R_group0)
	dual_group0					= R_group0(firmid_group0);
	
	clear R_group1 R_group0
	
	%Now return to me the variance of firm effects for those in group 1 
	NT							= size(y_group1,1);
	D							= sparse(1:NT,id_group1',1);
	N							= size(D,2);
	F							= sparse(1:NT,firmid_group1',1);
	F							= F(dual_group1,:); 
	A_right						= (1/sqrt(size(F,1)))*[sparse(size(F,1),N) F];
	A_left						= A_right;
	deMean						= 1;
	theta_1  					= kss_quadratic_form(y_group1,id_group1,firmid_group1,leave_out_level,[],0,0,type_of_algorithm,A_right,A_left,deMean,epsilon)			

	%Now return to me the variance of firm effects for those in group 0
	NT							= size(y_group0,1);
	D							= sparse(1:NT,id_group0',1);
	N							= size(D,2);
	F							= sparse(1:NT,firmid_group0',1);
	F							= F(dual_group0,:); 
	A_right						= (1/sqrt(size(F,1)))*[sparse(size(F,1),N) F];
	A_left						= A_right;
	deMean						= 1;
	theta_0  					= kss_quadratic_form(y_group0,id_group0,firmid_group0,leave_out_level,[],0,0,type_of_algorithm,A_right,A_left,deMean,epsilon)			



	%Export file to run graphical analysis
	fe							= [fe_group0(sel_group0); fe_group1(sel_group1)];
	out							= [y,id,firmid_orig,firmid,group,fe, theta_1.*ones(size(fe,1),1),theta_0.*ones(size(fe,1),1)];
	s							= [filename '_dual.csv'];
	dlmwrite(s, out, 'delimiter', '\t', 'precision', 16); 
	
	
	
	

end



