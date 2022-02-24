function [R2_full, R2_short1,R2_short2,R2_full_KSS,R2_short1_KSS,R2_short2_KSS]	= kss_rakim(id,firmid,y,controls,type_of_algorithm,STRINGA,cara,scale);

%Auxiliaries
	if nargin == 7 
	scale					= 1000;
	end
	
%Create auxiliaries
	NT				= size(y,1);
	count			= ones(NT,1);
	gcs 			= cell2mat(accumarray(id,count,[],@(x){cumsum(x)}));
	
%Define Fdelta for Network 1	
	sel				= (gcs == 1 | gcs ==2);
	firmid_orig1	= firmid(sel); %T=2.
	[~,~,firmid_use]= unique(firmid_orig1);
	gcs_use			= gcs(sel);
	firmid_lag_1	= firmid_use(gcs_use==1);
	firmid_con_1	= firmid_use(gcs_use==2);
	firmid_orig_lag_1= firmid_orig1(gcs_use==1);
	firmid_orig_con_1= firmid_orig1(gcs_use==2);
	
	F1= sparse((1:max(id))',firmid_lag_1',1,max(id),max(firmid_use));
	F2= sparse((1:max(id))',firmid_con_1',1,max(id),max(firmid_use));
	FDELTA1			= F2-F1; 
	
%Define Fdelta for Network 2			
	sel				= (gcs == 2 | gcs ==3);
	firmid_orig2	= firmid(sel);
	[~,~,firmid_use]= unique(firmid_orig2);
	gcs_use			= gcs(sel);
	firmid_lag_2	= firmid_use(gcs_use==2);
	firmid_con_2	= firmid_use(gcs_use==3);
	firmid_orig_lag_2= firmid_orig2(gcs_use==2);
	firmid_orig_con_2= firmid_orig2(gcs_use==3);
	F1				= sparse((1:max(id))',firmid_lag_2',1,max(id),max(firmid_use));
	F2				= sparse((1:max(id))',firmid_con_2',1,max(id),max(firmid_use));
	FDELTA2			= F2-F1; 
	ydelta			= y(gcs==3)-y(gcs==2);
	
	
%Summarize starting point 
	s=['Summarize Starting Point:'];
	disp(s)
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	disp(s)
	s=['Mean of Delta y: ' num2str(mean(ydelta))];
	disp(s);
	s=['Variance of Delta y: ' num2str(var(ydelta))];
	disp(s);
	s=['# of workers: ' int2str(max(id))];
	disp(s);
	s=['# of firms: ' int2str(max(firmid))];
	disp(s);
	s=['# of firms in Network 1: ' int2str(size(FDELTA1,2))];
	disp(s);
	s=['# of firms in Network 2: ' int2str(size(FDELTA2,2))];
	disp(s);
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s);
	s=['STARTING PRUNING.........'];
	disp(s)
	
	
%Save stuff for table 1
	N_individuals		= zeros(3,1);
	N_firms				= zeros(3,1);
	lagged_firms		= zeros(3,1);
	present_firms		= zeros(3,1);
	mean_wage			= zeros(3,1);
	var_wage			= zeros(3,1);

	N_individuals(1,1)	= max(id);
	N_firms(1,1)		= max(firmid);
	lagged_firms(1,1)	= size(FDELTA1,2);
	present_firms(1,1)	= size(FDELTA2,2);
	mean_wage(1,1)		= mean(ydelta);
	var_wage(1,1)		= var(ydelta);
	
%Run pruning algorithm
	tic
	[y,id,firmid,gcs,controls,ydelta,X,FDELTA1,FDELTA2,firmid_orig1,firmid_orig2,firmid_lag_1,firmid_con_1,firmid_lag_2,firmid_con_2,firmid_DICTIO1,firmid_DICTIO2,cara]=pruning_rakim(y,id,firmid,controls,cara);	
	toc
			
%Find the identified sets	
	identified1 					= identify_rakim(X,FDELTA1,FDELTA2,1);
	identified2 					= identify_rakim(X,FDELTA1,FDELTA2,2);
	
%Tell me which observation is 100% identified
	identified_lag_1				= identified1(firmid_lag_1);
	identified_con_1				= identified1(firmid_con_1);
	identified_lag_2				= identified2(firmid_lag_2);
	identified_con_2				= identified2(firmid_con_2);
	identified						= (identified_lag_1 == 1 & identified_con_1 == 1 & identified_lag_2 == 1 & identified_con_2 == 1);
	disp('# of identified movers')   	
	sum(identified)
	disp('share of identified movers')   	
	mean(identified)
	
	id_identified					= find(identified);
	
%Run model in levels
	[theta,sigma_i]				    = estimate_rakim_levels(y,id,firmid,controls,id_identified,scale,STRINGA,X);
	
	R2_full=0;
	R2_short1=0;
	R2_short2=0;
	R2_full_KSS=0;
	R2_short1_KSS=0;
	R2_short2_KSS=0;
	
	if 0 == 1	
	
%Backup which firms are identified and in which network
	[~,~,firmidNetwork1]			= unique(firmid_orig1);
	firmid_orig1_coll				= accumarray(firmidNetwork1,firmid_orig1,[],@(x)mean(x));	
	elist1							= [firmid_orig1_coll (1:size(identified1,1))' identified1];
	[~,~,firmidNetwork2]			= unique(firmid_orig2);
	firmid_orig2_coll				= accumarray(firmidNetwork2,firmid_orig2,[],@(x)mean(x));	
	elist2							= [firmid_orig2_coll (1:size(identified2,1))' identified2];
	LIST_BASE 						= array2table(elist1,...
    								'VariableNames',{'firmid_original','norm_id1','identified1'});
    LIST_SEL 						= array2table(elist2,...
        							'VariableNames',{'firmid_original','norm_id2','identified2'});
    merge							= outerjoin(LIST_BASE,LIST_SEL);
    merge 							= table2array(merge);
    merge							= merge(~any(isnan(merge),2),:); %not merged cases are firms present only in one network 
	identified_dual					= (merge(:,3)+merge(:,6)==2);
	dual_list						= [merge(:,1) merge(:,2) merge(:,5) identified_dual]; %original firm id; %firm id normalized in Network 1 % firm id normalized in Network 2 %identified in both networks
	
	clear merge identified_dual LIST_BASE LIST_SEL elist2 firmidNetwork2 elist1 firmidNetwork1

%Construct the A matrices for the model in differences
	NSEL							= sum(identified);
	A_diff_network1					= [FDELTA1(identified,:).*identified1' sparse(NSEL,size(X,2)-size(FDELTA1,2))];
	A_diff_network2					= [sparse(NSEL,size(FDELTA1,2)) FDELTA2(identified,:).*identified2' sparse(NSEL,size(X,2)-size(FDELTA1,2)-size(FDELTA2,2))];

%Stuff for normalization
	%LIsTa							= [firmid_DICTIO2 firmid_orig2_coll];
	%sel								= find(LIsTa(:,1)==firmid_MAIN);
	%firmid_PVR						= LIsTa(sel,2)
									
%Return me a list of firms that are DUAL CONNECTED AND IDENTIFIED 
	firm_size						= accumarray(firmid,1);
	firm_size						= [(1:size(firm_size))' firm_size];
	LIST_BASE 						= array2table(firm_size,...
    									'VariableNames',{'firmid_original','firm_size'});
    LIST_SEL 						= array2table(dual_list,...
        								'VariableNames',{'firmid_original','firmid_1','firmid_2','dual_identified'});
    merge							= outerjoin(LIST_BASE,LIST_SEL);
    merge 							= table2array(merge);
    merge							= merge(~any(isnan(merge),2),:);    																	
	dual_list						= [merge(:,1) merge(:,4) merge(:,5) merge(:,6) merge(:,2)]; %original firm id; %firm id normalized in Network 1 % firm id normalized in Network 2 %identified in both networks % size
	
	%sel							= dual_list(:,1)==firmid_PVR;
	%MYLIST							= dual_list(sel,:)
	%if size(MYLIST,1)== 0 | MYLIST(:,4) == 0 
	%error('normalization will not work')
	%end
	
	sel							    = (dual_list(:,4)==1);
	dual_list						= dual_list(sel,:);
	clear merge LIST_BASE  firm_size sel

%Export information where we report: (i) original firmid (ii) psi/lambda/psi_AKM (iii) whether the coefficient is identified in the PVR MODEL.
	xx					    		= X'*X;
	xy								= X'*ydelta;
	b								= pcg(xx,xy,1e-10,10000);	%overkill
	
	controls_effects				= b(size(FDELTA1,2)+size(FDELTA2,2)+1:end);
	out								= [controls_effects];
	s								= ['tables/controls_effects_' STRINGA '.csv'];
	dlmwrite(s, out, 'delimiter', '\t', 'precision', 16); 
	
	lambda							= b(1:size(FDELTA1,2));
	%lambda							= lambda-lambda(MYLIST(1,2)); %normalization
	out								= [firmid_DICTIO1 lambda identified1];
	s								= ['tables/RAKIM_LAMBDA_' STRINGA '.csv'];
	dlmwrite(s, out, 'delimiter', '\t', 'precision', 16); 
	
	psi								= b(size(FDELTA1,2)+1:size(FDELTA1,2)+size(FDELTA2,2));
	%psi							= psi-psi(MYLIST(1,3)); %normalization
	controls_1						= controls(gcs==2,:);
	controls_2						= controls(gcs==3,:);
	controls_DELTA  				= controls_2-controls_1;
	X_restrict						= [FDELTA2 controls_DELTA ones(size(FDELTA2,1),1)];
	xx					    		= X_restrict'*X_restrict;
	xy								= X_restrict'*ydelta;
	b								= pcg(xx,xy,1e-10,10000);	%overkill
	psi_AKM							= b(1:size(FDELTA2,2));
	%psi_AKM						= psi_AKM-psi_AKM(MYLIST(1,3));
	out								= [firmid_DICTIO2 psi identified2 psi_AKM];
	s								= ['tables/RAKIM_PSI_' STRINGA '.csv'];
	dlmwrite(s, out, 'delimiter', '\t', 'precision', 16);
	clear xx xy X_restrict controls_1 controls_2 controls_DELTA
 
%Now construct the corresponding A matrices for levels in dual connected
	disp('Total # of dual connected and identified firms')
	N_dual							= size(dual_list,1)
	sel								= dual_list(:,2)';		
	A_levels_network1				= sparse((1:N_dual)',sel,1,N_dual,size(X,2));
	A_levels_network1				= repelem(A_levels_network1,dual_list(:,5),1); %span it w.r.t. size of each firm
	sel								= size(FDELTA1,2) + dual_list(:,3)';
	A_levels_network2				= sparse((1:N_dual)',sel,1,N_dual,size(X,2));	
	A_levels_network2				= repelem(A_levels_network2,dual_list(:,5),1); %span it w.r.t. size of each firm	
	    
%Now find the leverages and corresponding Bii of the Full Model
	tic
	[R2_full, R2_full_KSS,Pii_unrestricted,theta,sigma_i]	= estimate_rakim(X,ydelta,A_diff_network1,A_diff_network2,A_levels_network1,A_levels_network2,type_of_algorithm,STRINGA,scale);
	disp('completed estimation of the full model')
	toc
	X_unrestrict					= X;
	
%Now on this main sample estimate a model that has firmid x gender
	[R2_interacted R2_interacted_KSS]= interacted_model(y,Pii_unrestricted,sigma_i,id,firmid,controls,cara)	

%Now find the leverages of the short regression (Network 1)
	controls_1						= controls(gcs==1,:);
	controls_2						= controls(gcs==2,:);
	controls_DELTA  				= controls_2-controls_1;
	X								= [FDELTA1 controls_DELTA ones(size(FDELTA1,1),1)];
	tic
	[R2_short1, R2_short1_KSS]		= estimate_rakim(X,ydelta,A_diff_network1,A_diff_network2,A_levels_network1,A_levels_network2,type_of_algorithm,STRINGA,scale);
    disp('completed estimation of the short model based on Network 1')
	toc

%Now find the leverages of the short regression (Network 2)
	controls_1						= controls(gcs==2,:);
	controls_2						= controls(gcs==3,:);
	controls_DELTA  				= controls_2-controls_1;
	X								= [FDELTA2 controls_DELTA ones(size(FDELTA1,1),1) ];
	tic
	[R2_short2, R2_short2_KSS,Pii_restricted] = estimate_rakim(X,ydelta,A_diff_network1,A_diff_network2,A_levels_network1,A_levels_network2,type_of_algorithm,STRINGA,scale);
	disp('completed estimation of the short model based on Network 2')
	toc
	X_restrict						= X;
	
%Check the # of overlapping firms
	present_in_network_1		   = (gcs == 1| gcs == 2);
	present_in_network_2		   = (gcs == 2| gcs == 3);
	present_in_network_1		   = accumarray(firmid,present_in_network_1,[],@(x)max(x));	
	present_in_network_2		   = accumarray(firmid,present_in_network_2,[],@(x)max(x));
	present_in_both_network		   = (present_in_network_1==1 & present_in_network_2 == 1);
	present_in_both_network		   = sum(present_in_both_network);
	
%Construct the sample of identified moves
	identified					  = identified(id);
	id_ident					  = id(identified);
	[~,~,id_ident]				  = unique(id_ident);	
	firmid_ident				  = firmid(identified);
	[~,~,firmid_ident]			  = unique(firmid_ident);
	y_ident						  = y(identified);
	count						  = ones(length(y_ident),1);
	gcs_iden 					  = cell2mat(accumarray(id_ident,count,[],@(x){cumsum(x)}));
	ydelta_ident				  = y_ident(gcs_iden==3)-y_ident(gcs_iden==2);
	firms_net1					  = firmid_ident(gcs_iden==2 | gcs_iden==1);
	[~,~,firms_net1]			  = unique(firms_net1);
	firms_net2					  = firmid_ident(gcs_iden==3 | gcs_iden==2);
	[~,~,firms_net2]			  = unique(firms_net2);
	
	
	s=['Summarizing Results: Observations with Pii<1'];
	disp(s)
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	disp(s)
	s=['Mean of Delta y: ' num2str(mean(ydelta))];
	disp(s);
	s=['Variance of Delta y: ' num2str(var(ydelta))];
	disp(s);
	s=['# of workers: ' int2str(max(id))];
	disp(s);
	s=['# of firms: ' int2str(max(firmid))];
	disp(s);
	s=['# of overlapping firms: ' int2str(present_in_both_network)];
	disp(s);
	s=['# of firms in Network 1: ' int2str(size(FDELTA1,2))];
	disp(s);
	s=['# of firms in Network 2: ' int2str(size(FDELTA2,2))];
	disp(s);
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	s=['Summarizing Results: Observations with Identified contemporaneous and lagged firm effects'];
	disp(s)
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	disp(s)
	s=['Mean of Delta y: ' num2str(mean(ydelta_ident))];
	disp(s);
	s=['Variance of Delta y: ' num2str(var(ydelta_ident))];
	disp(s);
	s=['# of workers: ' int2str(max(id_ident))];
	disp(s);
	s=['# of firms: ' int2str(max(firmid_ident))];
	disp(s);
	s=['# of firms in Network 1: ' int2str(max(firms_net1))];
	disp(s);
	s=['# of firms in Network 2: ' int2str(max(firms_net2))];
	disp(s);
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)	
	s=['Plug in R2:'];
	disp(s)
	s=['R2 full: ' num2str(R2_full)];
	disp(s);
	s=['R2 short - Network 1: ' num2str(R2_short1)];
	disp(s);	
	s=['R2 short - Network 2: ' num2str(R2_short2)];
	disp(s);		
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)	
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)	
	s=['KSS R2:'];
	disp(s)
	s=['R2 full: ' num2str(R2_full_KSS)];
	disp(s);
	s=['R2 short - Network 1: ' num2str(R2_short1_KSS)];
	disp(s);	
	s=['R2 short - Network 2: ' num2str(R2_short2_KSS)];
	disp(s);		
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	s=['Variance Decomposition --- MODEL IN FD -- PLUG IN:'];
	disp(s)
	s=['Variance of Lag Firm effects     ' num2str(theta(1,1))];
	disp(s);
	s=['Variance of Current Firm effects     ' num2str(theta(2,1))];
	disp(s);
	s=['Covariance of Lag, Current Firm effects     ' num2str(theta(3,1))];
	disp(s);
	s=['Correlation of Lag, Current Firm effects     ' num2str(theta(3,1)/(sqrt(theta(2,1))*sqrt(theta(1,1))))];
	disp(s);	
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	s=['Variance Decomposition --- MODEL IN FD -- KSS:'];
	disp(s)
	s=['Variance of Lag Firm effects     ' num2str(theta(1,2))];
	disp(s);
	s=['Variance of Current Firm effects     ' num2str(theta(2,2))];
	disp(s);
	s=['Covariance of Lag, Current Firm effects     ' num2str(theta(3,2))];
	disp(s);
	s=['Correlation of Lag, Current Firm effects     ' num2str(theta(3,2)/(sqrt(theta(2,2))*sqrt(theta(1,2))))];
	disp(s);	
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	s=['Variance Decomposition --- DUAL FIRMS -- PLUG IN:'];
	disp(s)
	s=['Variance of Lag Firm effects     ' num2str(theta(4,1))];
	disp(s);
	s=['Variance of Current Firm effects     ' num2str(theta(5,1))];
	disp(s);
	s=['Covariance of Lag, Current Firm effects     ' num2str(theta(6,1))];
	disp(s);
	s=['Correlation of Lag, Current Firm effects     ' num2str(theta(6,1)/(sqrt(theta(4,1))*sqrt(theta(5,1))))];
	disp(s);	
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	s=['Variance Decomposition --- DUAL FIRMS -- KSS:'];
	disp(s)
	s=['Variance of Lag Firm effects     ' num2str(theta(4,2))];
	disp(s);
	s=['Variance of Current Firm effects     ' num2str(theta(5,2))];
	disp(s);
	s=['Covariance of Lag, Current Firm effects     ' num2str(theta(6,2))];
	disp(s);
	s=['Correlation of Lag, Current Firm effects     ' num2str(theta(6,2)/(sqrt(theta(4,2))*sqrt(theta(5,2))))];
	disp(s);	
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	
	save(['tables/ALL_results' STRINGA])
	
	%Stuff for table 1
	N_individuals(2,1)	= max(id);
	N_firms(2,1)		= max(firmid);
	lagged_firms(2,1)	= size(FDELTA1,2);
	present_firms(2,1)	= size(FDELTA2,2);
	mean_wage(2,1)		= mean(ydelta);
	var_wage(2,1)		= var(ydelta);
	
	N_individuals(3,1)	= max(id_ident);
	N_firms(3,1)		= max(firmid_ident);
	lagged_firms(3,1)	= max(firms_net1);
	present_firms(3,1)	= max(firms_net2);
	mean_wage(3,1)		= mean(ydelta_ident);
	var_wage(3,1)		= var(ydelta_ident);
	
	%Export Summary stats
	out 				= [-9999999; -9999999; N_individuals(1,1); N_firms(1,1); lagged_firms(1,1); present_firms(1,1); mean_wage(1,1); var_wage(1,1); -9999999; N_individuals(2,1); N_firms(2,1); lagged_firms(2,1); present_firms(2,1); mean_wage(2,1); var_wage(2,1); -9999999; N_individuals(3,1); N_firms(3,1); mean_wage(3,1); var_wage(3,1)];
	s					= ['tables/TABLE_1' STRINGA '.csv'];
	dlmwrite(s, out, 'delimiter', '\t', 'precision', 16);
	
	%Stuff for table 2
	out 				= [var(ydelta); -9999999; R2_short1; R2_short2; R2_full; R2_interacted; -9999999; R2_short1_KSS; R2_short2_KSS; R2_full_KSS; R2_interacted_KSS];
	s					= ['tables/TABLE_2' STRINGA '.csv'];
	dlmwrite(s, out, 'delimiter', '\t', 'precision', 16);
	
	%Stuff for table 3
	out 				= [var(ydelta_ident); -9999999; theta(1,1); theta(2,1); 2*theta(3,1); -9999999; theta(1,2);; theta(2,2); 2*theta(3,2); ];
	s					= ['tables/TABLE_3' STRINGA '.csv'];
	dlmwrite(s, out, 'delimiter', '\t', 'precision', 16);
	
	%Stuff for table 4
	out 				= [max(firmid); N_dual; -9999999; theta(4,1); theta(5,1); theta(6,1); theta(6,1)/(sqrt(theta(4,1))*sqrt(theta(5,1))); -9999999; theta(4,2); theta(5,2); theta(6,2); theta(6,2)/(sqrt(theta(4,2))*sqrt(theta(5,2)));];
	s					= ['tables/TABLE_4' STRINGA '.csv'];
	dlmwrite(s, out, 'delimiter', '\t', 'precision', 16);
	end

end	
