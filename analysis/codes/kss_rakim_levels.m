function theta	= kss_rakim_levels(y,id,firmid,lagfirmid,firmid_orig,lagfirmid_orig,controls,cara,scale,STRINGA)
theta=0

%Now tell me which firm is dual identified (both as a current and a lag firm effect)
	firmid_orig1_coll				= accumarray(firmid,firmid_orig,[],@(x)mean(x)); %%cross walk: firmids
	firmi_size1						= accumarray(firmid,1);
	elist1							= [firmid_orig1_coll (1:size(firmi_size1,1))' firmi_size1];
	firmid_orig2_coll				= accumarray(lagfirmid,lagfirmid_orig,[],@(x)mean(x)); %%cross walk: lag firmids
	elist2							= [firmid_orig2_coll (1:size(firmid_orig2_coll,1))'];
	LIST_BASE 						= array2table(elist1,...
    								'VariableNames',{'firmid_original','norm_id1', 'size'});
    LIST_SEL 						= array2table(elist2,...
        							'VariableNames',{'firmid_original','norm_id2',  });
    merge							= outerjoin(LIST_BASE,LIST_SEL);
    merge 							= table2array(merge);
    merge							= merge(~any(isnan(merge),2),:); %not merged cases are firms present only in one network 
	dual_list						= [merge(:,1) merge(:,2) merge(:,5) ones(size(merge,1),1) merge(:,3)]; %original firm id; %firm id normalized in Network 1 % firm id normalized in Network 2 %identified in both networks %firm size
	
%Now merge with original firm size
	s								= ['../build/src/RAKIM/FIRM_SIZE.csv']
	list_size						= importdata(s);
	LIST_SIZE 						= array2table(list_size,...
    								'VariableNames',{'firmid_original','true_firm_size','settore'});
    LIST_DUAL 						= array2table(dual_list,...
    								'VariableNames',{'firmid_original','firmid_1','firmid_2','identified','old_firm_size'});
    merge							= outerjoin(LIST_DUAL,LIST_SIZE);
    merge 							= table2array(merge);
    merge							= merge(~any(isnan(merge),2),:); %don't consider non-merged cases
	dual_list						= [merge(:,1) merge(:,2) merge(:,3) ones(size(merge,1),1) merge(:,7) merge(:,8)]; %original firm id; %firm id normalized in Network 1 % firm id normalized in Network 2 %identified in both networks % real firm size %sector													
	N_dual							= size(dual_list,1)
	clear merge LIST_BASE LIST_SEL elist2 firmidNetwork2 elist1 firmid_orig1_coll firmid_orig2_coll firmi_size1
	
%This file assumes that the imported data corresponds to job-year observations for which firm effects are identified
	NT								= size(y,1);
	identified						= (1:NT)';
	identified						= identified>0;


%Run DAKM model in levels
	if 1 == 1
	[R2, R2_KSS,Pii_R2,theta,sigma_i,xb] = estimate_rakim_levels(3,y,id,firmid,lagfirmid,firmid_orig,lagfirmid_orig,controls,identified,dual_list,scale,STRINGA);
	out 					= [sigma_i(identified)];
	s						= ['../build/src/RAKIM/sigmaI' STRINGA '.csv']
	dlmwrite(s, out, 'precision', 16);
	end


%Analysis by Sector
	s						= ['../build/src/RAKIM/sigmaI' STRINGA '.csv']
	sigma_i					= importdata(s);
	theta_sector			= estimate_rakim_sector(y,sigma_i,id,firmid,lagfirmid,firmid_orig,lagfirmid_orig,controls,identified,dual_list,scale*2,STRINGA);
		 
%Run DAKM model interacting gender with firmids
	[R2_interacted R2_interacted_KSS]                             = interacted_model_levels(3,y,Pii_R2,sigma_i,id,firmid,lagfirmid,controls,cara);

%Run origin effects model only in levels (uncomment to print results)
	[R2_short1,R2_KSS_short1,Pii_R2_lagg,~,sigma_i_lagg]		  = estimate_rakim_levels(1,y,id,firmid,lagfirmid,firmid_orig,lagfirmid_orig,controls,identified,dual_list,scale,STRINGA);  
	[R2_interacted_lagg R2_interacted_KSS_lagg]	  				  = interacted_model_levels(1,y,Pii_R2_lagg,sigma_i_lagg,id,firmid,lagfirmid,controls,cara);
    

%Run AKM effects model in levels
	[R2_short2,R2_KSS_short2,Pii_R2_AKM,theta_AKM,sigma_i_AKM]	  = estimate_rakim_levels(2,y,id,firmid,lagfirmid,firmid_orig,lagfirmid_orig,controls,identified,dual_list,scale,STRINGA);
    [R2_interacted_AKM R2_interacted_KSS_AKM]	  			      = interacted_model_levels(2,y,Pii_R2_AKM,sigma_i_AKM,id,firmid,lagfirmid,controls,cara);
    

%Turn next lines on in order to obtain SEs in Figure 5, Panel(a).
if 	0 == 1
	s						= ['../build/src/RAKIM/sigmaI' STRINGA '.csv']
	sigma_i					= importdata(s); %get the KSS sigma_i
	
	s						= ['../build/src/RAKIM/LINCOM.csv'] %see "build_lincom"
	list_lincom				= importdata(s);
	LIST_LINCOM 			= array2table(list_lincom,...
    						  'VariableNames',{'firmid_original','true_firm_size','VA_L'});
    						  
    %Read cross-walk firmid
    for ll=1:2
    
    if ll==1 %firmid
    firmid_orig1_coll		= accumarray(firmid,firmid_orig,[],@(x)mean(x)); %%cross walk: firmids
	elist1					= [firmid_orig1_coll (1:size(firmid_orig1_coll,1))'];					  
    LIST_BASE 				= array2table(elist1,...
    								'VariableNames',{'firmid_original','norm_id1'});
    end
    
    if ll==2 %lagfirmid
	firmid_orig2_coll		= accumarray(lagfirmid,lagfirmid_orig,[],@(x)mean(x)); %%cross walk: lag firmids
	elist2					= [firmid_orig2_coll (1:size(firmid_orig2_coll,1))'];
    LIST_BASE 				= array2table(elist2,...
     							'VariableNames',{'firmid_original','norm_id2',  });
    end 
    														  						  
    merge					= outerjoin(LIST_BASE,LIST_LINCOM);
    merge 					= table2array(merge);
    merge					= merge(~any(isnan(merge),2),:); %don't consider non-merged cases
	list_lincom				= [merge(:,1) merge(:,2) merge(:,4) merge(:,5)]; %original firm id; %firm id normalized in Network 1 % firm id normalized in Network 2 %identified in both networks % real firm size %VA													
	[slope,KSS_SE]			= lincom_RAKIM(y,id,firmid,lagfirmid,controls,list_lincom,sigma_i,ll);
    end

end
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	disp(s)
	s=['AKM --- Plug-in Decomposition of Poaching Wages'];
	disp(s)	
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	disp(s)
	s=['Variance of Person Effects: ' num2str(theta_AKM(1,1))];
	disp(s);
	s=['Variance of Firm Effects: ' num2str(theta_AKM(2,1))];
	disp(s);
	s=['CoVariance of Person, Firm Effects: ' num2str(theta_AKM(4,1))];
	disp(s);
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	disp(s)
	s=['AKM --- KSS Decomposition of Poaching Wages'];
	disp(s)	
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	disp(s)
	s=['Variance of Person Effects: ' num2str(theta_AKM(1,2))];
	disp(s);
	s=['Variance of Firm Effects: ' num2str(theta_AKM(2,2))];
	disp(s);
	s=['CoVariance of Person, Firm Effects: ' num2str(theta_AKM(4,2))];
	disp(s);
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	
%Summarize Results
	s=['Estimation Sample for DAKM Model with Pii<1 and identified job-year observations'];
	disp(s)
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	disp(s)
	s=['Mean of y: ' num2str(mean(y))];
	disp(s);
	s=['Variance of y: ' num2str(var(y))];
	disp(s);
	s=['# of person x job observations: ' int2str(size(y,1))];
	disp(s);
	s=['# of workers: ' int2str(max(id))];
	disp(s);
	s=['# of cont.  firms: ' int2str(max(firmid))];
	disp(s);
	s=['# of lagged firms: ' int2str(max(lagfirmid))];
	disp(s);	
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	s=['Plug in R2:'];
	disp(s)
	s=['R2 short - Lagged: ' num2str(R2_short1)];
	disp(s);	
	s=['R2 short - Contemp: ' num2str(R2_short2)];
	disp(s);
	s=['R2 full: ' num2str(R2)];
	disp(s);	
	s=['R2 interacted: ' num2str(R2_interacted)];
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
	s=['R2 short - Lagged: ' num2str(R2_KSS_short1)];
	disp(s);	
	s=['R2 short - Lagged, Gender Inter: ' num2str(R2_interacted_KSS_lagg)];
	disp(s);	
	s=['R2 short - Contemp: ' num2str(R2_KSS_short2)];
	disp(s);
	s=['R2 short - Contemp, Gender Inter: ' num2str(R2_interacted_KSS_AKM)];
	disp(s);	
	s=['R2 full: ' num2str(R2_KSS)];
	disp(s);	
	s=['R2 interacted: ' num2str(R2_interacted_KSS)];
	disp(s);	
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	s=['Variance Decomposition -- IDENTIFIED JOB YEAR OBS -- PLUG IN:'];
	disp(s)
	s=['Variance of Person effects     ' num2str(theta(1,1))];
	disp(s);
	s=['Variance of Current Firm effects     ' num2str(theta(2,1))];
	disp(s);
	s=['Variance of Lag Firm effects     ' num2str(theta(3,1))];
	disp(s);
	s=['Covariance of Person, Current Firm effects     ' num2str(theta(4,1))];
	disp(s);
	s=['Covariance of Person, Lag Firm effects     ' num2str(theta(5,1))];
	disp(s);
	s=['Covariance of Lag, Current Firm effects     ' num2str(theta(6,1))];
	disp(s);
	s=['Correlation of Person, Current Firm effects     ' num2str(theta(4,1)/(sqrt(theta(1,1))*sqrt(theta(2,1))))];
	disp(s);
	s=['Correlation of Person, Lag Firm effects     ' num2str(theta(5,1)/(sqrt(theta(1,1))*sqrt(theta(3,1))))];
	disp(s);	
	s=['Correlation of Lag, Current Firm effects     ' num2str(theta(6,1)/(sqrt(theta(2,1))*sqrt(theta(3,1))))];
	disp(s);	
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	s=['Variance Decomposition --- -- IDENTIFIED JOB YEAR OBS -- -- KSS:'];
	disp(s)
	s=['Variance of Person effects     ' num2str(theta(1,2))];
	disp(s);
	s=['Variance of Current Firm effects     ' num2str(theta(2,2))];
	disp(s);
	s=['Variance of Lag Firm effects     ' num2str(theta(3,2))];
	disp(s);
	s=['Covariance of Person, Current Firm effects     ' num2str(theta(4,2))];
	disp(s);
	s=['Covariance of Person, Lag Firm effects     ' num2str(theta(5,2))];
	disp(s);
	s=['Covariance of Lag, Current Firm effects     ' num2str(theta(6,2))];
	disp(s);
	s=['Correlation of Person, Current Firm effects     ' num2str(theta(4,2)/(sqrt(theta(1,2))*sqrt(theta(2,2))))];
	disp(s);
	s=['Correlation of Person, Lag Firm effects     ' num2str(theta(5,2)/(sqrt(theta(1,2))*sqrt(theta(3,2))))];
	disp(s);	
	s=['Correlation of Lag, Current Firm effects     ' num2str(theta(6,2)/(sqrt(theta(2,2))*sqrt(theta(3,2))))];
	disp(s);
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	s=['Variance Decomposition --- DUAL FIRMS -- PLUG IN:'];
	disp(s)
	s=['Variance of Current Firm effects     ' num2str(theta(7,1))];
	disp(s);
	s=['Variance of Lag Firm effects     ' num2str(theta(8,1))];
	disp(s);
	s=['Covariance of Lag, Current Firm effects     ' num2str(theta(9,1))];
	disp(s);
	s=['Correlation of Lag, Current Firm effects     ' num2str(theta(9,1)/(sqrt(theta(7,1))*sqrt(theta(8,1))))];
	disp(s);	
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	s=['Variance Decomposition --- DUAL FIRMS -- KSS:'];
	disp(s)
	s=['Variance of Current Firm effects     ' num2str(theta(7,2))];
	disp(s);
	s=['Variance of Lag Firm effects     ' num2str(theta(8,2))];
	disp(s);
	s=['Covariance of Lag, Current Firm effects     ' num2str(theta(9,2))];
	disp(s);
	s=['Correlation of Lag, Current Firm effects     ' num2str(theta(9,2)/(sqrt(theta(7,2))*sqrt(theta(8,2))))];
	disp(s);	
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)	
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	
	NT_sum(2,1)			= size(y,1);
	N_individuals(2,1)	= max(id);
	[~,~,aux]			= unique([firmid_orig;lagfirmid_orig]);
	N_firms(2,1)		= max(aux);
	[~,~,aux]			= unique(lagfirmid_orig);
	lagged_firms(2,1)	= max(aux);
	[~,~,aux]			= unique(firmid_orig);
	present_firms(2,1)	= max(aux);
	mean_wage(2,1)		= mean(y);
	var_wage(2,1)		= var(y);

%Export the tables
	out 				= [-9999999; -9999999; NT_sum(1,1); N_individuals(1,1); N_firms(1,1); lagged_firms(1,1); present_firms(1,1); mean_wage(1,1); var_wage(1,1); -9999999; NT_sum(2,1); N_individuals(2,1); N_firms(2,1); lagged_firms(2,1); present_firms(2,1); mean_wage(2,1); var_wage(2,1)];
	s					= ['tables/TABLE_1' STRINGA '.csv'];
	dlmwrite(s, out, 'delimiter', '\t', 'precision', 16);
	
	%Stuff for table 2
	out 				= [-9999999; R2_short1; R2_short2; R2; R2_interacted; -9999999; -9999999; R2_KSS_short1; R2_KSS_short2; R2_KSS; R2_interacted_KSS];
	s					= ['tables/TABLE_2' STRINGA '.csv']; 
	dlmwrite(s, out, 'delimiter', '\t', 'precision', 16);

	%Stuff for table 4
	out 				= [var(y); -9999999; -9999999; theta(1,1); theta(2,1); theta(3,1); 2*theta(4,1); 2*theta(5,1); 2*theta(6,1); -9999999; -9999999; theta(1,2); theta(2,2); theta(3,2); 2*theta(4,2); 2*theta(5,2); 2*theta(6,2);];
	s					= ['tables/TABLE_5' STRINGA '.csv'];
	dlmwrite(s, out, 'delimiter', '\t', 'precision', 16);
	
	%Stuff for table 5
	out 				= [N_dual; -9999999; -9999999; theta(7,1); theta(8,1); theta(9,1); theta(9,1)/(sqrt(theta(7,1))*sqrt(theta(8,1))); -9999999; -9999999; theta(8,2); theta(9,2); theta(9,2)/(sqrt(theta(7,2))*sqrt(theta(8,2)));];
	s					= ['tables/TABLE_6' STRINGA '.csv'];
	dlmwrite(s, out, 'delimiter', '\t', 'precision', 16);		
	end		
