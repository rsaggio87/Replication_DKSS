function identified	= kss_rakim_build(id,firmid,lagfirmid,y,controls,cara,scale,STRINGA)


%Create auxiliaries
	NT				= size(y,1);
	count			= ones(NT,1);
	gcs 			= cell2mat(accumarray(id,count,[],@(x){cumsum(x)}));
		
%Save stuff for table 1
	NT_sum				= zeros(3,1);
	N_individuals		= zeros(3,1);
	N_firms				= zeros(3,1);
	lagged_firms		= zeros(3,1);
	present_firms		= zeros(3,1);
	mean_wage			= zeros(3,1);
	var_wage			= zeros(3,1);
	
	NT_sum(1,1)			= size(y,1);
	N_individuals(1,1)	= max(id);
	[~,~,aux]			= unique([firmid;lagfirmid]);
	N_firms(1,1)		= max(aux);
	[~,~,aux]			= unique(lagfirmid);
	lagged_firms(1,1)	= max(aux);
	[~,~,aux]			= unique(firmid);
	present_firms(1,1)	= max(aux);
	mean_wage(1,1)		= mean(y);
	var_wage(1,1)		= var(y);
	
%Run pruning algorithm
	tic
	[y,id,firmid,lagfirmid,controls,id_orig,firmid_orig,lagfirmid_orig,cara]=pruning_rakim_levels(y,id,firmid,lagfirmid,controls,cara);	
	toc
	
%Save the results of the pruning algorithm
	out 							= [y,id,firmid_orig,lagfirmid_orig,firmid,lagfirmid,cara,full(controls)];
	s								= ['../build/src/RAKIM/FIRST_PRUNED_' STRINGA '.csv']
	dlmwrite(s, out, 'precision', 16); 		
			
%Find the identified sets	
	identified_firmid				= identify_rakim_levels(id,firmid,lagfirmid,controls,1);
	identified_lagfirmid			= identify_rakim_levels(id,firmid,lagfirmid,controls,2);
	
%Tell me which person-job observation is identified
	identified_con					= identified_firmid(firmid);
	identified_lag					= identified_lagfirmid(lagfirmid);
	identified						= (identified_con==1 & identified_lag == 1);
	disp('share of identified person-job observations')   	
	mean(identified)
		
%Restrict the sample to job-year observations that are identified
	y								= y(identified);
	id								= id(identified);
	firmid							= firmid(identified);
	lagfirmid						= lagfirmid(identified);
	id_orig							= id_orig(identified);
	firmid_orig						= firmid_orig(identified);
	lagfirmid_orig					= lagfirmid_orig(identified);
	controls						= controls(identified,:);
	cara							= cara(identified,:);
	identified						= identified(identified);
	
	
%Re-normalize a bunch of things
	[~,~,id]						= unique(id);
	[~,~,firmid]					= unique(firmid);
	[~,~,lagfirmid]					= unique(lagfirmid);
    
%All leverages appears <1 from JLL algorithm so we do not run again initial algorithm.     
	
%Save the results of this first step
	out 					= [y,id,firmid_orig,lagfirmid_orig,firmid,lagfirmid,cara,full(controls)];
	s						= ['../build/src/RAKIM/PRUNED_' STRINGA '.csv']
	dlmwrite(s, out, 'precision', 16); 	
	
	end	
