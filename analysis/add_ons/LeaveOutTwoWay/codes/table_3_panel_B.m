function TABELLA=table_3_panel_B(filename_old,filename_young)
	
	
	%Start by loading old workers files
	results						= importdata([filename_old '.csv']);
	y_old						= results(:,1);
	firmid_old					= results(:,2);
	id_old						= results(:,3);
	old							= results(:,7);
	firmid_orig_old				= results(:,4);
	id_orig_old					= results(:,5);
	size_old					= accumarray(firmid_old,1);
	size_old					= size_old(firmid_old);
	clear results
	
	
	%Now load young workers
	results                     = importdata([filename_young '.csv']);
	y_young						= results(:,1);
	firmid_young				= results(:,2);
	id_young					= results(:,3);
	young						= results(:,7);
	firmid_orig_young			= results(:,4);
	id_orig_young				= results(:,5);
	size_young					= accumarray(firmid_young,1);
	size_young					= size_young(firmid_young);
	clear results
	
	%Find overlapping firms  
	sel    						= ismember(firmid_orig_old,firmid_orig_young);
	sel_young    				= ismember(firmid_orig_young,firmid_orig_old);
	
	%Summarize Dual Connected set of firms
	id							= [id_orig_old(sel); id_orig_young(sel_young)];
	[~,~,id]					= unique(id);
	firmid_orig					= [firmid_orig_old(sel); firmid_orig_young(sel_young)];
	[~,~,firmid]				= unique(firmid_orig);
	y							= [y_old(sel); y_young(sel_young)];
	old_dummy					= [old(sel); young(sel_young)];
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
	s=['Share young: ' num2str(mean(old_dummy))];
	disp(s);
	
	%Find largest firm in dual connected set
	key_firm					= find(SIZE==max(SIZE));
	key_firm					= mean((firmid_orig(key_firm)))
	
	%Renormalize the set of firm effects for old according to key_firm 
	to_norm						= mean(firmid_old(firmid_orig_old==key_firm));
	firmid_old 					= normalize_firm_effects(firmid_old,to_norm);
	
	%Renormalize the set of firm effects for young according to key_firm
	to_norm						= mean(firmid_old(firmid_orig_young==key_firm));
	firmid_young				= normalize_firm_effects(firmid_young,to_norm);
	
	%Estimate the two sets of firm effects again, save a .csv file that will be read by Stata.
	filename 					=['results/old_estimated_renormalized'];
	[~] 						= quick_AKM(y_old,id_old,firmid_old,id_orig_old,firmid_orig_old,filename);
	filename 					=['results/young_estimated_renormalized'];
	[~] 						= quick_AKM(y_young,id_young,firmid_young,id_orig_young,firmid_orig_young,filename);

end


