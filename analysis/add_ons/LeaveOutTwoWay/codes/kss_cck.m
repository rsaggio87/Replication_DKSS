function [ttest,V, student]= kss_cck(GROUP,y,id,firmid,year,controls, leave_out_level, type_of_algorithm,epsilon,filename)
	
	%NOTE: Currently works for T=2 and if group variable does not vary within individual. 
	
	%Residualize w.r.t to controls in full model
	if size(controls,2)>0
	y = resid_AKM(y,id,firmid,controls);
	end

	%Find leave one/two leave out CS for GROUP 0, run KSS
	filename0											= [filename 'GROUP_0']
	sel													= (GROUP == 0);
	[one]												= leave_out_COMPLETE(y(sel),id(sel),firmid(sel),leave_out_level,year(sel),GROUP(sel,:),2,0,0,0,0,0, type_of_algorithm,epsilon,filename0);
			 	
	%Find leave one/two leave out CS for GROUP 1
	filename1											= [filename 'GROUP_1']
	sel													= (GROUP == 1);
	[one]												= leave_out_COMPLETE(y(sel),id(sel),firmid(sel),leave_out_level,year(sel),GROUP(sel,:),2,0,0,0,0,0, type_of_algorithm,epsilon,filename1);
	
	%Find the dual connected set, normalize firm effects, create dataset used to plot firm effects group 1 against firm effects in group 0.  
	results_0											= importdata([filename 'GROUP_0.csv']);
	results_1											= importdata([filename 'GROUP_1.csv']);
	[key_firm,r]										= dual_connected_set(filename,results_0,results_1,leave_out_level,type_of_algorithm,epsilon);	

	%Now run test for equal firm effects
	[ttest,V, student]									= testing_cck(key_firm,r,filename,results_0,results_1,leave_out_level,type_of_algorithm,epsilon)

	
	
end


