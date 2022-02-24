function [pred_i pred_j] = prediction_DOUBLE_LEAVE(TYPE_TAXONOMY,taxonomy,TAXONOMY_I,TAXONOMY_J,MOVERS_i,MOVERS_j,FOCUS_I, FOCUS_J,Ppath_PROBLEM_i,Ppath_PROBLEM_j,y);



%Set up
prediction				= zeros(size(taxonomy,1),2);	
Npaths					= size(taxonomy,1); 

for ppp = 1:2
		
		if ppp 			== 1
		sel				= TAXONOMY_I == TYPE_TAXONOMY ;
		index_movers	= MOVERS_i(sel);
		Ppath			= Ppath_PROBLEM_i(sel); %1 or -1 to be used 
		FOCUS			= FOCUS_I(sel);
		end
		
		if ppp 			== 2
		sel				= TAXONOMY_J == TYPE_TAXONOMY ;
		index_movers	= MOVERS_j(sel);
		Ppath			= Ppath_PROBLEM_j(sel); %1 or -1 to be used 
		FOCUS			= FOCUS_J(sel);
		end
		
		pred			= y(index_movers).*Ppath;
		[aa,bb,cc]		= unique(FOCUS);

		
		pred 			= splitapply(@sum,pred,cc); %sum within i, this is going to be the predictor
				
		if ppp 			== 1				
		pred_i			= sparse(aa,1,pred,Npaths,1);
		end	
		
		if ppp 			== 2				
		pred_j			= sparse(aa,1,pred,Npaths,1);
		end				
end

end