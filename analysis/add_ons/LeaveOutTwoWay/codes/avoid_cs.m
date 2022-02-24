function [share_paths] = avoid_cs(list,list_path,auxiliary,yhat,eta_h,Pii,A_s)

	SOMMA 			= 0;
	SOMMA2 			= 0;
	Nmovers 		= size(list,1);
	y				= list(:,4);
	sigma_hat 		= y.*eta_h;
	dashboard_1		= y.*(y-yhat);
	dashboard_2		= (1/12)*y.*(sigma_hat+3*y.*(eta_h-y));
	movers_path 	= vertcat(list_path{:,1});
	index_movers 	= vertcat(auxiliary{:,1});
	share_paths		= zeros(Nmovers,1);
	
	for i = 1 : Nmovers
		
			
			%do the products
			if 0 == 1 
			coMovers = ismember(list(:,5),list(i,5)); %sharing a move?
			d2		 = dashboard_2(i)+dashboard_2+0.5*y(i)^2*y.^2;
			d1		 = dashboard_1(i).*dashboard_1;
			d3		 = y(i)^2*y.^2;
			end
			
			%now tell me where the products are associated with non-overlapping samples
			focus_list	  						= list_path{i};
			sel									= ismember(movers_path,focus_list);
			dudes_to_find_new_pred				= index_movers(sel);
			i_index								= dudes_to_find_new_pred==list(i,3);
			dudes_to_find_new_pred(i_index) 	= []; % leaving out i
			dudes_to_find_new_pred				= unique(dudes_to_find_new_pred);  %if more than the move overlap it will be counting twice, we don't need that) %if more than the move overlap it will be counting twice, we don't need that
			same_path							= ismember(list(:,3),dudes_to_find_new_pred);
			share_paths(i)						= mean(same_path);
			
			
			if 0 == 1
			%Update the sum, for the moment ignoring C
			not_i			= ~ismember(list(:,3),list(i,3));
			SOMMA 			= SOMMA  + sum(d2.*coMovers.*not_i) + sum(d1.*(1-coMovers).*(1-same_path).*not_i) + sum(d3.*(1-coMovers).*(same_path).*not_i);
			SOMMA2			= SOMMA2 + sum(d2.*coMovers.*not_i) + sum(d1.*(1-coMovers).*not_i);
			end
	end

end

    