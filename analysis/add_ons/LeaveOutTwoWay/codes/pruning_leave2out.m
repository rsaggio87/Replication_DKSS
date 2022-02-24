function [id,firmid,y,controls,id_old,firmid_old] = pruning_leave2out(id,firmid,y,controls,id_old,firmid_old)

	tot_art_point=1	
	while tot_art_point>0
		list 			= list_movers(firmid,id); %this gives me the list of movers and associated firms. list is nmovers x 3. First element is identifier of the mover, third element is the original one. 
		n_movers		= max(list(:,1));
		bad 			= zeros(n_movers,1);
		bad_workers 	= zeros(n_movers,1);
		artic_points 	= cell(n_movers,1);
	
		parfor 			i = 1:n_movers
						[n_of_bad_workers,bad_workers(i), artic_points{i}] = leave2out(i,list);
						bad(i)=n_of_bad_workers>0;
		end
						tot_art_point = sum(bad);
						disp('Total Articulation points found across Leave one Out graphs')
						artic_points=vertcat(artic_points{:,1});
						size(unique(artic_points),1)

		sel 			= ~ismember(id,artic_points);
		id				= id(sel);
		firmid			= firmid(sel);
		y				= y(sel);
		id_old			= id_old(sel);
		firmid_old		= firmid_old(sel);
		controls 		= controls(sel,:);
		
		%Reset id
		[~,~,n]			= unique(firmid);
		firmid			= n;
		[~,~,n]			= unique(id);
		id				= n; 
		
		%Leave out connected set after leaving the union of cut vertexes for leave out graphs
		[y,firmid,id,id_old,firmid_old,controls] = pruning_unbal_v3(y,firmid,id,id_old,firmid_old,controls);
						
	end
    
end    