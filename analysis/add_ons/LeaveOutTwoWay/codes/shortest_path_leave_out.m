function [list_path yhat dist] = shortest_path_leave_out(list,A_s);

	Nmovers 	= size(list,1);
	list_path 	= cell(Nmovers,2);
	yhat		= zeros(Nmovers,1);
	J			= max(list(:,1));
	list_orig 	= list;
	
	
	for i = 1 : Nmovers
			j=list_orig(i,1);
			jprime=list_orig(i,2);
			sel =~ismember(list_orig(:,3),list_orig(i,3));
			list = list_orig(sel,:);
			A 	= sparse([list(:,1);list(:,2)],[list(:,2);list(:,1)],1,J,J);
			[r,c,v] = find(A);
			A 	= sparse(r,c,1,J,J);
			[ DIST, PATH ]=graphkshortestpaths( A, j, jprime, 1 );
			dist(i)=DIST;
		
			%Reconstruct the connected path and the ids associated with such path
			PATH 			= cell2mat(PATH); 
			NPATHS  		= size(PATH,2)-1;	
			PATH 			= repelem(PATH,2);
			PATH			= PATH(2:end-1);
			
			if NPATHS > 1
					   PATH = reshape(PATH,[],NPATHS)'; 
			end
			
			edges			= A_s(sub2ind(size(A_s),PATH(:,1),PATH(:,2)));
			
			sel				= ismember(list(:,5),edges);
			aux				= list(sel,3);
			list_path{i,1}	= aux;
			list_path{i,2}	= list_orig(i,3).*ones(size(aux,1),1);
		
			%Build the predictor
			actual_path	    = [list(sel,1) list(sel,2)];
			NPATHS			= size(actual_path,1);
			[aa,bb,cc]		= unique(actual_path);
			cc				= reshape(cc,NPATHS,2);
			Nlist			= size(cc,1);
			Nodes			= max(max(cc));
			F				= sparse((1:Nlist)',cc(:,1),1,Nlist,Nodes) + sparse((1:Nlist)',cc(:,2),-1,Nlist,Nodes) ;
			outcomes		= list(sel,4);
			allRowsEqual 	= all(ismember(F, F(1,:), 'rows'));
			psi				= pinv(full(F'*F))*F'*outcomes;
			yhat(i) 		= psi(find(aa==j)) - psi(find(aa==jprime));

			
				
	end

end

    