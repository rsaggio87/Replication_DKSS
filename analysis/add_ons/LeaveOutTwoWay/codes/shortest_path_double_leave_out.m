function [yhat,P,list_path,list_edges,auxiliary,auxiliary_edges] = shortest_path_double_leave_out(first_run,j,jprime,dudes_to_leave_out,list_orig,A_s);
			
			
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			%j is the firm in period 2 of individual i.
			
			%jprime is the firm in period 1 of individual i.
			
			%list_orig contains the list of vertex, movers indicators, outcomes, etc. See mother function.
			
			%dudes_leave_out contains the list of workers to leave out from graph. In position 1, it must always report individual i.
			
			%A_s Adj matrix that reports the indexes of the moves.
			
			
			%yhat: prediction of y_i leaving out the movers contained in the list dudes_to_leave_out.
			
			%list_path:  movers associated with the shortest path that leads me to j and jprime after leaving out dudes_to_leave_out
			
			%list_edges: edges associated with the shortest path that leads me to j and jprime after leaving out dudes_to_leave_out
			
			%auxiliary: identifier of i.
			
			%auxiliary_edges: edge of individual i.
			
			%index_i: indexes in the list space  of i.
			
			%index_j: indexes in the list space  for individuals forming the path to get to i.
		
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			
			
			J=size(A_s,1);
			
			
			%Leave out movers contained in dudes_to_leave_out and draw corresponding Adj Matrix.
			index_orig = (1:size(list_orig,1))';
			list_orig  = [list_orig index_orig];
			sel =~ismember(list_orig(:,3),dudes_to_leave_out);
			list = list_orig(sel,:);
			A 	= sparse([list(:,1);list(:,2)],[list(:,2);list(:,1)],1,J,J);
			[r,c,v] = find(A);
			A 	= sparse(r,c,1,J,J);
			
			%Find Shortest Path
			[dist, PATH ]=graphkshortestpaths( A, j, jprime, 1 );
		
			%Reconstruct the connected path and the ids associated with such path
			PATH 			= cell2mat(PATH);
			if size(PATH,1)>0 
				NPATHS  		= size(PATH,2)-1;	
				PATH 			= repelem(PATH,2);
				PATH			= PATH(2:end-1);
			
				if NPATHS > 1
						   PATH = reshape(PATH,[],NPATHS)'; 
				end
				
				%Give me back the edges identifiers associated with the path
				edges			= A_s(sub2ind(size(A_s),PATH(:,1),PATH(:,2)));
			
			
				%Save the edges and *all* the movers travelling these edges
				sel				= ismember(list(:,5),edges);
				list_edges		= list(sel,5);
				list_path		= list(sel,3);
				randB 			= @(x)x(randsample(size(x,1),1));
				[~,~,list_edges] = unique(list_edges);
				
				%Now that we got the edges give me a unique mover transitioning over that edge
				if first_run == 1
				list_path 		= splitapply(randB,list_path,list_edges); 
				end
				list_path		= sort(list_path);%need to resort for MC only
				
				
				%Now that we got the movers identifiers, save a bunch of very important things
				sel 			= ismember(list(:,3),list_path);
				list_edges		= list(sel,3);
				list_edges		= list(sel,5);
				index_j			= list(sel,end);
				
				%Complimentary information related to individual i.
				auxiliary		= dudes_to_leave_out(1,1).*ones(size(list_path,1),1); 
				index			= find(list_orig(:,3)==dudes_to_leave_out(1,1));
				auxiliary_edges = list_orig(index,5).*ones(size(list_path,1),1); 
				index_i			= list_orig(index,end).*ones(size(list_path,1),1); 
			
			
				%Build the predictor
				actual_path	    = [list(sel,1) list(sel,2)];
				NPATHS			= size(actual_path,1);
				[aa,bb,cc]		= unique(actual_path);
				cc				= reshape(cc,NPATHS,2);
				Nlist			= size(cc,1);
				Nodes			= max(max(cc));
				F				= sparse((1:Nlist)',cc(:,1),1,Nlist,Nodes) + sparse((1:Nlist)',cc(:,2),-1,Nlist,Nodes) ;
				outcomes		= list(sel,4);
				psi				= pinv(full(F'*F))*F'*outcomes;
				indexJ			= find(aa==j);
				indexJprime		= find(aa==jprime);
				yhat	 		= psi(indexJ) - psi(indexJprime);
				aux 			= sparse(1,indexJ,1,1,Nodes) + sparse(1,indexJprime,-1,1,Nodes) ;
				P				= aux*(pinv(full(F'*F))*F');
				
			end		
			
			if size(PATH,1)==0
				yhat 		= NaN;
				list_path	= NaN;
				list_edges  = NaN;
				P			= NaN;
				auxiliary	= NaN;
				auxiliary_edges = NaN;
				
			end
			
end


    