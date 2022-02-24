function [share, C_share] = paths_network_double_leave(firmid,id,id_movers,firmid_delta,firmid_delta_f,ydelta,eta_h,Pii,L,Lchol,Fdelta,I_Lambda_P,L_P,Lambda_B,F,C)
% Author: Raffaele Saggio. raffaele.saggio@berkeley.edu

%  This function is an attempt to construct an unbiased estimator of the standard error
%  for the variance of firm effects. 

TOT_C = sum(sum(C.^2,1));


%Construct the list of movers, tell me who is a unique mover
vec 			 = C*ydelta;
list			 = [firmid_delta_f  firmid_delta id_movers ydelta];
list_aux 		 = [firmid_delta firmid_delta_f id_movers ydelta];
sel   			 = firmid_delta~=firmid_delta_f;
eta_h			 = eta_h(sel);
vec 			 = vec(sel);
list 			 = list(sel,:);
list_base		 = list;
list_aux 		 = list_aux(sel,:);
sel   			 = list(:,1) < list(:,2);
list(sel,:) 	 = list_aux(sel,:) ;
[e,~,edge_index] = unique(list(:,1:2),'rows','stable');
list 			 = [list edge_index];
summary			 = accumarray(edge_index,1);
single_edge 	 = summary==1;
single_edge      = single_edge(edge_index);
list 			 = [list_base edge_index single_edge]; %Vertex 1, Vertex 2, ID of the Movers, Ydelta, ID of the Move, Single Mover Indicator. 
sigma_hat		 = list(:,4).*eta_h;

%build unweighted adj matrix and the one of indexes
J 				 = max(firmid);
jj				 = max(edge_index);
A_s 			 = sparse([e(:,1);e(:,2)],[e(:,2);e(:,1)],[(1:jj)';(1:jj)'],J,J);


%Now tell me the shortest path and predictors to reach the firms of a given mover when leaving that mover out of the network.
Nmovers 		 = size(list,1);
yhat1			 = zeros(Nmovers,1);
yhat2			 = zeros(Nmovers,1);

tic
parfor i = 1:Nmovers
		
		%First estimator
		[yhat1(i),list_path,list_edges] = shortest_path_double_leave_out(list(i,1),list(i,2),list(i,3),list,A_s);
		
		%Unique mover within edge
		ids_to_drop = accumarray(list_edges, list_path, [], @min);
		
		%Second estimator, after dropping movers belonging to the path of the first estimator
		[yhat2(i)] = shortest_path_double_leave_out(list(i,1),list(i,2),[list(i,3);ids_to_drop],list,A_s);
		

end
disp('Time to find double leave one out path estimators')
toc

%Info on the Q-set.
problem = isnan(yhat2);
disp('share of workers for which is not possible to form two estimators of Delta y_i leaving y_i out')
share	= mean(problem)
sel 	= find(problem);
total 	= sum(vec.^2);
C_share = (sum(vec(sel).^2))/total

end

    