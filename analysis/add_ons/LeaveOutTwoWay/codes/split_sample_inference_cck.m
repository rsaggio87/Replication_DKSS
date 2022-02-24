function 	[sigma_1 sigma_2 D adjustment]= split_sample_inference_cck(id_movers,firmid_delta,firmid_delta_f,ydelta,vec,index_C,Ndelta,GROUP_delta,group_stack,firmid_stack,id_stack,firmid_delta_f_stack,firmid_delta_stack,Pii_unrestrict,Pii_restrict);

list			 = [firmid_delta_f  firmid_delta id_movers ydelta];
list_aux 		 = [firmid_delta firmid_delta_f id_movers ydelta];

sel				 = (firmid_delta~=firmid_delta_f); %movers belonging to a given group
id_movers		 = id_movers(sel);
GRANDE_MEDIA 	 = mean(ydelta(sel));
vec				 = vec(sel);
index_C			 = index_C(sel);
list 			 = list(sel,:);
list_base		 = list;
list_aux 		 = list_aux(sel,:);
sel   			 = list(:,1) < list(:,2);
list(sel,:) 	 = list_aux(sel,:) ; %sort
[e,~,edge_index] = unique(list(:,1:2),'rows','stable');
list 			 = [list edge_index];
summary			 = accumarray(edge_index,1);
single_edge 	 = summary==1;
single_edge      = single_edge(edge_index);
list 			 = [list_base edge_index single_edge index_C]; %Vertex 1, Vertex 2, ID of the Movers, Ydelta, ID of the Move, Single Mover Indicator, Indexes in the original sample space

%build unweighted adj matrix and report the edge identifier in each non zero element of the Adj matrix. 
J 				 = max([firmid_delta;firmid_delta_f]);
jj				 = max(edge_index);
A_s 			 = sparse([e(:,1);e(:,2)],[e(:,2);e(:,1)],[(1:jj)';(1:jj)'],J,J);

%Fill vectors
Nmovers 		 = size(list,1);
list_path 		 = cell(Nmovers,1);
list_path_2 	 = cell(Nmovers,1);
yhat1			 = zeros(Nmovers,1);
yhat2			 = zeros(Nmovers,1);
list_edges		 = cell(Nmovers,1);
auxiliary		 = cell(Nmovers,1);
auxiliary_2		 = cell(Nmovers,1);
focus_m			 = cell(Nmovers,1);
index_i			 = cell(Nmovers,1);
index_j			 = cell(Nmovers,1);
P_i1			 = cell(Nmovers,1);
P_j2			 = cell(Nmovers,1);
Ppath_1			 = cell(Nmovers,1);
Ppath_2			 = cell(Nmovers,1);
no_A2			 = zeros(Nmovers,1);
ind_paths	 	 = ones(Nmovers,1);

%Now we enter the important part of the code. Here we will try to form two predictors of ydelta_i after leaving i out.
%We will use a shortest path algorithm to find the first prediction, which always exists because the input data is from a leave out connected set.
%The second prediction might exist or not exist. We'll record whether that is the case.

%Along the way we will also be filling in the elements of the auxiliary matrix D which will be useful for the next part

%Now we enter the important part of the code. Here we will try to form two predictors of ydelta_i after leaving i out.
%We will use a shortest path algorithm to find the first prediction, which always exists because the input data is from a leave out connected set.
%The second prediction might exist or not exist. We'll record whether that is the case.

%Along the way we will also be filling in the elements of the auxiliary matrix D which will be useful for the next part

tic
parfor i = 1:Nmovers
		
		%Find shortest leave one out path, this will be our first estimator. 
		first_run 		= 1;
		[yhat1(i),P1,list_1,~,auxiliary{i}] = shortest_path_double_leave_out(first_run,list(i,1),list(i,2),list(i,3),list,A_s);
		
		list_path{i} 	= list_1;
		Ppath_1{i}		= P1';
		type_of_mover(i)= 1;
		
		%Second estimator of yhat, after leaving indexes associated with first estimator.  
		first_run 		= 0;
		[yhat2(i),P2,list_2,~,auxiliary_2{i}] = shortest_path_double_leave_out(first_run,list(i,1),list(i,2),[list(i,3);list_1],list,A_s);
		
		list_path_2{i} 	= list_2;
		Ppath_2{i}		= P2';
		
		%Normalize things if A_2 is empty
		if isnan(yhat2(i)) == 1 
			no_A2(i)		= 1;
			yhat2(i)		= NaN;
			list_path_2{i}	= list(i,3);
			auxiliary_2{i}	= list(i,3);
			focus_m{i}		= [];
			index_i{i}		= [];
			index_j{i}		= [];
			P_i1{i}			= [];
			P_j2{i}			= [];
			Ppath_2{i}		= 1;	
		end
		
		%Update the first estimator
		if isnan(yhat2(i)) == 0
		
		tentativo_fallito = 0;
		ind_paths(i)= 2;
		conto=1;
		Nsize=1;
		
		
		while tentativo_fallito < 1  & ind_paths(i)<=4 &  Nsize<100 
		
				%update the prediction in subsample 1
				[tentativo,~,list_path_aux_tent] = shortest_path_double_leave_out(first_run,list(i,1),list(i,2),[list(i,3);list_1;list_2],list,A_s);
				
				if isnan(tentativo) == 0 
					list_1			= [list_1;list_path_aux_tent];
					list_1			= sort(list_1);
					[yhat1(i) P1]	= subsample_akm(i,list,list_1); 
					list_path{i} 	= list_1;
					Ppath_1{i}		= P1';
					type_of_mover(i)= 3;
					auxiliary{i}	= list(i,3).*(ones(size(list_1,1),1));
					ind_paths(i)    = ind_paths(i) +  1;
					Nsize			= min([size(list_1,1);size(list_2,1)]);
				end
				
				if isnan(tentativo) == 1 
					tentativo_fallito=1;
				end	
				
				
				%update the prediction in subsample 2
				[tentativo,~,list_path_aux_tent] = shortest_path_double_leave_out(first_run,list(i,1),list(i,2),[list(i,3);list_1;list_2],list,A_s);
				
				if isnan(tentativo) == 0 
					list_2			= [list_2;list_path_aux_tent];
					list_2			= sort(list_2);
					[yhat2(i) P2]	= subsample_akm(i,list,list_2); 	
					type_of_mover(i)= 4;
					list_path_2{i} 	= list_2;
					Ppath_2{i}		= P2';
					auxiliary_2{i}	= list(i,3).*(ones(size(list_2,1),1));
					ind_paths(i)    = ind_paths(i) +  1;
					Nsize			= min([size(list_1,1);size(list_2,1)]);
				end
				
				if isnan(tentativo) == 1
					tentativo_fallito=1;
				end
								
		
		conto=conto+1;
		end
		end
		
		%Reshape the indexes, this is needed for calculations of the matrix D.
		if 	no_A2(i) 		== 0
			index_j_C1 		= list_1;
			index_j_C2 		= list_2; 
			
			P1		   		= P1';
			P2		   		= P2';
			index_P1		= (1:size(index_j_C1,1))';
			index_P2		= (1:size(index_j_C2,1))';
			
		
			%Form all the pairwise combinations across indices that were used to create the two estimators 
			[A,B] 			= meshgrid(index_j_C1,index_j_C2);
			c				= cat(2,A',B');
			index_C			= reshape(c,[],2); %this tells me which elements of D I need to upgrade (so I have to look in the Ndelta space)
			
			
			[A,B] 			= meshgrid(index_P1,index_P2);
			c				= cat(2,A',B');
			index_P			= reshape(c,[],2); %this essentially reshapes the elements of C and P can be multiplied together
			
			
			%Close this animal.
			size_list		= size(index_C,1);
			focus_m{i}		= ones(size_list,1).*list(i,3);
			index_i{i}		= index_C(:,1);
			index_j{i}		= index_C(:,2);
			P_i1{i}			= P1(index_P(:,1));
			P_j2{i}			= P2(index_P(:,2));
		
		end

					
end
disp('Time to find leave one out path estimators')
toc
disp('Tabulation on # of independent paths found for a given i')
tabulate(ind_paths)


%Summarize what we have learned thus far
sel   = no_A2==1;
disp('share of workers for which is not possible to form two non-overlapping estimators of Delta y_i')
Q_share = mean(sel)


%Build \tilde{\sigma}
sigma_1			= (list(:,4)-yhat1).*list(:,4);
sigma_2 		= (list(:,4)-yhat1).*(list(:,4)-yhat2);
sigma_2(sel)	= (list(sel,4)-GRANDE_MEDIA).^2;

%Update \sigma_i for dudes for which we cannot find a second estimator
sigma_i(sel)	=  (list(sel,4)-GRANDE_MEDIA).^2;
%sigma_i(sel)	=  (list(sel,4)).^2;

%Reshuffle the outcomes in correct space
y				= sparse(id_movers,1,list(:,4),Ndelta,1);
vec				= sparse(id_movers,1,vec,Ndelta,1); 	%reshape in full Ndelta space (counting stayers as well);
sigma_1			= sparse(id_movers,1,sigma_1,Ndelta,1); %reshape in full Ndelta space (counting stayers as well);
sigma_2			= sparse(id_movers,1,sigma_2,Ndelta,1); %reshape in full Ndelta space (counting stayers as well);
no_A2			= sparse(id_movers,1,no_A2,Ndelta,1);   %reshape in full Ndelta space (counting stayers as well);
yhat1			= sparse(id_movers,1,yhat1,Ndelta,1);   %reshape in full Ndelta space (counting stayers as well);
yhat2			= sparse(id_movers,1,yhat2,Ndelta,1);   %reshape in full Ndelta space (counting stayers as well); 

%Save key stuff
focus_m 		= vertcat(focus_m{:,1});
index_i 		= vertcat(index_i{:,1});
index_j 		= vertcat(index_j{:,1});
P_i1 			= vertcat(P_i1{:,1});
P_j2 			= vertcat(P_j2{:,1});
Ppath_1 		= vertcat(Ppath_1{:,1});
Ppath_2 		= vertcat(Ppath_2{:,1});


%Now get cross products of the C matrix. The idea is to first compute C*b_i where b_i is the basis vector in position i. Once that is done, we then query the resulting nx1 vector according to individuals that built the prediction for i.
tic
mega_list		= [focus_m index_i index_j];
NMEGA			= size(mega_list,1);
C_vec1			= cell(Nmovers,1);
C_vec2			= cell(Nmovers,1);
parfor i = 1:Nmovers
	   focal_worker  = list(i,3)
	   focal_worker	 = sparse(focal_worker,1,1,Ndelta,1);
	   C_vec		 = C_build_cck(focal_worker,GROUP_delta,group_stack,firmid_stack,id_stack,firmid_delta_f_stack,firmid_delta_stack,Pii_unrestrict,Pii_restrict);
	   sel			 = (mega_list(:,1)==list(i,3));
	   my_pals	 	 = mega_list(sel,2)
	   C_vec1{i}	 = C_vec(my_pals);
	   my_pals	 	 = mega_list(sel,3)
	   C_vec2{i}	 = C_vec(my_pals);	 	
end
C_vec1		   = vertcat(C_vec1{:,1});
C_vec2		   = vertcat(C_vec2{:,1});
numero		   = C_vec1.*C_vec2.*P_i1.*P_j2;
D			   = sparse(index_i,index_j,numero,Ndelta,Ndelta);
disp('Time to build the matrix D')
toc


%Reshape things into a new world where for each i we provide the list of indexes in A_i1. Similarly for A_i2
index_movers 	= vertcat(list_path{:,1});    %this is A_i1.
index_movers_i 	= vertcat(auxiliary{:,1}); 

index_movers_2	= vertcat(list_path_2{:,1});  %this is A_i2.
index_movers_i_2= vertcat(auxiliary_2{:,1}); 


%Now find the problematic pairs.
paths_problem_i = cell(Nmovers,1);
paths_problem_j = cell(Nmovers,1);
A_i2_in_Aj1		= cell(Nmovers,1);
A_i1_in_Aj2		= cell(Nmovers,1);
A_i2_in_Aj2		= cell(Nmovers,1);

tic
%Now find the problematic pairs.
paths_problem_i = cell(Nmovers,1);
paths_problem_j = cell(Nmovers,1);
A_i2_in_Aj1		= cell(Nmovers,1);
A_i1_in_Aj2		= cell(Nmovers,1);
A_i2_in_Aj2		= cell(Nmovers,1);
mygroupP 		= @(x)x(randsample(size(x,1),1));

tic
parfor i = 1:Nmovers
	   sel			 	= ismember(index_movers_i_2,list(i,3));
	   A_i2	 	    	= [index_movers_2(sel)]; 
	
	   sel			 	= ismember(index_movers_i,list(i,3));
	   A_i1	 		 	= [index_movers(sel)]; %focus worker and all the movers used to derive the first estimator of y for the focus worker.
	   others_path	 	= [index_movers(~sel)];
	   others_index	 	= [index_movers_i(~sel)];
	   
	   problem1 		= ismember(others_path,A_i1);   	%a given path in A_ell,1 overlaps with A_i,1
	   problem2 		= ismember(others_index,A_i1);  	%ell is in A_i,
	   problem3 		= ismember(others_path,list(i,3));  %i is in A_ell,1
		
	   sel 				= (problem1+problem2+problem3>0);
	   
	   Nsel				= sum(sel);
	   list_prob		= [list(i,3).*ones(Nsel,1) others_index(sel)];
	   
	   %Output the problematic pairs
	   list_prob    	 = unique(list_prob,'rows','stable') %suppose A_i1={j,3,5,6} and A_j1 = {i,3,5,6}. Then (i,j) are clearly problematic and they will have "multiple" problems, as many paths in A_j1 overlap with A_i1. I don't need to keep track of all the (3,3), (5,5), etc, so just tell me that (i,j) are indeed problematic. Also double counting because of the normalization impose above (if i is in Aj1 and j is in Ai1 then the corresponding will have two rows with i,l where i<l)
	   paths_problem_i{i}= list_prob(:,1);
	   paths_problem_j{i}= list_prob(:,2);	    	
	
	   %Check overlapping paths
	   intersect=zeros(size(list_prob,1),3);
	  
	   %Second line: does A_i2 intersect A_j1?
	   sel				= ismember(index_movers_i,list_prob(:,2));
	   A_j1_i			= index_movers_i(sel);
	   A_j1				= [index_movers(sel)];
	   sel				= ismember(A_j1,A_i2);
	   [~,~,indici]		= unique(A_j1_i);
	   intersect(:,1)	= splitapply(@max,sel,indici);
	  		
	   %Third line: does A_i1 intersect A_j2?
	   sel				= ismember(index_movers_i_2,list_prob(:,2));
	   A_j2_i			= index_movers_i_2(sel);
	   A_j2				= [index_movers_2(sel)];
	   sel				= ismember(A_j2,A_i1);
	   [~,~,indici]		= unique(A_j2_i);
	   intersect(:,2)	= splitapply(@max,sel,indici);
	  		
	   %Fourth line: does A_i2 intersect A_j2?
	   sel				= ismember(A_j2,A_i2);
	   intersect(:,3)	= splitapply(@max,sel,indici);		
	   
	   
	   %Collect the results	
	   A_i2_in_Aj1{i}			= intersect(:,1);		
	   A_i1_in_Aj2{i}			= intersect(:,2);
	   A_i2_in_Aj2{i}			= intersect(:,3);
	   
end
disp('Time to find the problematic pairs')
toc
prob_pairs 				= [vertcat(paths_problem_i{:,1}) vertcat(paths_problem_j{:,1})];
PR						= sparse(prob_pairs(:,1),prob_pairs(:,2),1,max(list(:,3)),max(list(:,3)));
disp('checking, must report one')
issymmetric(PR)

%Bring important information in the space of problematic pairs
y_ii					= y(prob_pairs(:,1));
y_jj					= y(prob_pairs(:,2));

sigma_1_i				= sigma_1(prob_pairs(:,1));
sigma_1_j				= sigma_1(prob_pairs(:,2));

sigma_2_i				= sigma_2(prob_pairs(:,1));
sigma_2_j				= sigma_2(prob_pairs(:,2));

yhat_1_ii				= yhat1(prob_pairs(:,1));
yhat_1_jj				= yhat1(prob_pairs(:,2));
yhat_2_ii				= yhat2(prob_pairs(:,1));
yhat_2_jj				= yhat2(prob_pairs(:,2));

no_A2_ii				= no_A2(prob_pairs(:,1));
no_A2_jj				= no_A2(prob_pairs(:,2));



%Pin down C_ij of the problematic pairs (only need to look at upper triangular part)
[r c ]					= find(triu(PR'));
Cij						= cell(Nmovers,1);
tic
parfor i = 1:Nmovers
	   focal_worker  = list(i,3)
	   focal_worker	 = sparse(focal_worker,1,1,Ndelta,1);
	   C_vec		 = C_build_cck(focal_worker,GROUP_delta,group_stack,firmid_stack,id_stack,firmid_delta_f_stack,firmid_delta_stack,Pii_unrestrict,Pii_restrict);
	   sel			 = (r==list(i,3));
	   my_pals	 	 = c(sel,1);
	   Cij{i}=C_vec(my_pals);
end
disp('Time to build the Cij of problematic pairs')
toc
Cij		    			= vertcat(Cij{:,1});
Cij						= sparse(r,c,Cij,Ndelta,Ndelta); %make it a sparse matrix.
Cij						= Cij'+triu(Cij,1); %make it symmetric
Cij						= Cij(sub2ind(size(Cij),prob_pairs(:,1),prob_pairs(:,2))); %back to be a vector but now properly indexed
magical_number			= Cij.^2 +  2*(D(sub2ind(size(D),prob_pairs(:,1),prob_pairs(:,2)))) + 2*(D(sub2ind(size(D),prob_pairs(:,2),prob_pairs(:,1))))  ;

%Now figure out the taxonomy of the problematic pair: Start by looking whether \ell is in A_i1 or A_i2
j_in_Ai1				= ismember(prob_pairs,[index_movers_i index_movers],'rows');
j_in_Ai2				= ismember(prob_pairs,[index_movers_i_2 index_movers_2],'rows');

tabulate(j_in_Ai1+j_in_Ai2)

%Now whether i is in A_j1 or Aj2
i_in_Aj1				= ismember([prob_pairs(:,2) prob_pairs(:,1)],[index_movers_i index_movers],'rows');
i_in_Aj2				= ismember([prob_pairs(:,2) prob_pairs(:,1)],[index_movers_i_2 index_movers_2],'rows');
tabulate(i_in_Aj1+i_in_Aj2)

%The intersection in paths now
A_i2_in_Aj1 			= vertcat(A_i2_in_Aj1{:,1});
A_i1_in_Aj2 			= vertcat(A_i1_in_Aj2{:,1});
A_i2_in_Aj2 			= vertcat(A_i2_in_Aj2{:,1});

%Classify the problematic pairs now
taxonomy 				= 1.*(no_A2_ii== 0 & no_A2_jj == 0) + 2.*(no_A2_ii==1 | no_A2_jj== 1);
tabulate(taxonomy)


%Now i need to correction factor for this particular cut. Begin by finding the sigma_product
sigma_3_i				= y_ii.*(y_ii-yhat_2_ii);
sigma_3_j				= y_jj.*(y_jj-yhat_2_jj);

sigma_4_i				= (y_ii-GRANDE_MEDIA).^2;
sigma_4_j				= (y_jj-GRANDE_MEDIA).^2;

%Compute the adjusted sigma_product now
[sigma_product bad_pairs tag] = compute_sigma_product(taxonomy,no_A2_jj,no_A2_ii,j_in_Ai1,j_in_Ai2,i_in_Aj1,i_in_Aj2,A_i2_in_Aj1,A_i1_in_Aj2,A_i2_in_Aj2,sigma_1_i,sigma_2_i,sigma_3_i,sigma_4_i,sigma_1_j,sigma_2_j,sigma_3_j,sigma_4_j,magical_number);


%Adjustment factor
adjustment				= 2*sum((magical_number).*(sigma_product-(sigma_1_i.*sigma_1_j))); 


end




    