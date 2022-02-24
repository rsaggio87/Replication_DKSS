function  [V_sym V_nosym Q_share C_share_1 C_share_2 C_share_3 C_share_4 C_share_5 Dsym vec index_movers index_movers_i index_movers_2 index_movers_i_2 Ppath_1 Ppath_2 no_A2 list_problem taxonomy magical_number TAXONOMY_I TAXONOMY_J MOVERS_i MOVERS_j FOCUS_I FOCUS_J Ppath_PROBLEM_i Ppath_PROBLEM_j taxonomyALT TAXONOMY_ALT_I TAXONOMY_ALT_J MOVERS_ALT_i MOVERS_ALT_j FOCUS_ALT_I FOCUS_ALT_J Ppath_ALT_PROBLEM_i Ppath_ALT_PROBLEM_j] = paths_network_v2(id_movers,firmid_delta,firmid_delta_f,ydelta,vec,Fdelta,L,F,tol,epsilon,type_,pfun_,Lambda_P,Lambda_B)


%Construct the list of movers, tell me who is a unique mover
Ndelta			 = size(ydelta,1);
index_C			 = (1:Ndelta)';
list			 = [firmid_delta_f  firmid_delta id_movers ydelta];
list_aux 		 = [firmid_delta firmid_delta_f id_movers ydelta];
sel   			 = firmid_delta~=firmid_delta_f; %focus on movers only
id_movers		 = id_movers(sel);
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
tic
Nmovers 		 = size(list,1);
list_path 		 = cell(Nmovers,1);
list_path_2 	 = cell(Nmovers,1);
yhat1			 = zeros(Nmovers,1);
yhat2			 = zeros(Nmovers,1);
list_edges		 = cell(Nmovers,1);
auxiliary		 = cell(Nmovers,1);
auxiliary_2		 = cell(Nmovers,1);
auxiliary_edges	 = cell(Nmovers,1);
focus_m			 = cell(Nmovers,1);
index_i			 = cell(Nmovers,1);
index_j			 = cell(Nmovers,1);
P_i1			 = cell(Nmovers,1);
P_j2			 = cell(Nmovers,1);
Ppath_1			 = cell(Nmovers,1);
Ppath_2			 = cell(Nmovers,1);
no_A2			 = zeros(Nmovers,1);
toc

%Now we enter the important part of the code. Here we will try to form two predictors of ydelta_i after leaving i out.
%We will use a shortest path algorithm to find the first prediction, which always exists because the input data is from a leave out connected set.
%The second prediction might exist or not exist. We'll record whether that is the case.

%Along the way we will also be filling in the elements of the auxiliary matrix D which will be useful for the next part

tic
parfor i = 1:Nmovers
		
		%Find shortest leave one out path, this will be our first estimator. 
		[yhat1(i),P1,list_path_aux,list_edges_aux,auxiliary{i},auxiliary_edges{i}] = shortest_path_double_leave_out(list(i,1),list(i,2),list(i,3),list,A_s);
		
		list_path{i} 	= list_path_aux;
		list_edges{i}	= list_edges_aux;
		Ppath_1{i}		= P1';
		
		%Second estimator of yhat, after leaving indexes associated with first estimator.  
		[yhat2(i),P2,list_path_aux2,~,auxiliary_2{i}] = shortest_path_double_leave_out(list(i,1),list(i,2),[list(i,3);list_path_aux],list,A_s);
		
		list_path_2{i} 				 = list_path_aux2;
		
		%Normalize things if A_2 is empty
		if isnan(yhat2(i)) == 1 
			no_A2(i)		= 1;
			yhat2(i)		= 1;
			list_path_2{i}	= list(i,3);
			auxiliary_2{i}	= list(i,3);
			focus_m{i}		= [];
			index_i{i}		= [];
			index_j{i}		= [];
			P_i1{i}			= [];
			P_j2{i}			= [];
			Ppath_2{i}		= 1;	
		end
		
		
		%Reshape the indexes, this is needed for calculations of the matrix D.
		if 	no_A2(i) 		== 0
			index_j_C1 		= list_path_aux;
			index_j_C2 		= list_path_aux2;
			Ppath_2{i}		= P2'; 
			
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
			
			
			%Close this bitch.
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


%Summarize what we have learned thus far
sel   = no_A2==1;
disp('share of workers for which is not possible to form two non-overlapping estimators of Delta y_i')
Q_share = mean(sel)

%Build \tilde{\sigma}
sigma_1			= (list(:,4)-yhat1).*list(:,4);
sigma_2 		= (list(:,4)-yhat1).*(list(:,4)-yhat2);
sigma_2(sel)	= (list(sel,4)).^2;

%Reshuffle the outcomes in correct space
y					= sparse(id_movers,1,list(:,4),Ndelta,1);
vec				    = sparse(id_movers,1,vec,Ndelta,1); 	%reshape in full Ndelta space (counting stayers as well);
sigma_1				= sparse(id_movers,1,sigma_1,Ndelta,1); %reshape in full Ndelta space (counting stayers as well);
sigma_2				= sparse(id_movers,1,sigma_2,Ndelta,1); %reshape in full Ndelta space (counting stayers as well);
no_A2				= sparse(id_movers,1,no_A2,Ndelta,1);   %reshape in full Ndelta space (counting stayers as well);
yhat1				= sparse(id_movers,1,yhat1,Ndelta,1);   %reshape in full Ndelta space (counting stayers as well);
yhat2				= sparse(id_movers,1,yhat2,Ndelta,1);   %reshape in full Ndelta space (counting stayers as well);


%Save key stuff
focus_m 				= vertcat(focus_m{:,1});
index_i 				= vertcat(index_i{:,1});
index_j 				= vertcat(index_j{:,1});
P_i1 					= vertcat(P_i1{:,1});
P_j2 					= vertcat(P_j2{:,1});
Ppath_1 				= vertcat(Ppath_1{:,1}); %this will only be used for the Montecarlo
Ppath_2 				= vertcat(Ppath_2{:,1});

%Now get cross products of the C matrix
tic
P						= diag(Lambda_P);
B						= diag(Lambda_B);
mega_list				= [focus_m index_i; focus_m index_j];
NMEGA					= size(mega_list,1)
[mega_list,~,index_mega]= unique(mega_list,'rows','stable'); %individual cross products that I need to query.
elist					= [firmid_delta(mega_list(:,1)) firmid_delta_f(mega_list(:,1)) firmid_delta(mega_list(:,2)) firmid_delta_f(mega_list(:,2))];
[elist,~,index_unique]	= unique(elist,'rows');
[Mij, Bij] 				= eff_res_FAST_FE_ONLY(elist,Fdelta,L,F,tol,epsilon,type_,pfun_); %Apply function
Mij						= -Mij(index_unique);%Reshape into the original elist space where M=I-P;
Bij						= Bij(index_unique);%Reshape into the original elist space
Mii						= 1-P(mega_list(:,1)); %notice T=2 requirement here.
Mjj						= 1-P(mega_list(:,2));
Bii						= B(mega_list(:,1));
Bjj						= B(mega_list(:,2));
Cij						= Bij - 0.5*Mij.*(Mii.\Bii) -0.5*Mij.*(Mjj.\Bjj); %top of page 20
Cij						= Cij(index_mega);
Cij						= reshape(Cij,[NMEGA/2,2]);
numero					= Cij(:,1).*Cij(:,2).*P_i1.*P_j2;
D						= sparse(index_i,index_j,numero,Ndelta,Ndelta);
disp('Time to build the matrix D')
toc

%Reshape things into a new world where for each i we provide the list of indexes in A_i1. Similarly for A_i2
index_movers 	= vertcat(list_path{:,1});    %this is A_i1.
index_movers_i 	= vertcat(auxiliary{:,1}); 

index_movers_2	= vertcat(list_path_2{:,1});  %this is A_i2.
index_movers_i_2= vertcat(auxiliary_2{:,1}); 


%Now we find which pairs (i,l) are problematic.
paths_problem_i = cell(Nmovers,1);
paths_problem_j = cell(Nmovers,1);

parfor i = 1:Nmovers			%nota: questo parfor puo' essere chiaramente evitato con uno uso intellingente di meshgrid 
	   sel			 		= ismember(index_movers_i,list(i,3));
	   A_i1	 		 		= [index_movers(sel)]; %focus worker and all the movers used to derive the first estimator of y for the focus worker.
	   
	   others_path	 		= [index_movers(~sel)];
	   others_index	 		= [index_movers_i(~sel)];
	   
	   problem1 			= ismember(others_path,A_i1);   %a given path in A_ell,1 overlaps with A_i,1
	   problem2 			= ismember(others_index,A_i1);  %ell is in A_i,
	   problem3 			= ismember(others_path,list(i,3));  %i is in A_ell,1
		
	   sel 					= (problem1+problem2+problem3>0);
	   
	   Nsel				 	= sum(sel);
	   paths_problem_i{i}   = list(i,3).*ones(Nsel,1);
	   paths_problem_j{i}	= others_index(sel);	   
	   	    	
end


%Summarize what we have learned on these problematic pairs.
list_problem 	= [vertcat(paths_problem_i{:,1}) vertcat(paths_problem_j{:,1})];
list_problem    = unique(list_problem,'rows','stable'); %suppose A_i1={j,3,5,6} and A_j1 = {i,3,5,6}. Then (i,j) are clearly problematic and they will have "multiple" problems, as many paths in A_j1 overlap with A_i1. I don't need to keep track of all the (3,3), (5,5), etc, so just tell me that (i,j) are indeed problematic.

%Reshaping the list of problematic pairs to unique pairs. So only (i,j) are in \mathcal{P} iff i>l (lower triangular matrix) 
sel 			= (list_problem(:,1)>list_problem(:,2));
list_problem	= list_problem(sel,:);

%I know some of these calculations might be redundant, will improve in later version
elist					= [firmid_delta(list_problem(:,1)) firmid_delta_f(list_problem(:,1)) firmid_delta(list_problem(:,2)) firmid_delta_f(list_problem(:,2))];
[elist,~,index_unique]	= unique(elist,'rows');
[Mij, Bij] 				= eff_res_FAST_FE_ONLY(elist,Fdelta,L,F,tol,epsilon,type_,pfun_); %Apply function
Mij						= -Mij(index_unique);%Reshape into the original elist space where M=I-P;
Bij						= Bij(index_unique);%Reshape into the original elist space
Mii						= 1-P(list_problem(:,1)); %notice T=2 requirement here.
Mjj						= 1-P(list_problem(:,2));
Bii						= B(list_problem(:,1));
Bjj						= B(list_problem(:,2));
Cij						= Bij - 0.5*Mij.*(Mii.\Bii) -0.5*Mij.*(Mjj.\Bjj); %top of page 20


%Bring important information in the correct space now (the \mathcal{P} set)
y_ii			= y(list_problem(:,1));
y_jj			= y(list_problem(:,2));

sigma_1_i		= sigma_1(list_problem(:,1));
sigma_1_j		= sigma_1(list_problem(:,2));

sigma_2_i		= sigma_2(list_problem(:,1));
sigma_2_j		= sigma_2(list_problem(:,2));

no_A2_i			= no_A2(list_problem(:,1));
no_A2_j			= no_A2(list_problem(:,2));

yhat_1_ii		= yhat1(list_problem(:,1));
yhat_1_jj		= yhat1(list_problem(:,2));
yhat_2_ii		= yhat2(list_problem(:,1));
yhat_2_jj		= yhat2(list_problem(:,2));


%Key numbers that we want to calculate out of these problematic pairs
Npaths			= size(list_problem,1);
taxonomy		= -9.*ones(Npaths,1);
taxonomyALT		= -9.*ones(Npaths,1);
sigma_product	= zeros(Npaths,1);
sigma_product_ALT=zeros(Npaths,1);
magical_number 	 = Cij.^2 +  2*(D(sub2ind(size(D),list_problem(:,1),list_problem(:,2)))) + 2*(D(sub2ind(size(D),list_problem(:,2),list_problem(:,1))))  ;

Ppath_PROBLEM_i	 = cell(Npaths,1);
Ppath_PROBLEM_j	 = cell(Npaths,1);
MOVERS_i		 = cell(Npaths,1);
MOVERS_j		 = cell(Npaths,1);
TAXONOMY_I		 = cell(Npaths,1);
TAXONOMY_J		 = cell(Npaths,1);
FOCUS_I			 = cell(Npaths,1);
FOCUS_J			 = cell(Npaths,1);


Ppath_ALT_PROBLEM_i	 = cell(Npaths,1);
Ppath_ALT_PROBLEM_j	 = cell(Npaths,1);
MOVERS_ALT_i	 	 = cell(Npaths,1);
MOVERS_ALT_j	 	 = cell(Npaths,1);
TAXONOMY_ALT_I		 = cell(Npaths,1);
TAXONOMY_ALT_J		 = cell(Npaths,1);
FOCUS_ALT_I			 = cell(Npaths,1);
FOCUS_ALT_J			 = cell(Npaths,1);

PREDIZIONE_i	 = zeros(Npaths,1);
PREDIZIONE_j	 = zeros(Npaths,1);


tic
parfor i = 1:Npaths
	
	%Which individuals we are looking at
	focus_i 	= list_problem(i,1);
	focus_j 	= list_problem(i,2);
	
	%Their outcomes
	y_i			= y_ii(i);
	y_j			= y_jj(i);
	
	%Predictors for i
	yhat_1_i	= yhat_1_ii(i);
	yhat_2_i	= yhat_2_ii(i);
	
	%Predictors for j
	yhat_1_j	= yhat_1_jj(i);
	yhat_2_j	= yhat_2_jj(i);
	
	%Do they have an A2?
	i_no_A2	 	= no_A2_i(i);
	j_no_A2	 	= no_A2_j(i);
	
	%Is j in A_i1, A_i2?
	PATH_i		= index_movers(index_movers_i==focus_i);
	sel			= ismember(PATH_i,focus_j);
	j_in_Ai1	= max(sel);
	PATH_i		= index_movers_2(index_movers_i_2==focus_i);
	sel			= ismember(PATH_i,focus_j);
	j_in_Ai2	= max(sel);
	
	%Is i in A_j1,A_j2
	PATH_j		= index_movers(index_movers_i==focus_j);
	sel			= ismember(PATH_j,focus_i);
	i_in_Aj1	= max(sel);
	PATH_i		= index_movers_2(index_movers_i_2==focus_j);
	sel			= ismember(PATH_i,focus_i);
	i_in_Aj2	= max(sel);
	
	%Solve prediction problem for i
	index		 				= find(list(:,3)==focus_i); %%to find the set of firms associated with i
	[pred_i,PP_i_orig,patH_i]  	= shortest_path_double_leave_out(list(index,1),list(index,2),[focus_i;focus_j],list,A_s);
	
	
	%Solve prediction problem for j
	index		 				= find(list(:,3)==focus_j); %%to find the set of firms associated with j
	[pred_j,PP_j_orig,patH_j] 	= shortest_path_double_leave_out(list(index,1),list(index,2),[focus_i;focus_j],list,A_s);
	
	PP_i_orig					= PP_i_orig';
	PP_j_orig					= PP_j_orig';
	
	
	if isnan(pred_i) == 1 & isnan(pred_j) == 1 %%this pair breaks connectedness
			taxonomy(i) = 6;
			
			if magical_number(i)<=0
				sigma_product(i) = (y_i^2)*(y_j^2);
			end	
			
			if magical_number(i)>0
				sigma_product(i) = 0;
			end
	end
	
	if ~isnan(pred_j) == 1 & ~isnan(pred_i) == 1 %%this pair does not break connectedness
				
				%see if I can solve prediction problem for j, after dropping individuals used to build prediction for i, as well as i himself
				index						 		= find(list(:,3)==focus_j);
				[new_pred_j,PP_j,new_patH_j] 		= shortest_path_double_leave_out(list(index,1),list(index,2),[focus_i;focus_j;patH_i],list,A_s);
				PP_j								= PP_j';
				
				%This estimator worked!
				if isnan(new_pred_j) == 0	
					taxonomy(i) 					= 1;
					sigma_product(i) 				= y_i*(y_i-pred_i)*(y_j-new_pred_j)*y_j;
					MOVERS_i{i}						= patH_i;
					Ppath_PROBLEM_i{i}	 			= PP_i_orig;
					MOVERS_j{i}						= new_patH_j;
					Ppath_PROBLEM_j{i}	 			= PP_j;
					TAXONOMY_I{i}					= taxonomy(i).*ones(size(PP_i_orig,1),1);
					TAXONOMY_J{i}					= taxonomy(i).*ones(size(PP_j,1),1);
					FOCUS_I{i}						= i.*ones(size(PP_i_orig,1),1);
					FOCUS_J{i}						= i.*ones(size(PP_j,1),1);
					PREDIZIONE_i(i)	 				= pred_i;
					PREDIZIONE_j(i)	 				= new_pred_j;
					
				
				end
				
			    if isnan(new_pred_j) == 1  	%This estimator did not work!
						
					%try the other direction
					index						    = find(list(:,3)==focus_i);
					[new_pred_i,PP_i,new_patH_i]    = shortest_path_double_leave_out(list(index,1),list(index,2),[focus_i;focus_j;patH_j],list,A_s);
					PP_i							= PP_i';
						
					if isnan(new_pred_i) == 0 %it worked
								taxonomy(i) 	 	= 1;
		   						sigma_product(i) 	= y_j*(y_j-pred_j)*(y_i-new_pred_i)*y_i; 
		   						Ppath_PROBLEM_i{i}	= PP_i;
								Ppath_PROBLEM_j{i}	= PP_j_orig;
								MOVERS_j{i}			= patH_j;
								MOVERS_i{i}			= new_patH_i;
								TAXONOMY_I{i}		= taxonomy(i).*ones(size(PP_i,1),1);
								TAXONOMY_J{i}		= taxonomy(i).*ones(size(PP_j_orig,1),1);
								FOCUS_I{i}			= i.*ones(size(PP_i,1),1);
								FOCUS_J{i}			= i.*ones(size(PP_j_orig,1),1);
								PREDIZIONE_i(i)	 	= new_pred_i;
								PREDIZIONE_j(i)	 	= pred_j;									
					end
							
				    if isnan(new_pred_i) == 1	%double leave two out did not work, try other tricks now
								go_forward = 1;					
				
	%%%%%%%%%%%%%%%%%%%CASE 1.1 (w.r.t. i)		
										
								if i_no_A2 == 0 & j_no_A2 == 0 & j_in_Ai1 == 0 & j_in_Ai2 == 0 
								
										taxonomy(i) 						= 2;
								
										if i_in_Aj1 == 0
											new_pred_j						= yhat_1_j;
											taxonomyALT(i)					= 1;
										end
								
										if i_in_Aj1 == 1
											new_pred_j						= yhat_2_j;
											taxonomyALT(i)					= 2;
										end
									
										sigma_product(i) 					= sigma_2_i(i)*(y_j-new_pred_j)*y_j; 
										
				
										go_forward = 0;
				
								end
				
				
	%%%%%%%%%%%%%%%%%%%CASE 1.1 (w.r.t. j)
								if i_no_A2 == 0 & j_no_A2 == 0 & i_in_Aj1 == 0 & i_in_Aj2 == 0 & go_forward == 1
								
											taxonomy(i) 					= 2;
								
											if j_in_Ai1 == 0
												new_pred_i					= yhat_1_i;
												taxonomyALT(i)				= 3;
											end
								
											if j_in_Ai1 == 1
												new_pred_i					= yhat_2_i;
												taxonomyALT(i)				= 4;
											end
									
											sigma_product(i) 				= sigma_2_j(i)*(y_i-new_pred_i)*y_i; 
												
								end
				
				
	%%%%%%%%%%%%%%%%%%%CASE 1.2(only applies if symmetry holds)		
								if i_no_A2 == 0 & j_in_Ai1 == 1 | i_no_A2 == 0  & j_in_Ai2 == 1 %now we are in the world where the double leave out did not work. Try with the symmetry trick, making sure that j is actually in the union of A_i1 and A_i
					
									index						 			= find(list(:,3)==focus_j);
									[new_pred_j,PP_j,new_patH_j] 			= shortest_path_double_leave_out(list(index,1),list(index,2),[focus_j;patH_i],list,A_s);
									PP_j									= PP_j';
					
					
									if isnan(new_pred_j) == 0 %it worked!
											taxonomy(i) 					= 3;
											sigma_product(i) 				= 0.5*sigma_2_i(i)*(y_j-new_pred_j)*y_j; 
											
											
											if magical_number(i)			<=0
												sigma_product_ALT(i) 		= 0.5*(y_i^2)*(y_j-pred_j)*y_j + 0.5*(y_j^2)*(y_i-pred_i)*y_i; %this is in case we don't want to use symmetry
											end

											if magical_number(i)			>0
												sigma_product_ALT(i) 		= 0; %this is in case we don't want to use symmetry
											end	
											
											Ppath_PROBLEM_i{i}	 			= PP_i_orig;
											Ppath_PROBLEM_j{i}	 			= PP_j;		
											MOVERS_j{i}						= new_patH_j;
											MOVERS_i{i}						= patH_i;	
											TAXONOMY_I{i}					= taxonomy(i).*ones(size(PP_i_orig,1),1);
											TAXONOMY_J{i}					= taxonomy(i).*ones(size(PP_j,1),1); 
											FOCUS_I{i}						= i.*ones(size(PP_i_orig,1),1);
											FOCUS_J{i}						= i.*ones(size(PP_j,1),1);
											
											PREDIZIONE_i(i)	 				= pred_i;
											PREDIZIONE_j(i)	 				= new_pred_j;
											
											Ppath_ALT_PROBLEM_i{i}	 		= PP_i_orig;
											Ppath_ALT_PROBLEM_j{i}	 		= PP_j_orig;		
											MOVERS_ALT_j{i}					= patH_j;
											MOVERS_ALT_i{i}					= patH_i;
											TAXONOMY_ALT_I{i}				= taxonomy(i).*ones(size(PP_i_orig,1),1);
											TAXONOMY_ALT_J{i}				= taxonomy(i).*ones(size(PP_j_orig,1),1); 
											FOCUS_ALT_I{i}					= i.*ones(size(PP_i_orig,1),1);
											FOCUS_ALT_J{i}					= i.*ones(size(PP_j_orig,1),1);	
											
											
									end
						
									if isnan(new_pred_j) == 1 %it did not work!
								
													%try the other direction
													if j_no_A2 == 0 & i_in_Aj1 == 1 |  j_no_A2 == 0  & i_in_Aj2 == 1
												
																				index						 				= find(list(:,3)==focus_i);
																				[new_pred_i,PP_i,new_patH_i] 				= shortest_path_double_leave_out(list(index,1),list(index,2),[focus_i;patH_j],list,A_s);
																				PP_i										= PP_i';
																			
																				if isnan(new_pred_i) == 1 %it did not work!
																			
																								taxonomy(i) 				= 5;
							
																								if magical_number(i)		<=0
																									sigma_product(i) 		= 0.5*(y_i^2)*(y_j-pred_j)*y_j + 0.5*(y_j^2)*(y_i-pred_i)*y_i;
																								end
			
																								if magical_number(i)		>0
																									sigma_product(i) 		= 0;
																								end
																								
																								Ppath_PROBLEM_i{i}			= PP_i_orig;
																								Ppath_PROBLEM_j{i}			= PP_j_orig;
																								MOVERS_j{i}					= patH_j;
																								MOVERS_i{i}					= patH_i;
																								TAXONOMY_I{i}				= taxonomy(i).*ones(size(PP_i_orig,1),1);
																								TAXONOMY_J{i}				= taxonomy(i).*ones(size(PP_j_orig,1),1);
																								FOCUS_I{i}					= i.*ones(size(PP_i_orig,1),1);
																								FOCUS_J{i}					= i.*ones(size(PP_j_orig,1),1);
																								
																								
																				end
																			
																				if isnan(new_pred_i) == 0 %it worked
																			
																								taxonomy(i) 					= 4;
																								sigma_product(i) 				= 0.5*sigma_2_j(i)*(y_i-new_pred_i)*y_i; 
																								
																								if magical_number(i)			<=0
																									sigma_product_ALT(i) 		= 0.5*(y_i^2)*(y_j-pred_j)*y_j + 0.5*(y_j^2)*(y_i-pred_i)*y_i; %this is in case we don't want to use symmetry
																								end

																								if magical_number(i)			>0
																									sigma_product_ALT(i) 		= 0; %this is in case we don't want to use symmetry
																								end	
																								
																								Ppath_PROBLEM_i{i}				= PP_i;
																								Ppath_PROBLEM_j{i}				= PP_j_orig;
																								MOVERS_j{i}						= patH_j;
																								MOVERS_i{i}						= new_patH_i;
																								TAXONOMY_I{i}					= taxonomy(i).*ones(size(PP_i,1),1);
																								TAXONOMY_J{i}					= taxonomy(i).*ones(size(PP_j_orig,1),1);	
																								FOCUS_I{i}						= i.*ones(size(PP_i,1),1);
																								FOCUS_J{i}						= i.*ones(size(PP_j_orig,1),1);
																								PREDIZIONE_i(i)	 				= new_pred_i;
																								PREDIZIONE_j(i)	 				= pred_j;
																								
																								Ppath_ALT_PROBLEM_i{i}	 		= PP_i_orig;
																								Ppath_ALT_PROBLEM_j{i}	 		= PP_j_orig;		
																								MOVERS_ALT_j{i}					= patH_j;
																								MOVERS_ALT_i{i}					= patH_i;
																								TAXONOMY_ALT_I{i}				= taxonomy(i).*ones(size(PP_i_orig,1),1);
																								TAXONOMY_ALT_J{i}				= taxonomy(i).*ones(size(PP_j_orig,1),1); 
																								FOCUS_ALT_I{i}					= i.*ones(size(PP_i_orig,1),1);
																								FOCUS_ALT_J{i}					= i.*ones(size(PP_j_orig,1),1);		
																								
																				end
												
													end %chiudo altra direzione fallita as well.
									end %chiudo simmetria fallita 
								end %chiudo simmetria
						end	 %chiudo secondo tentativo fallito double  
					end	%chiudo primo tentativo fallito double	    	
								
		if taxonomy(i) == -9 % all remaining pairs that do not break connectivity but for which all our remaining tricks did not work
								
								taxonomy(i) 				= 5;
		
								if magical_number(i)		<=0
									sigma_product(i) 		= 0.5*(y_i^2)*(y_j-pred_j)*y_j + 0.5*(y_j^2)*(y_i-pred_i)*y_i;
								end

								if magical_number(i)		>0
									sigma_product(i) 		= 0;
								end
								
								Ppath_PROBLEM_i{i}			= PP_i_orig;
								Ppath_PROBLEM_j{i}			= PP_j_orig;
								MOVERS_j{i}					= patH_j;
								MOVERS_i{i}					= patH_i;
								TAXONOMY_I{i}				= taxonomy(i).*ones(size(PP_i_orig,1),1);
								TAXONOMY_J{i}				= taxonomy(i).*ones(size(PP_j_orig,1),1);
								FOCUS_I{i}					= i.*ones(size(PP_i_orig,1),1);
								FOCUS_J{i}					= i.*ones(size(PP_j_orig,1),1);
								
								
		end %chiudo remaining pairs
   end %chiudo connectivity		
   
if isnan(sigma_product(i))==1
   
 	error('something is wrong')
	   
end
 																		  				
end	%chiudo parfor	
disp('Time to find joint leave 2 out path estimators')
toc


%Stuff for MC
MOVERS_i		= vertcat(MOVERS_i{:,1});
MOVERS_j		= vertcat(MOVERS_j{:,1});
Ppath_PROBLEM_i = vertcat(Ppath_PROBLEM_i{:,1});
Ppath_PROBLEM_j = vertcat(Ppath_PROBLEM_j{:,1});
TAXONOMY_I 		= vertcat(TAXONOMY_I{:,1});
TAXONOMY_J 		= vertcat(TAXONOMY_J{:,1});
FOCUS_I 		= vertcat(FOCUS_I{:,1});
FOCUS_J 		= vertcat(FOCUS_J{:,1});


MOVERS_ALT_i	= vertcat(MOVERS_ALT_i{:,1});
MOVERS_ALT_j	= vertcat(MOVERS_ALT_j{:,1});
Ppath_ALT_PROBLEM_i = vertcat(Ppath_ALT_PROBLEM_i{:,1});
Ppath_ALT_PROBLEM_j = vertcat(Ppath_ALT_PROBLEM_j{:,1});
TAXONOMY_ALT_I 		= vertcat(TAXONOMY_ALT_I{:,1});
TAXONOMY_ALT_J 		= vertcat(TAXONOMY_ALT_J{:,1});
FOCUS_ALT_I 		= vertcat(FOCUS_ALT_I{:,1});
FOCUS_ALT_J 		= vertcat(FOCUS_ALT_J{:,1});

%Summarize now
C_share = zeros(1,6);
for pp =1:6
				C_share(:,pp)=  mean(taxonomy==pp);
end

disp('Relative weight of (i,l) pairs for which it was possible to create a joint leave 2 out estimator')
				C_share_1=C_share(:,1)

disp('Relative weight of (i,l) pairs for which it was possible to apply the yet unexplored idea')
				C_share_2=C_share(:,2)
				
disp('Relative weight of (i,l) pairs for which imposing symmetry worked')
				C_share_3=C_share(:,3)
				
disp('Relative weight of (i,l) pairs for which imposing symmetry worked')
				C_share_4=C_share(:,4)
				
disp('Relative weight of (i,l) pairs for which none of the above worked (but (i,l) do not break connectivity)')
				C_share_5=C_share(:,5)
				
disp('Relative weight of (i,l) pairs that break connectivity)')
				C_share_6=C_share(:,6)

%Now building the four pieces that compose the sampling variance estimator of the firm effects.
first_piece 		= 4*sum((vec.^2).*sigma_2);

%Use Hutchison trick and this (https://en.wikipedia.org/wiki/Hadamard_product_(matrices)) to get the second piece
	I_Lambda_P		= (speye(Ndelta,Ndelta)-Lambda_P);
	L_P				= ichol(I_Lambda_P,struct('type','ict','droptol',1e-2,'diagcomp',2));
    VCM 			= parallel.pool.Constant(diag(sigma_1));
    Fdelta  		= parallel.pool.Constant(Fdelta);
    L  				= parallel.pool.Constant(L);
    Lambda_B		= parallel.pool.Constant(Lambda_B);
    I_Lambda_P		= parallel.pool.Constant(I_Lambda_P);
    L_P  			= parallel.pool.Constant(L_P);
    F 				= parallel.pool.Constant(F);   
    %NSIM			= floor(20*0.05^(-2)*log(2/0.05)); %w.p 95% the absolute difference b/w the trace estimator and the actual trace is going to be less than 5% - see table 1 of http://www.cs.tau.ac.il/~stoledo/Bib/Pubs/trace3.pdf 
    NSIM			= 5000;
    aux_SIM			= zeros(NSIM,1);
    parfor s=1:NSIM    
                   xsimul			= randn(Ndelta,1);       
				   left				= C_build_FD(xsimul,Fdelta.Value,L.Value,pfun_,Lambda_B.Value,I_Lambda_P.Value,L_P.Value,F.Value);
				   xsimul			= VCM.Value*xsimul;
				   right			= C_build_FD(xsimul,Fdelta.Value,L.Value,pfun_,Lambda_B.Value,I_Lambda_P.Value,L_P.Value,F.Value);
				   aux_SIM(s)		= left'*VCM.Value*right;
    end
	second_piece		= 2*mean(aux_SIM);
	Dsym				= 2*D + 2*D';
	third_piece			= 2*sigma_1'*Dsym*sigma_1;
	fourth_piece		= 4*sum((magical_number).*(sigma_product-(sigma_1_i.*sigma_1_j))); %multiply by 4 instead of 2 because I only considered upper triangular part of the matrix P.

%Finished (this is all assuming symmetry thus far)
	V_sym				= first_piece - second_piece - third_piece - fourth_piece;

%If you are not willing to assume symmetry, then you need to replace the fourth piece (associated with sigma product)
	sel					= (taxonomy==3 | taxonomy==4); %this is where we used the symmetry trick
	sigma_product(sel)	= sigma_product_ALT(sel);
	fourth_piece		= 4*sum((magical_number).*(sigma_product-(sigma_1_i.*sigma_1_j))); %multiply by 4 instead of 2 because I only considered upper triangular part of the matrix P.
	V_nosym				= first_piece - second_piece - third_piece - fourth_piece;



    