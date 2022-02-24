function  [V, share_bp, COV_R1, A_i2_in_Aj1, A_i1_in_Aj2, A_i2_in_Aj2 Ppath_1, Ppath_2, no_A2, prob_pairs, j_in_Ai1, j_in_Ai2, i_in_Aj1, i_in_Aj2, Cij, D, magical_number, index_movers,index_movers_i,index_movers_2,index_movers_i_2,taxonomy,Dq,magical_number_q] = paths_network_v3(id_movers,firmid_delta,firmid_delta_f,ydelta,vec,Fdelta,L,F,tol,epsilon,type_,pfun_,Lambda_P,Lambda_B,x1bar,lambda_1,vec2,sigma_i,SUM_EIG)



%Read whether the function has to calculate weak id VCM
if nargin > 14
		do_weak 		= 1;
end

if nargin <= 14
		do_weak 		= 0;
end	

if do_weak == 0
		COV_R1 			= 0;
		Lambda_B2		= 0;
		lambda_1		= 0;
		x1bar			= 0;
end

if do_weak == 1
			Lambda_B2=Lambda_B;
			Q_use=size(x1bar,2);
			for qqq=1:Q_use
				norm_fact		= lambda_1(qqq)*(x1bar(:,qqq).^2);
				norm_fact		= spdiags(norm_fact,0,size(ydelta,1),size(ydelta,1));
				Lambda_B2		= Lambda_B2-norm_fact;
			end
end


%Construct the list of movers, tell me who is a unique mover
Ndelta			 = size(ydelta,1);
index_C			 = (1:Ndelta)';
list			 = [firmid_delta_f  firmid_delta id_movers ydelta];
list_aux 		 = [firmid_delta firmid_delta_f id_movers ydelta];
sel   			 = firmid_delta~=firmid_delta_f; %focus on movers only
id_movers		 = id_movers(sel);
%GRANDE MEDIA
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
tic
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


toc

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

%Now get cross products of the C matrix
tic
P				= diag(Lambda_P);
B				= diag(Lambda_B);
mega_list		= [focus_m index_i; focus_m index_j];
NMEGA			= size(mega_list,1);
[mega_list,~,index_mega]= unique(mega_list,'rows','stable'); %individual cross products that I need to query.
elist			= [firmid_delta(mega_list(:,1)) firmid_delta_f(mega_list(:,1)) firmid_delta(mega_list(:,2)) firmid_delta_f(mega_list(:,2))];
[elist,~,index_unique]	= unique(elist,'rows');
[Mij, Bij] 		= eff_res_FAST_FE_ONLY(elist,Fdelta,L,F,tol,epsilon,type_,pfun_); %Apply function
Mij				= -Mij(index_unique);%Reshape into the original elist space where M=I-P;
Bij				= Bij(index_unique);%Reshape into the original elist space
Mii				= 1-P(mega_list(:,1)); %notice T=2 requirement here.
Mjj				= 1-P(mega_list(:,2));
Bii				= B(mega_list(:,1));
Bjj				= B(mega_list(:,2));
Cij				= Bij - 0.5*Mij.*(Mii.\Bii) -0.5*Mij.*(Mjj.\Bjj); %top of page 20
C_ij_fix		= Cij;
Cij				= Cij(index_mega);
Cij				= reshape(Cij,[NMEGA/2,2]);
numero			= Cij(:,1).*Cij(:,2).*P_i1.*P_j2;
D				= sparse(index_i,index_j,numero,Ndelta,Ndelta);
disp('Time to build the matrix D')
toc

%Update the C matrix for Weak id
if do_weak == 1
	Cijq=C_ij_fix;
	for qqq=1:Q_use
	wij			= x1bar(mega_list(:,1),qqq).*x1bar(mega_list(:,2),qqq);
	wii			= x1bar(mega_list(:,1),qqq).*x1bar(mega_list(:,1),qqq);
	wjj			= x1bar(mega_list(:,2),qqq).*x1bar(mega_list(:,2),qqq);
	Cijq		= Cijq -wij*lambda_1(qqq) + 0.5*lambda_1(qqq)*Mij.*(Mii.\wii) + 0.5*lambda_1(qqq)*Mij.*(Mjj.\wjj) ;
	end
	Cijq		= Cijq(index_mega);
	Cijq		= reshape(Cijq,[NMEGA/2,2]);
	numero		= Cijq(:,1).*Cijq(:,2).*P_i1.*P_j2;
	Dq			= sparse(index_i,index_j,numero,Ndelta,Ndelta);
	clear C_ij_fix
end

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
elist					= [firmid_delta(r) firmid_delta_f(r) firmid_delta(c) firmid_delta_f(c)];
[elist,~,index_unique]	= unique(elist,'rows');
[Mij, Bij] 				= eff_res_FAST_FE_ONLY(elist,Fdelta,L,F,tol,epsilon,type_,pfun_); %Apply function
Mij						= -Mij(index_unique);%Reshape into the original elist space where M=I-P;
Bij						= Bij(index_unique);%Reshape into the original elist space
Mii						= 1-P(r); %notice T=2 requirement here.
Mjj						= 1-P(c);
Bii						= B(r);
Bjj						= B(c);
Cij						= Bij - 0.5*Mij.*(Mii.\Bii) -0.5*Mij.*(Mjj.\Bjj); %top of page 20

if do_weak == 1
		Cijq=Cij;
		for qqq=1:Q_use
		wij			= x1bar(r,qqq).*x1bar(c,qqq);
		wii			= x1bar(r,qqq).*x1bar(r,qqq);
		wjj			= x1bar(c,qqq).*x1bar(c,qqq);
		Cijq		= Cijq -wij*lambda_1(qqq) + 0.5*lambda_1(qqq)*Mij.*(Mii.\wii) + 0.5*lambda_1(qqq)*Mij.*(Mjj.\wjj) ;
		end
end

Cij						= sparse(r,c,Cij,max(list(:,3)),max(list(:,3))); %make it a sparse matrix.
Cij						= Cij'+triu(Cij,1); %make it symmetric
Cij						= Cij(sub2ind(size(Cij),prob_pairs(:,1),prob_pairs(:,2))); %back to be a vector but now properly indexed
magical_number			= Cij.^2 +  2*(D(sub2ind(size(D),prob_pairs(:,1),prob_pairs(:,2)))) + 2*(D(sub2ind(size(D),prob_pairs(:,2),prob_pairs(:,1))))  ;

if do_weak == 1
		Cijq			= sparse(r,c,Cijq,max(list(:,3)),max(list(:,3))); %make it a sparse matrix.
		Cijq			= Cijq'+triu(Cijq,1); %make it symmetric
		Cijq			= Cijq(sub2ind(size(Cijq),prob_pairs(:,1),prob_pairs(:,2))); %back to be a vector but now properly indexed
		magical_number_q= Cijq.^2 +  2*(Dq(sub2ind(size(Dq),prob_pairs(:,1),prob_pairs(:,2)))) + 2*(Dq(sub2ind(size(Dq),prob_pairs(:,2),prob_pairs(:,1))));
end


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

%Finish up
[V, share_bp, COV_R1]	= finish_SE(ydelta,taxonomy,A_i2_in_Aj1,A_i1_in_Aj2,A_i2_in_Aj2,vec2,pfun_,prob_pairs,j_in_Ai1,j_in_Ai2,i_in_Aj1,i_in_Aj2,index_movers_i,index_movers,index_movers_i_2,index_movers_2,no_A2_ii,no_A2_jj,Q_share,y_ii,y_jj,sigma_1_i,sigma_1_j,sigma_2_i,sigma_2_j,yhat_1_ii,yhat_2_ii,yhat_1_jj,yhat_2_jj,vec,sigma_2,Ndelta,Lambda_P,sigma_1,Fdelta,L,Lambda_B,do_weak,Lambda_B2,F,x1bar,lambda_1,D,magical_number,magical_number_q,Dq,sigma_i,GRANDE_MEDIA);


end




    