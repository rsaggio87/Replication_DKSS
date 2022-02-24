function  [V_sym, V_nosym, share_bp, COV_R1, Ppath_1, Ppath_2, no_A2, prob_pairs, j_in_Ai1, j_in_Ai2, i_in_Aj1, i_in_Aj2, Cij, D, magical_number, index_movers,index_movers_i,index_movers_2,index_movers_i_2,index_circle1,index_circle2,focus_circle1_i,focus_circle2_i,focus_circle1_j,focus_circle2_j,PGIGA_CIRCLE1,PGIGA_CIRCLE2,Dq,magical_number_q] = paths_network_v4(C,Cq,id_movers,firmid_delta,firmid_delta_f,ydelta,vec,Fdelta,L,F,tol,epsilon,type_,pfun_,Lambda_P,Lambda_B,x1bar,lambda_1,vec2,sigma_i,SUM_EIG)

%Read whether the function has to calculate weak id VCM
if nargin > 14
		do_weak 		= 1;
end

if nargin <= 14
		do_weak 		= 0;
end	

if do_weak == 0
		COV_R1 			= 0;
end

if do_weak == 1
		
		Lambda_B2		= lambda_1*(x1bar.^2);
		Lambda_B2		= spdiags(Lambda_B2,0,size(ydelta,1),size(ydelta,1));
		Lambda_B2		= Lambda_B-Lambda_B2;
		
		norm_fact		= (lambda_1*x1bar)*sum(x1bar.*ydelta);
	
end


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
ycircle1		 = cell(Nmovers,1);
ycircle2		 = cell(Nmovers,1);
Pcircle1		 = cell(Nmovers,1);
Pcircle2		 = cell(Nmovers,1);
index_circle1	 = cell(Nmovers,1);
index_circle2	 = cell(Nmovers,1);
focus_circle1_i	 = cell(Nmovers,1);
focus_circle1_j	 = cell(Nmovers,1);
focus_circle2_i	 = cell(Nmovers,1);
focus_circle2_j	 = cell(Nmovers,1);
PGIGA_CIRCLE1	 = cell(Nmovers,1);
PGIGA_CIRCLE2	 = cell(Nmovers,1);
no_A2			 = zeros(Nmovers,1);


toc

%Now we enter the important part of the code. Here we will try to form two predictors of ydelta_i after leaving i out.
%We will use a shortest path algorithm to find the first prediction, which always exists because the input data is from a leave out connected set.
%The second prediction might exist or not exist. We'll record whether that is the case.

%Along the way we will also be filling in the elements of the auxiliary matrix D which will be useful for the next part

tic
parfor i = 1:Nmovers
		
		%Find shortest leave one out path, this will be our first estimator. 
		[yhat1(i),P1,list_path_aux,list_edges_aux,auxiliary{i},auxiliary_edges{i},~,~,ycircle1{i},Pcircle1{i},index_circle1{i},focus_circle1_i{i},focus_circle1_j{i},PGIGA_CIRCLE1{i}] = shortest_path_double_leave_out(list(i,1),list(i,2),list(i,3),list,A_s);
		
		list_path{i} 	= list_path_aux;
		list_edges{i}	= list_edges_aux;
		Ppath_1{i}		= P1';
		
		%Second estimator of yhat, after leaving indexes associated with first estimator.  
		[yhat2(i),P2,list_path_aux2,~,auxiliary_2{i},~,~,~,ycircle2{i},Pcircle2{i},index_circle2{i},focus_circle2_i{i},focus_circle2_j{i},PGIGA_CIRCLE2{i}] = shortest_path_double_leave_out(list(i,1),list(i,2),[list(i,3);list_path_aux],list,A_s);
		
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
			ycircle2{i}		= 0;
			Pcircle2{i}		= 0;	
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

%Summarize what we have learned thus far
sel   = no_A2==1;
disp('share of workers for which is not possible to form two non-overlapping estimators of Delta y_i')
Q_share = mean(sel)

%Build \tilde{\sigma}
sigma_1			= (list(:,4)-yhat1).*list(:,4);
sigma_2 		= (list(:,4)-yhat1).*(list(:,4)-yhat2);
sigma_2(sel)	= list(sel,4).^2;

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
ycircle1 		= vertcat(ycircle1{:,1}); 
ycircle2 		= vertcat(ycircle2{:,1});
Pcircle1 		= vertcat(Pcircle1{:,1}); 
Pcircle2 		= vertcat(Pcircle2{:,1});
index_circle1 	= vertcat(index_circle1{:,1}); 
index_circle2 	= vertcat(index_circle2{:,1});
focus_circle1_i = vertcat(focus_circle1_i{:,1}); 
focus_circle2_i = vertcat(focus_circle2_i{:,1});
focus_circle1_j = vertcat(focus_circle1_j{:,1}); 
focus_circle2_j = vertcat(focus_circle2_j{:,1});
PGIGA_CIRCLE1 	= vertcat(PGIGA_CIRCLE1{:,1}); 
PGIGA_CIRCLE2 	= vertcat(PGIGA_CIRCLE2{:,1});


%Now get cross products of the C matrix
tic
Cim				= C(sub2ind(size(C),focus_m,index_i));
Cjm				= C(sub2ind(size(C),focus_m,index_j));
numero			= Cim.*Cjm.*P_i1.*P_j2;
D				= sparse(index_i,index_j,numero,Ndelta,Ndelta);
disp('Time to build the matrix D')
toc

%Update the C matrix for Weak id
if do_weak == 1
Cim				= Cq(sub2ind(size(Cq),focus_m,index_i));
Cjm				= Cq(sub2ind(size(Cq),focus_m,index_j));
numero			= Cim.*Cjm.*P_i1.*P_j2;
Dq				= sparse(index_i,index_j,numero,Ndelta,Ndelta);
end

%Reshape things into a new world where for each i we provide the list of indexes in A_i1. Similarly for A_i2
index_movers 	= vertcat(list_path{:,1});    %this is A_i1.
index_movers_i 	= vertcat(auxiliary{:,1}); 

index_movers_2	= vertcat(list_path_2{:,1});  %this is A_i2.
index_movers_i_2= vertcat(auxiliary_2{:,1}); 



%Now find the problematic pairs
paths_problem_i 			= cell(Nmovers,1);
paths_problem_j 			= cell(Nmovers,1);

tic
parfor i = 1:Nmovers	
	   sel			 		= ismember(index_movers_i,list(i,3));
	   A_i1	 		 		= [index_movers(sel)]; %focus worker and all the movers used to derive the first estimator of y for the focus worker.
	   
	   others_path	 		= [index_movers(~sel)];
	   others_index	 		= [index_movers_i(~sel)];
	   
	   problem1 			= ismember(others_path,A_i1);   %a given path in A_ell,1 overlaps with A_i,1
	   problem2 			= ismember(others_index,A_i1);  %ell is in A_i,
	   problem3 			= ismember(others_path,list(i,3));  %i is in A_ell,1
		
	   sel 					= (problem1+problem2+problem3>0);
	   
	   Nsel				 	= sum(sel);
	   list_prob			= [list(i,3).*ones(Nsel,1) others_index(sel)];
	   sel					= list_prob(:,1)>list_prob(:,2);
	   list_prob(sel,:)		= [list_prob(sel,2) list_prob(sel,1)]; %notice the normalization here, essentially only looking at upper triangular part of the (hypothetical) binary nxn matrix that reports 1 whenever (i,l) are problematic
	   	
	   paths_problem_i{i}	= list_prob(:,1);
	   paths_problem_j{i}	= list_prob(:,2);	    	
	
end
prob_pairs 					= [vertcat(paths_problem_i{:,1}) vertcat(paths_problem_j{:,1})];
prob_pairs    				= unique(prob_pairs,'rows','stable'); %suppose A_i1={j,3,5,6} and A_j1 = {i,3,5,6}. Then (i,j) are clearly problematic and they will have "multiple" problems, as many paths in A_j1 overlap with A_i1. I don't need to keep track of all the (3,3), (5,5), etc, so just tell me that (i,j) are indeed problematic. Also double counting because of the normalization impose above (if i is in Aj1 and j is in Ai1 then the corresponding will have two rows with i,l where i<l)
disp('Time to find the problematic pairs')
toc

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


%Pin down C_ij of the problematic pairs
Cij						= C(sub2ind(size(C),prob_pairs(:,1),prob_pairs(:,2)));
Cijq					= Cq(sub2ind(size(Cq),prob_pairs(:,1),prob_pairs(:,2)));

%Now figure out the taxonomy of the problematic pair: Start by looking whether \ell is in A_i1 or A_i2
j_in_Ai1				= ismember(prob_pairs,[index_movers_i index_movers],'rows');
j_in_Ai2				= ismember(prob_pairs,[index_movers_i_2 index_movers_2],'rows');
tabulate(j_in_Ai1+j_in_Ai2)

%Now whether i is in A_j1 or Aj2
i_in_Aj1				= ismember([prob_pairs(:,2) prob_pairs(:,1)],[index_movers_i index_movers],'rows');
i_in_Aj2				= ismember([prob_pairs(:,2) prob_pairs(:,1)],[index_movers_i_2 index_movers_2],'rows');
tabulate(i_in_Aj1+i_in_Aj2)

%This is very important: for those problematic pairs because j is in either A_i1 or A_i2, update the ycircle
NPROB					= size(prob_pairs,1);
ycircle					= sparse(NPROB,1);
index					= find(j_in_Ai1);
bad_pairs				= prob_pairs(j_in_Ai1==1,:);
sel						= ismember([index_movers_i index_movers],bad_pairs,'rows');
ycircle					= ycircle + sparse(index,1, ycircle1(sel),NPROB,1);
index					= find(j_in_Ai2);
bad_pairs				= prob_pairs(j_in_Ai2==1,:);
sel						= ismember([index_movers_i_2 index_movers_2],bad_pairs,'rows');
ycircle					= ycircle + sparse(index,1, ycircle2(sel),NPROB,1);

%Ready to fill the product of sigmas for the problematic pairs now
sigma_product			= sparse(NPROB,1);

%Now update the good pairs
bad_pairs 				= j_in_Ai1 + j_in_Ai2;
sel	 					= find(bad_pairs==0 & i_in_Aj1==0);
numero					= sigma_2_i.*(y_jj).*(y_jj-yhat_1_jj);
sigma_product			= sigma_product + sparse(sel,1,numero(sel),NPROB,1);
sel	 					= find(bad_pairs==0 & i_in_Aj1==1);
numero					= sigma_2_i.*(y_jj).*(y_jj-yhat_2_jj);
sigma_product			= sigma_product + sparse(sel,1,numero(sel),NPROB,1);

%Now update the bad pairs
index					= find(j_in_Ai1);
numero					= 0.5*(sigma_2_i(j_in_Ai1).*(y_jj(j_in_Ai1)).*(y_jj(j_in_Ai1)-ycircle(j_in_Ai1)));
sigma_product			= sigma_product + sparse(index,1,numero,NPROB,1);

index					= find(j_in_Ai2);
numero					= 0.5*(sigma_2_i(j_in_Ai2).*(y_jj(j_in_Ai2)).*(y_jj(j_in_Ai2)-ycircle(j_in_Ai2)));
sigma_product			= sigma_product + sparse(index,1,numero,NPROB,1);
sigma_product_fixed		= sigma_product;


%Now building the four pieces that compose the sampling variance estimator of the firm effects.
first_piece 			= 4*sum((vec.^2).*sigma_2);

%Second piece
second_piece			= 2*sigma_1'*(C.^2)*sigma_1;

%The last two pieces
magical_number			= Cij.^2 +  2*(D(sub2ind(size(D),prob_pairs(:,1),prob_pairs(:,2)))) + 2*(D(sub2ind(size(D),prob_pairs(:,2),prob_pairs(:,1))))  ;
Dsym					= 2*D + 2*D';
third_piece				= 2*sigma_1'*Dsym*sigma_1;
fourth_piece			= 4*sum((magical_number).*(sigma_product-(sigma_1_i.*sigma_1_j))); %multiply by 4 instead of 2 because I only considered upper triangular part of the that identifies problematic pairs

%Finished
V_sym					= first_piece - second_piece - third_piece - fourth_piece;

%Estimator that does not rely on symmetry
index					= find(j_in_Ai2);
numero					= (magical_number(index)<0).*(sigma_1_i(index)).*(y_jj(index).^2)+0.*(magical_number(index)>=0);
sigma_product(index)	= numero;

index					= find(j_in_Ai1);
numero					= (magical_number(index)<0).*(y_ii(index).*(y_ii(index)-yhat_2_ii(index))).*(y_jj(index).^2)+0.*(magical_number(index)>=0);
sigma_product(index)	= numero;

fourth_piece			= 4*sum((magical_number).*(sigma_product-(sigma_1_i.*sigma_1_j))); %multiply by 4 instead of 2 because I only considered upper triangular part of the that identifies problematic pairs

V_nosym					= first_piece - second_piece - third_piece - fourth_piece;


%VCM for weak id
if do_weak == 1
	first_piece 		= 4*sum((vec2.^2).*sigma_2);
	second_piece		= 2*sigma_1'*(Cq.^2)*sigma_1;
	Dsym				= 2*Dq + 2*Dq';
	third_piece			= 2*sigma_1'*Dsym*sigma_1;
	magical_number_q	= Cijq.^2 +  2*(Dq(sub2ind(size(Dq),prob_pairs(:,1),prob_pairs(:,2)))) + 2*(Dq(sub2ind(size(Dq),prob_pairs(:,2),prob_pairs(:,1))))  ;
	fourth_piece		= 4*sum((magical_number_q).*(sigma_product_fixed-(sigma_1_i.*sigma_1_j)));
	COV_R1(2,2)			= first_piece - second_piece - third_piece - fourth_piece;
	if COV_R1(2,2)<0
				COV_R1(2,2) = first_piece;
	end
	COV_R1(1,1)			= sum(x1bar.^2.*sigma_i);
	COV_R1(1,2)			= 2*sum(x1bar.*sigma_2.*vec2);
	COV_R1(2,1)			= 2*sum(x1bar.*sigma_2.*vec2);
end


%Share of bad pairs
share_bp				= sum(sigma_product_fixed(bad_pairs==1))/V_sym





    