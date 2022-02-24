function  [V_sym, V_nosym, share_bp, COV_R1, Ppath_1, Ppath_2, no_A2, prob_pairs, j_in_Ai1, j_in_Ai2, i_in_Aj1, i_in_Aj2, Cij, D, magical_number, index_movers,index_movers_i,index_movers_2,index_movers_i_2,index_circle1,index_circle2,focus_circle1_i,focus_circle2_i,focus_circle1_j,focus_circle2_j,PGIGA_CIRCLE1,PGIGA_CIRCLE2,taxonomy,Dq,magical_number_q] = paths_network_v3(id_movers,firmid_delta,firmid_delta_f,ydelta,vec,Fdelta,L,F,tol,epsilon,type_,pfun_,Lambda_P,Lambda_B,x1bar,lambda_1,vec2,sigma_i,SUM_EIG)

%GRANDE MEDIA
	GRANDE_MEDIA = mean(ydelta);


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
		Lambda_B2		= lambda_1*(x1bar.^2);
		Lambda_B2		= spdiags(Lambda_B2,0,size(ydelta,1),size(ydelta,1));
		Lambda_B2		= Lambda_B-Lambda_B2;	
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
			yhat2(i)		= NaN;
			list_path_2{i}	= list(i,3);
			auxiliary_2{i}	= list(i,3);
			focus_m{i}		= [];
			index_i{i}		= [];
			index_j{i}		= [];
			P_i1{i}			= [];
			P_j2{i}			= [];
			Ppath_2{i}		= 1;
			ycircle2{i}		= NaN;
			index_circle2{i}= list(i,3);
			focus_circle2_i{i}=list(i,3);
			focus_circle2_j{i}=list(i,3);
			PGIGA_CIRCLE2{i}= list(i,3);	
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
Cij				= Cij(index_mega);
Cij				= reshape(Cij,[NMEGA/2,2]);
numero			= Cij(:,1).*Cij(:,2).*P_i1.*P_j2;
D				= sparse(index_i,index_j,numero,Ndelta,Ndelta);
disp('Time to build the matrix D')
toc

%Update the C matrix for Weak id
if do_weak == 1
	wij			= x1bar(mega_list(:,1)).*x1bar(mega_list(:,2));
	wii			= x1bar(mega_list(:,1)).*x1bar(mega_list(:,1));
	wjj			= x1bar(mega_list(:,2)).*x1bar(mega_list(:,2));
	Cijq		= Bij - 0.5*Mij.*(Mii.\Bii) -0.5*Mij.*(Mjj.\Bjj) -wij*lambda_1 + 0.5*lambda_1*Mij.*(Mii.\wii) + 0.5*lambda_1*Mij.*(Mjj.\wjj) ;
	Cijq		= Cijq(index_mega);
	Cijq		= reshape(Cijq,[NMEGA/2,2]);
	numero		= Cijq(:,1).*Cijq(:,2).*P_i1.*P_j2;
	Dq			= sparse(index_i,index_j,numero,Ndelta,Ndelta);
end

%Reshape things into a new world where for each i we provide the list of indexes in A_i1. Similarly for A_i2
index_movers 	= vertcat(list_path{:,1});    %this is A_i1.
index_movers_i 	= vertcat(auxiliary{:,1}); 

index_movers_2	= vertcat(list_path_2{:,1});  %this is A_i2.
index_movers_i_2= vertcat(auxiliary_2{:,1}); 



%Now find the problematic pairs. This is where we will end up splitting the code
paths_problem_i = cell(Nmovers,1);
paths_problem_j = cell(Nmovers,1);

tic
parfor i = 1:Nmovers	
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
	   
	   if Q_share == 0
	   sel				= list_prob(:,1)>list_prob(:,2);
	   list_prob(sel,:)	= [list_prob(sel,2) list_prob(sel,1)]; %notice the normalization here, essentially only looking at upper triangular part of the (hypothetical) binary nxn matrix that reports 1 whenever (i,l) are problematic
	   end
	   	
	   paths_problem_i{i}= list_prob(:,1);
	   paths_problem_j{i}= list_prob(:,2);	    	
	
end
prob_pairs 				= [vertcat(paths_problem_i{:,1}) vertcat(paths_problem_j{:,1})];
prob_pairs    			= unique(prob_pairs,'rows','stable'); %suppose A_i1={j,3,5,6} and A_j1 = {i,3,5,6}. Then (i,j) are clearly problematic and they will have "multiple" problems, as many paths in A_j1 overlap with A_i1. I don't need to keep track of all the (3,3), (5,5), etc, so just tell me that (i,j) are indeed problematic. Also double counting because of the normalization impose above (if i is in Aj1 and j is in Ai1 then the corresponding will have two rows with i,l where i<l)
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

no_A2_ii				= no_A2(prob_pairs(:,1));
no_A2_jj				= no_A2(prob_pairs(:,2));



%Pin down C_ij of the problematic pairs
elist					= [firmid_delta(prob_pairs(:,1)) firmid_delta_f(prob_pairs(:,1)) firmid_delta(prob_pairs(:,2)) firmid_delta_f(prob_pairs(:,2))];
[elist,~,index_unique]	= unique(elist,'rows');
[Mij, Bij] 				= eff_res_FAST_FE_ONLY(elist,Fdelta,L,F,tol,epsilon,type_,pfun_); %Apply function
Mij						= -Mij(index_unique);%Reshape into the original elist space where M=I-P;
Bij						= Bij(index_unique);%Reshape into the original elist space
Mii						= 1-P(prob_pairs(:,1)); %notice T=2 requirement here.
Mjj						= 1-P(prob_pairs(:,2));
Bii						= B(prob_pairs(:,1));
Bjj						= B(prob_pairs(:,2));
Cij						= Bij - 0.5*Mij.*(Mii.\Bii) -0.5*Mij.*(Mjj.\Bjj); %top of page 20
magical_number			= Cij.^2 +  2*(D(sub2ind(size(D),prob_pairs(:,1),prob_pairs(:,2)))) + 2*(D(sub2ind(size(D),prob_pairs(:,2),prob_pairs(:,1))))  ;
mean(magical_number)

if do_weak == 1
		wij				= x1bar(prob_pairs(:,1)).*x1bar(prob_pairs(:,2));
		wii				= x1bar(prob_pairs(:,1)).*x1bar(prob_pairs(:,1));
		wjj				= x1bar(prob_pairs(:,2)).*x1bar(prob_pairs(:,2));
		Cijq			= Bij - 0.5*Mij.*(Mii.\Bii) -0.5*Mij.*(Mjj.\Bjj) -wij*lambda_1 + 0.5*lambda_1*Mij.*(Mii.\wii) + 0.5*lambda_1*Mij.*(Mjj.\wjj) ;
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


%Classify the problematic pairs now
taxonomy 				= 1.*(no_A2_ii==0) + 2.*(j_in_Ai1==0 & no_A2_ii== 1) + 3.*(j_in_Ai1== 1 & no_A2_ii== 1);
tabulate(taxonomy)

%Ready to fill the product of sigmas for the problematic pairs now
sigma_product			= sparse(NPROB,1);

if Q_share == 0 
	%Now update the good-bad pairs
	bad_pairs 			= j_in_Ai1 + j_in_Ai2;
	sel	 				= find(bad_pairs==0 & i_in_Aj1==0);
	numero				= sigma_2_i.*(y_jj).*(y_jj-yhat_1_jj);
	sigma_product		= sigma_product + sparse(sel,1,numero(sel),NPROB,1);
	sel	 				= find(bad_pairs==0 & i_in_Aj1==1);
	numero				= sigma_2_i.*(y_jj).*(y_jj-yhat_2_jj);
	sigma_product		= sigma_product + sparse(sel,1,numero(sel),NPROB,1);

	%Now update the bad-bad pairs
	index				= find(j_in_Ai1);
	numero				= 0.5*(sigma_2_i(j_in_Ai1).*(y_jj(j_in_Ai1)).*(y_jj(j_in_Ai1)-ycircle(j_in_Ai1))); %1/2 in front because I am taking a single individual per node, so code still needs 
	sigma_product		= sigma_product + sparse(index,1,numero,NPROB,1);

	index				= find(j_in_Ai2);
	numero				= 0.5*(sigma_2_i(j_in_Ai2).*(y_jj(j_in_Ai2)).*(y_jj(j_in_Ai2)-ycircle(j_in_Ai2))); %1/2 in front because I am taking a single individual per node, it would be different otherwise.
	sigma_product		= sigma_product + sparse(index,1,numero,NPROB,1);
	sigma_product_fixed	= sigma_product;
	
	disp('full mean')
	mean(sigma_product)
	disp('bad pairs =1')
	mean(sigma_product(bad_pairs==1))
	disp('bad pairs =0')
	mean(sigma_product(bad_pairs==0))
	
end

if Q_share > 0 
	bad_pairs 			= (j_in_Ai1== 1 & taxonomy==1 | j_in_Ai2 == 1 & taxonomy==1);
	
	%Case #1 (bad-good pairs)
	sel	 				= find(j_in_Ai1== 0 & j_in_Ai2 == 0 & taxonomy==1 & i_in_Aj1==0);
	numero				= sigma_2_i.*(y_jj).*(y_jj-yhat_1_jj);
	sigma_product		= sigma_product + sparse(sel,1,numero(sel),NPROB,1);
	
	sel	 				= find(j_in_Ai1== 0 & j_in_Ai2 == 0 & taxonomy==1 & i_in_Aj1==1);
	numero				= sigma_2_i.*(y_jj).*(y_jj-yhat_2_jj);
	sigma_product		= sigma_product + sparse(sel,1,numero(sel),NPROB,1);

	%Case #1 (bad-bad pairs)
	sel					= (j_in_Ai1==1 & taxonomy==1);
	index				= find(sel);
	numero				= 0.5*(sigma_2_i(sel).*(y_jj(sel)).*(y_jj(sel)-ycircle(sel))); %1/2 in front because I am taking a single individual per node, so code still needs 
	sigma_product		= sigma_product + sparse(index,1,numero,NPROB,1);

	sel					= (j_in_Ai2==1 & taxonomy==1);
	index				= find(sel);
	numero				= 0.5*(sigma_2_i(sel).*(y_jj(sel)).*(y_jj(sel)-ycircle(sel))); %1/2 in front because I am taking a single individual per node, it would be different otherwise.
	sigma_product		= sigma_product + sparse(index,1,numero,NPROB,1);
	
	%Case #2
	sel					= taxonomy == 2;
	index				= find(sel);
	numero				= (sigma_1_i.*((y_jj-GRANDE_MEDIA).^2)).*(magical_number<0)+ 0.*(magical_number>=0);
	sigma_product		= sigma_product + sparse(index,1,numero(sel),NPROB,1);
	
	%Case #3
	sel					= taxonomy == 3;
	index				= find(sel);
	numero				= (((y_ii-GRANDE_MEDIA).^2).*((y_jj-GRANDE_MEDIA).^2)).*(magical_number<0)+ 0.*(magical_number>=0);
	sigma_product		= sigma_product + sparse(index,1,numero(sel),NPROB,1);
	
	
	%Is it (i,l) that is truly problematic or (l,i)?
	sel 				= prob_pairs(:,1)<prob_pairs(:,2);
	prob_pairs(sel,:)	= [prob_pairs(sel,2) prob_pairs(sel,1)];
	[~,~,index_prob]	= unique(prob_pairs,'rows','stable');
	findMAx 			= @(x)min(x(x(:,1)==max(x(:,1)),2));
	index				= (1:NPROB)';
	sel			 		= splitapply(findMAx,[taxonomy index],index_prob);
	
	
	%Bring everything in the correct space
	prob_pairs			= prob_pairs(sel,:);
	
	sigma_product		= sigma_product(sel);
	magical_number		= magical_number(sel);
	
	taxonomy			= taxonomy(sel);
	
	j_in_Ai1			= j_in_Ai1(sel);
	j_in_Ai2			= j_in_Ai2(sel);

	y_ii				= y_ii(sel);
	y_jj				= y_jj(sel);

	sigma_1_i			= sigma_1_i(sel);
	sigma_1_j			= sigma_1_j(sel);

	sigma_2_i			= sigma_2_i(sel);
	sigma_2_j			= sigma_2_j(sel);

	yhat_1_ii			= yhat_1_ii(sel);
	yhat_1_jj			= yhat_1_jj(sel);
	yhat_2_ii			= yhat_2_ii(sel);
	yhat_2_jj			= yhat_2_jj(sel);

	magical_number_q	= magical_number_q(sel);
	bad_pairs			= bad_pairs(sel);
	sigma_product_fixed	= sigma_product;
end


%Now building the four pieces that compose the sampling variance estimator of the firm effects.
first_piece 			= 4*sum((vec.^2).*sigma_2)

%Use Hutchison trick and this (https://en.wikipedia.org/wiki/Hadamard_product_(matrices)) to get the second piece
I_Lambda_P				= (speye(Ndelta,Ndelta)-Lambda_P);
L_P						= ichol(I_Lambda_P,struct('type','ict','droptol',1e-2,'diagcomp',2));
VCM 					= parallel.pool.Constant(diag(sigma_1));
Fdelta  				= parallel.pool.Constant(Fdelta);
L  						= parallel.pool.Constant(L);
Lambda_B				= parallel.pool.Constant(Lambda_B);
if do_weak == 1
 	Lambda_B2			= parallel.pool.Constant(Lambda_B2);
end
I_Lambda_P				= parallel.pool.Constant(I_Lambda_P);
L_P  					= parallel.pool.Constant(L_P);
F 						= parallel.pool.Constant(F);   
%NSIM					= floor(20*0.05^(-2)*log(2/0.05)); %w.p 95% the absolute difference b/w the trace estimator and the actual trace is going to be less than 5% - see table 1 of http://www.cs.tau.ac.il/~stoledo/Bib/Pubs/trace3.pdf 
NSIM					= 1000;
aux_SIM					= zeros(NSIM,1);
aux_SIM_2				= zeros(NSIM,1);
rng(1234)
XSIMUL					= randn(Ndelta,NSIM);

    parfor s=1:NSIM    
        %xsimul			= randn(Ndelta,1);
        xsimul			= XSIMUL(:,s);       
	 	left			= C_build_FD(xsimul,Fdelta.Value,L.Value,pfun_,Lambda_B.Value,I_Lambda_P.Value,L_P.Value,F.Value);
	 	
	 	if do_weak == 1
	 		left_2		= C_build_FD(xsimul,Fdelta.Value,L.Value,pfun_,Lambda_B2.Value,I_Lambda_P.Value,L_P.Value,F.Value);
	 		norm_fact	= (lambda_1*x1bar)*sum(x1bar.*xsimul);
	 		left_2 		= left_2 - norm_fact;
	 	end
	 	
		xsimul			= VCM.Value*xsimul;
		right			= C_build_FD(xsimul,Fdelta.Value,L.Value,pfun_,Lambda_B.Value,I_Lambda_P.Value,L_P.Value,F.Value);
		
		if do_weak == 1
	 		right_2		= C_build_FD(xsimul,Fdelta.Value,L.Value,pfun_,Lambda_B2.Value,I_Lambda_P.Value,L_P.Value,F.Value);
	 		norm_fact	= (lambda_1*x1bar)*sum(x1bar.*xsimul);
	 		right_2 	= right_2 - norm_fact;
	 	end
		
		aux_SIM(s)		= left'*VCM.Value*right;
		
		if do_weak == 1
			aux_SIM_2(s)= left_2'*VCM.Value*right_2;
		end
    end  

second_piece			= 2*mean(aux_SIM)

%The last two pieces
Dsym					= 2*D + 2*D';
third_piece				= 2*sigma_1'*Dsym*sigma_1
fourth_piece			= 4*sum((magical_number).*(sigma_product-(sigma_1_i.*sigma_1_j))) %multiply by 4 instead of 2 because I only considered upper triangular part of the that identifies problematic pairs
%lista_old =[prob_pairs magical_number sigma_product sigma_1_i.*sigma_1_j (magical_number).*(sigma_product-(sigma_1_i.*sigma_1_j))];
%save('old')
%Finished
V_sym					= first_piece - second_piece - third_piece - fourth_piece;

%Estimator that does not rely on symmetry
index					= find(j_in_Ai2==1 & taxonomy==1);
numero					= (magical_number(index)<0).*(sigma_1_i(index).*((y_jj(index)-GRANDE_MEDIA).^2))+0.*(magical_number(index)>=0);
sigma_product(index)	= numero;

index					= find(j_in_Ai1==1	 & taxonomy==1);
numero					= (magical_number(index)<0).*(y_ii(index).*(y_ii(index)-yhat_2_ii(index))).*((y_jj(index)-GRANDE_MEDIA).^2)+0.*(magical_number(index)>=0);
sigma_product(index)	= numero;

fourth_piece			= 4*sum((magical_number).*(sigma_product-(sigma_1_i.*sigma_1_j))); %multiply by 4 instead of 2 because I only considered upper triangular part of the that identifies problematic pairs

V_nosym					= first_piece - second_piece - third_piece - fourth_piece;


%VCM for weak id
if do_weak == 1
	first_piece 		= 4*sum((vec2.^2).*sigma_2);
	second_piece		= 2*mean(aux_SIM_2);
	Dsym				= 2*Dq + 2*Dq';
	third_piece			= 2*sigma_1'*Dsym*sigma_1;
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
%save('old')
end




    