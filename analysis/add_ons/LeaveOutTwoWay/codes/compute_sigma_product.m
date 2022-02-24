function [sigma_product bad_pairs tag] = compute_sigma_product(taxonomy,no_A2_jj,no_A2_ii,j_in_Ai1,j_in_Ai2,i_in_Aj1,i_in_Aj2,A_i2_in_Aj1,A_i1_in_Aj2,A_i2_in_Aj2,sigma_1_i,sigma_2_i,sigma_3_i,sigma_4_i,sigma_1_j,sigma_2_j,sigma_3_j,sigma_4_j,magical_number);

NPROB					= size(taxonomy,1);

sigma_product			= sparse(NPROB,1);
tag					    = sparse(NPROB,1);

%Line 2
sel 					= find(taxonomy == 1 & tag == 0 & j_in_Ai1==1 & i_in_Aj1 == 0 & A_i2_in_Aj1 == 0);
numero					= sigma_3_i.*sigma_1_j;
sigma_product			= sigma_product + sparse(sel,1,numero(sel),NPROB,1);
tag						= tag + sparse(sel,1,2,NPROB,1);

%Line 3
sel 					= find(taxonomy == 1 & tag == 0 & i_in_Aj1 == 1 & j_in_Ai1==0 & A_i1_in_Aj2 == 0);
numero					= sigma_3_j.*sigma_1_i;	
sigma_product			= sigma_product + sparse(sel,1,numero(sel),NPROB,1);
tag						= tag + sparse(sel,1,3,NPROB,1);

%Line 4
sel 					= find(taxonomy == 1 & tag == 0 & i_in_Aj1 == 1 & j_in_Ai1==1 & A_i2_in_Aj2 == 0);
numero					= sigma_3_j.*sigma_3_i;
sigma_product			= sigma_product + sparse(sel,1,numero(sel),NPROB,1);
tag						= tag + sparse(sel,1,4,NPROB,1);

%Line 5
sel 					= find(taxonomy == 1 & tag == 0 & j_in_Ai1== 0 & j_in_Ai2 == 0 & i_in_Aj1==0);
numero					= sigma_2_i.*sigma_1_j;
sigma_product			= sigma_product + sparse(sel,1,numero(sel),NPROB,1);
tag						= tag + sparse(sel,1,5,NPROB,1);

%Line 6
sel 					= find(taxonomy == 1 & tag == 0 & j_in_Ai1== 0 & j_in_Ai2 == 0 & i_in_Aj1==1);
numero					= sigma_2_i.*sigma_3_j;
sigma_product			= sigma_product + sparse(sel,1,numero(sel),NPROB,1);
tag						= tag + sparse(sel,1,6,NPROB,1);

%Line 7
sel 					= find(taxonomy == 1 & tag == 0 & i_in_Aj1== 0 & i_in_Aj2 == 0 & j_in_Ai1==0);
numero					= sigma_1_i.*sigma_2_j;
sigma_product			= sigma_product + sparse(sel,1,numero(sel),NPROB,1);
tag						= tag + sparse(sel,1,7,NPROB,1);

%Line 8
sel 					= find(taxonomy == 1 & tag == 0 & i_in_Aj1== 0 & i_in_Aj2 == 0 & j_in_Ai1==1);
numero					= sigma_3_i.*sigma_2_j;
sigma_product			= sigma_product + sparse(sel,1,numero(sel),NPROB,1);
tag						= tag + sparse(sel,1,8,NPROB,1);

%Line 9 
sel 					= find(taxonomy == 1 & tag == 0 & j_in_Ai1== 0);
numero					= (sigma_1_i.*sigma_4_j).*(magical_number<0);
sigma_product			= sigma_product + sparse(sel,1,numero(sel),NPROB,1);
tag						= tag + sparse(sel,1,9,NPROB,1);

%Line 10 
sel 					= find(taxonomy == 1 & tag == 0 & j_in_Ai1== 1);
numero					= (sigma_3_i.*sigma_4_j).*(magical_number<0);
sigma_product			= sigma_product + sparse(sel,1,numero(sel),NPROB,1);
tag						= tag + sparse(sel,1,10,NPROB,1);
sel						= find(tag>=9);
bad_pairs				= sparse(sel,1,1,NPROB,1);

 
%%%%%%% Taxonomy==2, i.e. cases where either Q_i or Q_j >0
tag					    = sparse(NPROB,1);
%Line 2
sel 					= find(taxonomy == 2 & tag == 0 & no_A2_ii== 0 & j_in_Ai1==1 & i_in_Aj2 == 0 & A_i2_in_Aj1 == 0);
numero					= sigma_3_i.*sigma_1_j;
sigma_product			= sigma_product + sparse(sel,1,numero(sel),NPROB,1);
tag						= tag + sparse(sel,1,1,NPROB,1);
	
%Line 3
sel 					= find(taxonomy == 2 & tag == 0 & no_A2_jj == 0 &  i_in_Aj1 == 1 & j_in_Ai1==0 & A_i1_in_Aj2 == 0);
numero					= sigma_3_j.*sigma_1_i;
sigma_product			= sigma_product + sparse(sel,1,numero(sel),NPROB,1);
tag						= tag + sparse(sel,1,1,NPROB,1);

%Line 4
sel 					= find(taxonomy == 2 & tag == 0 & no_A2_ii== 0 & j_in_Ai1==0 & j_in_Ai2==0 & i_in_Aj1==0);
numero					= sigma_2_i.*sigma_1_j;
sigma_product			= sigma_product + sparse(sel,1,numero(sel),NPROB,1);
tag						= tag + sparse(sel,1,1,NPROB,1);

%Line 5
sel 					= find(taxonomy == 2 & tag == 0 & no_A2_jj== 0 & i_in_Aj1==0 & i_in_Aj2==0 & j_in_Ai1==0);
numero					= sigma_2_j.*sigma_1_i;
sigma_product			= sigma_product + sparse(sel,1,numero(sel),NPROB,1);
tag						= tag + sparse(sel,1,1,NPROB,1);

%Line 6
sel 					= find(taxonomy == 2 & tag == 0 & j_in_Ai1==0);
numero					= (sigma_1_i.*sigma_4_j).*(magical_number<0)+0.*(magical_number>=0);
sigma_product			= sigma_product + sparse(sel,1,numero(sel),NPROB,1);
tag						= tag + sparse(sel,1,-9,NPROB,1);

%Line 7
sel 					= find(taxonomy == 2 & tag == 0 & no_A2_ii== 0  & j_in_Ai1==1);
numero					= (sigma_3_i.*sigma_4_j).*(magical_number<0)+0.*(magical_number>=0);
sigma_product			= sigma_product + sparse(sel,1,numero(sel),NPROB,1);
tag						= tag + sparse(sel,1,-9,NPROB,1);

%Line 8
sel 					= find(taxonomy == 2 & tag == 0 & i_in_Aj1==0);
numero					= (sigma_4_i.*sigma_1_j).*(magical_number<0)+0.*(magical_number>=0);
sigma_product			= sigma_product + sparse(sel,1,numero(sel),NPROB,1);
tag						= tag + sparse(sel,1,-9,NPROB,1);

%Line 9
sel 					= find(taxonomy == 2 & tag == 0 & no_A2_jj == 0 &  i_in_Aj1==1);
numero					= (sigma_4_i.*sigma_3_j).*(magical_number<0)+0.*(magical_number>=0);
sigma_product			= sigma_product + sparse(sel,1,numero(sel),NPROB,1);
tag						= tag + sparse(sel,1,-9,NPROB,1);

%Line 10
sel 					= find(taxonomy == 2 & tag == 0);
numero					= (sigma_4_i.*sigma_4_j).*(magical_number<0)+0.*(magical_number>=0);
sigma_product			= sigma_product + sparse(sel,1,numero(sel),NPROB,1);
tag						= tag + sparse(sel,1,-9,NPROB,1);
tabulate(tag)

%Calculate bad pairs
sel						= find(tag>=6);
bad_pairs				= bad_pairs+sparse(sel,1,1,NPROB,1);
tabulate(bad_pairs)

end    


