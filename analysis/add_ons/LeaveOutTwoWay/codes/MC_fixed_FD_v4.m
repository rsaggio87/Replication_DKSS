function  [V_sym V_nosym COV_R1] = MC_fixed_FD_v4(C,Cq,y,Fdelta,L,F,tol,epsilon,type_algorithm,pfun_,Lambda_P,Lambda_B,Cij,D,vec,id_movers,movers_sel,index_movers,index_movers_i,index_movers_2,index_movers_i_2,Ppath_1,Ppath_2,j_in_Ai1,j_in_Ai2,i_in_Aj1,i_in_Aj2,no_A2,magical_number,prob_pairs,index_circle1,index_circle2,focus_circle1_i,focus_circle2_i,focus_circle1_j,focus_circle2_j,PGIGA_CIRCLE1,PGIGA_CIRCLE2,vec2,Dq,magical_number_q,sigma_i,x1bar,Lambda_B2,lambda_1)

do_weak			= 1;
Ndelta			= size(y,1);
ymovers			= y(movers_sel);
	
%Update the first leave one out estimator
y_others		= y(index_movers).*Ppath_1;
[aa,bb,cc]		= unique(index_movers_i);
yhat1 			= splitapply(@sum,y_others,cc);

%Update the second leave one out estimator
y_others		= y(index_movers_2).*Ppath_2;
[aa,bb,cc]		= unique(index_movers_i_2);
yhat2 			= splitapply(@sum,y_others,cc);

%Update the first circle estimator
y_others		= y(index_circle1).*PGIGA_CIRCLE1;
[aa,bb,cc]		= unique([focus_circle1_i focus_circle1_j],'rows','stable');
ycircle1 		= splitapply(@sum,y_others,cc);

%Update the second circle estimator
y_others		= y(index_circle2).*PGIGA_CIRCLE2;
[aa,bb,cc]		= unique([focus_circle2_i focus_circle2_j],'rows','stable');
ycircle2 		= splitapply(@sum,y_others,cc);

%Now get the sigmas and additional information updated for this new vector of y.
sigma_1			= (ymovers-yhat1).*ymovers;
sigma_2			= (ymovers-yhat1).*(ymovers-yhat2);
sel				= no_A2(movers_sel)==1;
sigma_2(sel)	= ymovers(sel).^2;
sigma_1			= sparse(id_movers,1,sigma_1,Ndelta,1); %reshape in full Ndelta space (counting stayers as well);
sigma_2			= sparse(id_movers,1,sigma_2,Ndelta,1); %reshape in full Ndelta space (counting stayers as well);
yhat1			= sparse(id_movers,1,yhat1,Ndelta,1);   %reshape in full Ndelta space (counting stayers as well);
yhat2			= sparse(id_movers,1,yhat2,Ndelta,1);   %reshape in full Ndelta space (counting stayers as well);

%Bring important information in the correct space now (the \mathcal{P} set)
y_ii			= y(prob_pairs(:,1));
y_jj			= y(prob_pairs(:,2));

sigma_1_i		= sigma_1(prob_pairs(:,1));
sigma_1_j		= sigma_1(prob_pairs(:,2));

sigma_2_i		= sigma_2(prob_pairs(:,1));
sigma_2_j		= sigma_2(prob_pairs(:,2));

yhat_1_ii		= yhat1(prob_pairs(:,1));
yhat_1_jj		= yhat1(prob_pairs(:,2));
yhat_2_ii		= yhat2(prob_pairs(:,1));
yhat_2_jj		= yhat2(prob_pairs(:,2));


%This is very important: for those problematic pairs because j is in either A_i1 or A_i2, update the ycircle
NPROB			= size(prob_pairs,1);
ycircle			= sparse(NPROB,1);
index			= find(j_in_Ai1);
bad_pairs		= prob_pairs(j_in_Ai1==1,:);
sel				= ismember([index_movers_i index_movers],bad_pairs,'rows');
ycircle			= ycircle + sparse(index,1, ycircle1(sel),NPROB,1);
index			= find(j_in_Ai2);
bad_pairs		= prob_pairs(j_in_Ai2==1,:);
sel				= ismember([index_movers_i_2 index_movers_2],bad_pairs,'rows');
ycircle			= ycircle + sparse(index,1, ycircle2(sel),NPROB,1);

%Ready to fill the product of sigmas for the problematic pairs now
sigma_product			= sparse(NPROB,1);

%Now update the good pairs
good_pairs 				= j_in_Ai1 + j_in_Ai2;
sel	 					= find(good_pairs==0 & i_in_Aj1==0);
numero					= sigma_2_i.*(y_jj).*(y_jj-yhat_1_jj);
sigma_product			= sigma_product + sparse(sel,1,numero(sel),NPROB,1);
sel	 					= find(good_pairs==0 & i_in_Aj1==1);
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
	first_piece 		= 4*sum((vec.^2).*sigma_2);
	second_piece		= 2*sigma_1'*(C.^2)*sigma_1;

%The last two pieces
	Dsym				= 2*D + 2*D';
	third_piece			= 2*sigma_1'*Dsym*sigma_1;
	fourth_piece		= 4*sum((magical_number).*(sigma_product-(sigma_1_i.*sigma_1_j))); %multiply by 4 instead of 2 because I only considered upper triangular part of the that identifies problematic pairs, see the command [r,c,~]			= find(triu(prob_pairs));.

%Finished
	V_sym				= first_piece - second_piece - third_piece - fourth_piece;

%Estimator that does not rely on symmetry
	index				= find(j_in_Ai2);
	numero				= (magical_number(index)<0).*(sigma_1_i(index)).*(y_jj(index).^2)+0.*(magical_number(index)>=0);
	sigma_product(index)= numero;

	index				= find(j_in_Ai1);
	numero				= (magical_number(index)<0).*(y_ii(index).*(y_ii(index)-yhat_2_ii(index))).*(y_jj(index).^2)+0.*(magical_number(index)>=0);
	sigma_product(index)= numero;

	fourth_piece		= 4*sum((magical_number).*(sigma_product-(sigma_1_i.*sigma_1_j))); %multiply by 4 instead of 2 because I only considered upper triangular part of the that identifies problematic pairs

	V_nosym				= first_piece - second_piece - third_piece - fourth_piece;
	
%Weak ID Confidence Intervals
	first_piece 		= 4*sum((vec2.^2).*sigma_2);
	second_piece		= 2*sigma_1'*(Cq.^2)*sigma_1;
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