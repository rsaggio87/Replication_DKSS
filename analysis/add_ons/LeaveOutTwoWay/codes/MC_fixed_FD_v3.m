function  [V share_bp COV_R1] = MC_fixed_FD_v3(ydelta,A_i2_in_Aj1,A_i1_in_Aj2,A_i2_in_Aj2,y,Fdelta,L,F,tol,epsilon,type_algorithm,pfun_,Lambda_P,Lambda_B,Cij,D,vec,id_movers,movers_sel,index_movers,index_movers_i,index_movers_2,index_movers_i_2,Ppath_1,Ppath_2,j_in_Ai1,j_in_Ai2,i_in_Aj1,i_in_Aj2,no_A2,magical_number,prob_pairs,taxonomy,vec2,Dq,magical_number_q,sigma_i,x1bar,Lambda_B2,lambda_1)

%GRANDE MEDIA
GRANDE_MEDIA 	= mean(y);

%Q-share
Q_share 		= mean(no_A2==1);

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


%Now get the sigmas and additional information updated for this new vector of y.
sigma_1			= (ymovers-yhat1).*ymovers;
sigma_2			= (ymovers-yhat1).*(ymovers-yhat2);
sel				= no_A2(movers_sel)==1;
sigma_2(sel)	= ymovers(sel).^2;



%Update \sigma_i for dudes for which we cannot find a second estimator
sigma_i(sel)	=  (ymovers(sel)-GRANDE_MEDIA).^2;

%Reshape
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

no_A2_ii		= no_A2(prob_pairs(:,1));
no_A2_jj		= no_A2(prob_pairs(:,2));

%Finish up
[V, share_bp, COV_R1]	= finish_SE(ydelta,taxonomy,A_i2_in_Aj1,A_i1_in_Aj2,A_i2_in_Aj2,vec2,pfun_,prob_pairs,j_in_Ai1,j_in_Ai2,i_in_Aj1,i_in_Aj2,index_movers_i,index_movers,index_movers_i_2,index_movers_2,no_A2_ii,no_A2_jj,Q_share,y_ii,y_jj,sigma_1_i,sigma_1_j,sigma_2_i,sigma_2_j,yhat_1_ii,yhat_2_ii,yhat_1_jj,yhat_2_jj,vec,sigma_2,Ndelta,Lambda_P,sigma_1,Fdelta,L,Lambda_B,do_weak,Lambda_B2,F,x1bar,lambda_1,D,magical_number,magical_number_q,Dq,sigma_i,GRANDE_MEDIA);
end