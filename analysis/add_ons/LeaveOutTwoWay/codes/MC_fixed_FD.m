function  [V_sym V_nosym]	 = MC_fixed_FD(y,Chamard,Dsym,vec,id_movers,movers_sel,index_movers,index_movers_i,index_movers_2,index_movers_i_2,Ppath_1,Ppath_2,no_A2,list_problem,taxonomy,magical_number,TAXONOMY_I,TAXONOMY_J,MOVERS_i,MOVERS_j,FOCUS_I, FOCUS_J,Ppath_PROBLEM_i,Ppath_PROBLEM_j,taxonomyALT,TAXONOMY_ALT_I,TAXONOMY_ALT_J,MOVERS_ALT_i,MOVERS_ALT_j,FOCUS_ALT_I, FOCUS_ALT_J,Ppath_ALT_PROBLEM_i,Ppath_ALT_PROBLEM_j)

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
sigma_1			= sparse(id_movers,1,sigma_1,Ndelta,1); %reshape in full Ndelta space (counting stayers as well);
sigma_2			= sparse(id_movers,1,sigma_2,Ndelta,1); %reshape in full Ndelta space (counting stayers as well);
yhat1			= sparse(id_movers,1,yhat1,Ndelta,1);   %reshape in full Ndelta space (counting stayers as well);
yhat2			= sparse(id_movers,1,yhat2,Ndelta,1);   %reshape in full Ndelta space (counting stayers as well);

%Bring important information in the correct space now (the \mathcal{P} set)
y_i				= y(list_problem(:,1));
y_j				= y(list_problem(:,2));

sigma_1_i		= sigma_1(list_problem(:,1));
sigma_1_j		= sigma_1(list_problem(:,2));

sigma_2_i		= sigma_2(list_problem(:,1));
sigma_2_j		= sigma_2(list_problem(:,2));

yhat_1_ii		= yhat1(list_problem(:,1));
yhat_1_jj		= yhat1(list_problem(:,2));
yhat_2_ii		= yhat2(list_problem(:,1));
yhat_2_jj		= yhat2(list_problem(:,2));
Npaths			= size(list_problem,1);

%Updating the sigma_products
sigma_product	= zeros(Npaths,1);

%Now Taxonomy equal to 6
sel				= find(taxonomy == 6);
to_use			= (y_i.^2).*(y_j.^2).*(magical_number<=0)+0.*(magical_number>0);
sigma_product	= sigma_product + sparse(sel,1,to_use(sel),Npaths,1);


%Now Taxonomy equal to 5
sel					= find(taxonomy == 5);
if size(sel,1)>0
	[pred_i pred_j] = prediction_DOUBLE_LEAVE(5,taxonomy,TAXONOMY_I,TAXONOMY_J,MOVERS_i,MOVERS_j,FOCUS_I, FOCUS_J,Ppath_PROBLEM_i,Ppath_PROBLEM_j,y);
	to_use			= (0.5*(y_i.^2).*(y_j-pred_j).*y_j + 0.5*(y_j.^2).*(y_i-pred_i).*y_i).*(magical_number<=0) + 0.*(magical_number>0);
	sigma_product	= sigma_product + sparse(sel,1,to_use(sel),Npaths,1);
end	

%Now Taxonomy equal to 4
sel					= find(taxonomy == 4);
if size(sel,1)>0
	[pred_i pred_j] = prediction_DOUBLE_LEAVE(4,taxonomy,TAXONOMY_I,TAXONOMY_J,MOVERS_i,MOVERS_j,FOCUS_I, FOCUS_J,Ppath_PROBLEM_i,Ppath_PROBLEM_j,y);
	to_use			= .5*sigma_2_j.*(y_i-pred_i).*y_i; 
	sigma_product	= sigma_product + sparse(sel,1,to_use(sel),Npaths,1);
end

%Now Taxonomy equal to 3
sel					= find(taxonomy == 3);
if size(sel,1)>0
	[pred_i pred_j] = prediction_DOUBLE_LEAVE(3,taxonomy,TAXONOMY_I,TAXONOMY_J,MOVERS_i,MOVERS_j,FOCUS_I, FOCUS_J,Ppath_PROBLEM_i,Ppath_PROBLEM_j,y);	
	to_use			= 0.5*sigma_2_i.*(y_j-pred_j).*y_j;
	sigma_product	= sigma_product + sparse(sel,1,to_use(sel),Npaths,1); %can't upgrade to zero because sigma_product reads the conservative estimator that does not impose symmetry
end

%Now Taxonomy equal to 2
to_use			= sigma_2_i.*(y_j-yhat_1_jj).*y_j;
sel				= find(taxonomy == 2 & taxonomyALT == 1 );
sigma_product	= sigma_product + sparse(sel,1,to_use(sel),Npaths,1);

to_use			= sigma_2_i.*(y_j-yhat_2_jj).*y_j;
sel				= find(taxonomy == 2 & taxonomyALT == 2 );
sigma_product	= sigma_product + sparse(sel,1,to_use(sel),Npaths,1);

to_use			= sigma_2_j.*(y_i-yhat_1_ii).*y_i;
sel				= find(taxonomy == 2 & taxonomyALT == 3 );
sigma_product	= sigma_product + sparse(sel,1,to_use(sel),Npaths,1);

to_use			= sigma_2_j.*(y_i-yhat_2_ii).*y_i;
sel				= find(taxonomy == 2 & taxonomyALT == 4 );
sigma_product	= sigma_product + sparse(sel,1,to_use(sel),Npaths,1);

%Now Taxonomy equal to 1
sel					= find(taxonomy == 1);
if size(sel,1)>0
	[pred_i pred_j] = prediction_DOUBLE_LEAVE(1,taxonomy,TAXONOMY_I,TAXONOMY_J,MOVERS_i,MOVERS_j,FOCUS_I, FOCUS_J,Ppath_PROBLEM_i,Ppath_PROBLEM_j,y);
	to_use			= y_i.*(y_i-pred_i).*(y_j-pred_j).*y_j;
	sigma_product	= sigma_product + sparse(sel,1,to_use(sel),Npaths,1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Finish part that imposes symmetry

%Close
first_piece 		= 4*sum((vec.^2).*sigma_2);
second_piece		= 2*sigma_1'*Chamard*sigma_1;
third_piece			= 2*sigma_1'*Dsym*sigma_1;
fourth_piece		= 4*sum((magical_number).*(sigma_product-(sigma_1_i.*sigma_1_j))); %multiply by 4 instead of 2 because I only considered upper triangular part of the matrix P.

%Finished (this is all assuming symmetry thus far)
V_sym				= first_piece - second_piece - third_piece - fourth_piece;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Now reset the predictions that impose symmetry
sel					= find(taxonomy == 4 | taxonomy==3);
sigma_product(sel) 	= 0;


sel					= find(taxonomy == 4);
if size(sel,1)>0
	[pred_i pred_j] = prediction_DOUBLE_LEAVE(4,taxonomy,TAXONOMY_ALT_I,TAXONOMY_ALT_J,MOVERS_ALT_i,MOVERS_ALT_j,FOCUS_ALT_I, FOCUS_ALT_J,Ppath_ALT_PROBLEM_i,Ppath_ALT_PROBLEM_j,y);
	to_use			= (0.5*(y_i.^2).*(y_j-pred_j).*y_j + 0.5*(y_j.^2).*(y_i-pred_i).*y_i).*(magical_number<=0) + 0.*(magical_number>0);
	sigma_product	= sigma_product + sparse(sel,1,to_use(sel),Npaths,1);
end

sel					= find(taxonomy == 3);

if size(sel,1)>0
	[pred_i pred_j] = prediction_DOUBLE_LEAVE(3,taxonomy,TAXONOMY_ALT_I,TAXONOMY_ALT_J,MOVERS_ALT_i,MOVERS_ALT_j,FOCUS_ALT_I, FOCUS_ALT_J,Ppath_ALT_PROBLEM_i,Ppath_ALT_PROBLEM_j,y);	
	to_use			= (0.5*(y_i.^2).*(y_j-pred_j).*y_j + 0.5*(y_j.^2).*(y_i-pred_i).*y_i).*(magical_number<=0) + 0.*(magical_number>0);
	sigma_product	= sigma_product + sparse(sel,1,to_use(sel),Npaths,1); 
	
end


fourth_piece		= 4*sum((magical_number).*(sigma_product-(sigma_1_i.*sigma_1_j))); %multiply by 4 instead of 2 because I only considered upper triangular part of the matrix P.
V_nosym				= first_piece - second_piece - third_piece - fourth_piece;

    