function [Lambda_P, Lambda_B_fe] = eff_res_inputA(A_right,A_left,deMean,Q,X,xx,Lchol,N,J,K,elist,clustering_level,movers,T,type_algorithm,id,firmid,match_id,epsilon)
%This function calculates, using parallel coding, (Pii,Bii) for a general 
%two way fixed effects model. The code is likely to be slow on large
%datasets.




%Dimensions
NT=size(X,1);
M=size(elist,1);

%PreCreate
Pii=zeros(M,1);
Bii_fe=zeros(M,1);
Bii_cov=zeros(M,1);
Bii_pe=zeros(M,1);

%Options for solver
numIterations = 300; %iteration for the pcg solver
tol=1e-6; %tol for pcg
tolProb=0.5;

%Objects for parfor 
    xx_c = parallel.pool.Constant(xx);
    X_c =  parallel.pool.Constant(X);
    Xright=parallel.pool.Constant(X(elist(:,2),:));
    Xleft=parallel.pool.Constant(X(elist(:,1),:));

%Special case of Laplacian design matrix    
    Nmatches=max(match_id);
    match_id_movers=match_id(movers);
    firmid_movers=firmid(movers);
    id_movers=id(movers);
    [~,sel,~]=unique(match_id_movers,'rows');
    match_id_movers=match_id_movers(sel);
    firmid_movers=firmid_movers(sel);
    id_movers=id_movers(sel);
    maxT=max(T(~movers));
    count=ones(NT,1);
    gcs = cell2mat(accumarray(id,count,[],@(x){cumsum(x)}));
    sel_stayers=(gcs==1).*(~movers);
    sel_stayers=sel_stayers>0;
    stayers_matches_sel=match_id(sel_stayers);
    Tinv=1./T;
    elist_JLL=[id_movers N+firmid_movers id_movers N+firmid_movers];
    M=size(elist_JLL,1);
    Pii_movers=zeros(M,1);
    Bii_fe_movers=zeros(M,1);
    Bii_cov_movers=zeros(M,1);
    Bii_pe_movers=zeros(M,1);

    
    if  strcmp(type_algorithm,'exact')
        Xright = sparse((1:M)',elist_JLL(:,1),1,M,N+J);
        Xright =  Xright+sparse((1:M)',elist_JLL(:,2),-1,M,N+J);
        Xright=parallel.pool.Constant(Xright);
        
        parfor i=1:M
        [xtilde, flag]= pcg(xx_c.Value,Xright.Value(i,:)',tol,numIterations,Lchol);

        %Statistical Leverage
        Pii_movers(i)=Xright.Value(i,:)*xtilde;

        
        %Bii
        aux_right=A_right*xtilde;
        aux_left=A_left*xtilde;
        
        if deMean == 1
        	aux_right 			= aux_right - mean(aux_right);
        	aux_left  			=	aux_left -  mean(aux_left);
        	Bii_fe_movers(i)	=aux_left'*aux_right;
        end
        
        if deMean == 0
        	Bii_fe_movers(i)	=aux_left'*Q*aux_right;
        end
        
        end
    	
    end	
    
   if  strcmp(type_algorithm,'JLL') 
       
       %Number of random draws to implement Random Projection Algorithm.
       disp('# of Simulated Projections for JLL:')
       
       %scale = ceil(log(N+J)/epsilon^2); Actual bound provided by S-S
       %(2011). 
       
       %In my experience, the formula above is too conservative. One can obtain very
       %good approximation in very reasonable computation time using the
       %following alternative found here: http://www.cs.cmu.edu/~jkoutis/EffectiveResistances/EffectiveResistances.m
       
       scale = ceil(log2(NT))/epsilon
       
       %Auxiliary components
       X_right=parallel.pool.Constant(A_right);
       X_left=parallel.pool.Constant(A_left);  
       X_pe=parallel.pool.Constant([X(:,1:N) sparse(NT,J)]);
       elist_1=parallel.pool.Constant(elist_JLL(:,1));
       elist_2=parallel.pool.Constant(elist_JLL(:,2));
       elist_3=parallel.pool.Constant(elist_JLL(:,3));
       elist_4=parallel.pool.Constant(elist_JLL(:,4));        
       
       parfor i=1:scale
                
                %To avoid redundant warnings from each worker.
                Z=0;
                ZB_right=0;
                ZB_left	=0;
                
                %Random Projection Matrix (Rademacher)
                ons = (rand(1,NT) > tolProb);
                ons = ons - not(ons);
                ons = ons./sqrt(scale);
                
                %Lambda_P
                [Z, flag]= pcg(xx_c.Value,(ons*(X_c.Value))',tol,numIterations,Lchol);
                
                %Lambda_B
                ons = (rand(1,size(A_right,1)) > tolProb);
                ons = ons - not(ons);
                ons = ons./sqrt(scale);
                if deMean == 1
                ons=ons-mean(ons);
                end
                
                [ZB_right, flag]	  = pcg(xx_c.Value,(ons*X_right.Value)',tol,numIterations,Lchol);
                [ZB_left,  flag]      = pcg(xx_c.Value,(ons*X_left.Value)',tol,numIterations,Lchol);
                
                %Collect results of the given draw
                Pii_movers=Pii_movers+(((Z(elist_1.Value)-Z(elist_2.Value))).*(((Z(elist_3.Value)-Z(elist_4.Value)))));
                Bii_fe_movers=Bii_fe_movers+(((ZB_left(elist_1.Value)-ZB_left(elist_2.Value))).*(((ZB_right(elist_3.Value)-ZB_right(elist_4.Value)))));
       end
                
   end
     
%Assign step
Pii_movers=sparse(match_id_movers,1,Pii_movers,Nmatches,1);
Pii_stayers=sparse(stayers_matches_sel,1,Tinv(sel_stayers),Nmatches,1);
Pii=Pii_movers+Pii_stayers;
Bii_fe=sparse(match_id_movers,1,Bii_fe_movers,Nmatches,1);
  
%Create the matrices.
rows=elist(:,1);
column=elist(:,2);
index_cluster=elist(:,3);

Pii=Pii(index_cluster);
Bii_fe=Bii_fe(index_cluster);

%Lambda P
Lambda_P=sparse(rows,column,Pii,NT,NT);
Lambda_P=Lambda_P+triu(Lambda_P,1)'; %make it symmetric.

%Lambda B var(fe)
Lambda_B_fe=sparse(rows,column,Bii_fe,NT,NT);
Lambda_B_fe=Lambda_B_fe+triu(Lambda_B_fe,1)'; %make it symmetric.

end

