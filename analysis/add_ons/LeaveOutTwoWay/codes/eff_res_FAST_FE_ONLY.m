function [Pii, Bii] = eff_res_FAST_FE_ONLY(elist,Fdelta,L,F,tol,epsilon,type_,pfun_)
%Computation of Pii, Bii for variance of firm effects in two-way fixed
%effects model in a model that partials out worker effects by taking first
%differences.
    
%% Auxiliary Step
    m= size(Fdelta,1); %dimension of the overall system.
    J = size(Fdelta,2);%number of nodes in the graph
    numIterations = 1000; %iteration for the pcg solver
    tolProb = 0.5;%use to create the johnson lindestraus projection matrix
    M= size(elist,1);
	NT=size(F,1);
    
    
%Preallocate 
tic   
    output 		= sparse(M,1);
    output_B	= sparse(M,1);  
    %L_c 		= parallel.pool.Constant(L);
    %F_c 		= parallel.pool.Constant(F);
    %Fdelta_c 	= parallel.pool.Constant(Fdelta);
    %LCHOL 		= parallel.pool.Constant(@() pfun_);
    %elist_c		= parallel.pool.Constant(elist); 
toc   
    
    
%% Finding leverages and cross products using S-S (2008)
if strcmp(type_,'JLL')   
    disp('Number of Simulated Projections of the Matrix')
    scale = ceil(log2(m))/epsilon
    optimset('display','off');
    tic
    parfor s=1:scale
                ons = (rand(1,m) > tolProb);
                ons = ons - not(ons);
                ons = ons./sqrt(scale);   
                %Lambda_P
                [Z, flag]= pcg(L,(ons*(Fdelta))',tol,numIterations,pfun_);
            	output=output+(((Z(elist(:,1))-Z(elist(:,2)))).*(((Z(elist(:,3))-Z(elist(:,4))))));
            	Z=[]; 
                %Lambda_B
                ons = (rand(1,NT) > tolProb);
                ons = ons - not(ons);
                ons = ons./sqrt(scale); 
                ons=ons-mean(ons);
                [ZB, flag]= pcg(L,(ons*F)',tol,numIterations,pfun_);
                ons=[];
                ZB=ZB-sum(F*ZB)/(size(F,1));
                output_B=output_B+(((ZB(elist(:,1))-ZB(elist(:,2)))).*(((ZB(elist(:,3))-ZB(elist(:,4))))));
    end  
    disp('Time to Compute Hat Matrices (JL)')
    toc
end
%% Finding leverages and cross products using naive but exact parfor.
if strcmp(type_,'exact')
    tic  
        parfor m=1:M
            Bleft=sparse(elist(m,1:2),1,[-1 1],J,1);
            Bright=sparse(elist(m,3:4),1,[-1 1],J,1);
            [xtilde flag]=pcg(L,Bright,tol,numIterations,pfun_);
            %Lambda_P
            output(m)=Bleft'*xtilde;
            %Lambda_B
            [xtilde_left flag]=pcg(L,Bleft,tol,numIterations,pfun_);
            aux=cov(F*xtilde_left,F*xtilde)*(size(F,1)-1);
            output_B(m)=aux(1,2);
        end
    disp('Time to Compute Leverage Scores (Exact)')
    toc
end
%% Export
Pii=output;
Bii=output_B;
end 
   
    
