function [R2, R2_KSS,Pii_R2,theta,sigma_i] = estimate_rakim(X,y,A_diff_network1,A_diff_network2,A_levels_network1,A_levels_network2,type_of_algorithm,STRINGA,scale);

%Auxiliaries
	numIterations			= 10000;
	tol						= 1e-6;
	if nargin <= 7
	scale					= 1000;
	end
	tolProb					= 0.5;
	variance_decomp			= 0;
	
	if nargout 	> 3
	variance_decomp = 1;
	end

%First step
	xx					    = X'*X;
	K						= size(xx,2);	
	
%Estimate the model
	xy						= X'*y;
	b						= pcg(xx,xy,tol,numIterations);	
	xb						= X*b;	

%Plug-in R2
	ESS					    = var(xb);
	TSS					    = var(y);
	R2					    = ESS/TSS;
	eta					    = y-xb;
	dof						= size(y,1)-size(X,2)-1;
	NT						= size(y,1);
	R2						= 1-(1-R2)*((NT-1)/dof);
	
%Plug-in variance decompositions		
if variance_decomp == 1

	for parametro = 1:6
		for metodo 	  = 1:1

			if parametro == 1
			left		= A_diff_network1*b;
			right		= A_diff_network1*b;
			end

			if parametro == 2
			left		= A_diff_network2*b;
			right		= A_diff_network2*b;
			end
			
			if parametro == 3
			left		= A_diff_network2*b;
			right		= A_diff_network1*b;
			end
			
			if parametro == 4
			left		= A_levels_network1*b;
			right		= A_levels_network1*b;
			end

			if parametro == 5
			left		= A_levels_network2*b;
			right		= A_levels_network2*b;
			end
			
			if parametro == 6
			left		= A_levels_network2*b;
			right		= A_levels_network1*b;
			s			= ['tables/RAKIM_RAW_' STRINGA '.csv'];
			out			= [left right];
			dlmwrite(s, out, 'delimiter', '\t', 'precision', 16); 
			end

			COV		 	= cov(left,right);
			theta(parametro,metodo)=COV(1,2);
			
		end
	end		
end	

%Calculate Pii
if strcmp(type_of_algorithm,'exact')
	
	% Count the unique transitions
	[new_X,bbb,index_X]		= unique(X,'rows','stable');
	N						= size(new_X,1);
	
	% Precreate
	Pii						= zeros(N,1);
	Bii_net1				= zeros(N,1);
	Bii_net2				= zeros(N,1);
	Bii_net12				= zeros(N,1);
	Bii_dual_net1			= zeros(N,1);
	Bii_dual_net2			= zeros(N,1);
	Bii_dual_net12			= zeros(N,1);
	
	
	% Set in parallel
	xx						= parallel.pool.Constant(xx);
	new_X					= parallel.pool.Constant(new_X);
	A_diff_network1			= parallel.pool.Constant(A_diff_network1);
	A_diff_network2			= parallel.pool.Constant(A_diff_network2);
	A_levels_network1		= parallel.pool.Constant(A_levels_network1);
	A_levels_network2		= parallel.pool.Constant(A_levels_network2);

	parfor	 i=1:N
	[aux flag]	  	 		= pcg(xx.Value,new_X.Value(i,:)',tol,numIterations);
	
	%Pii
	aux						= aux-mean(aux);
	Pii(i)					= new_X.Value(i,:)*aux;
	
		if variance_decomp == 1
			
			%Bii in FD
			net1			= A_diff_network1.Value*aux;
			net2			= A_diff_network2.Value*aux;
			COV				= cov(net1,net2);
			Bii_net1(i)		= COV(1,1);
			Bii_net2(i)		= COV(2,2);
			Bii_net12(i)	= COV(1,2);
	
			%Bii in Levels
			net1			= A_levels_network1.Value*aux;
			net2			= A_levels_network2.Value*aux;
			COV				= cov(net1,net2);
			Bii_dual_net1(i)= COV(1,1);
			Bii_dual_net2(i)= COV(2,2);
			Bii_dual_net12(i)= COV(1,2);
			
		end

	end
	
	% Collect
	Pii						= Pii(index_X);
	if variance_decomp == 1
	Bii_net1				= Bii_net1(index_X);
	Bii_net2				= Bii_net2(index_X);
	Bii_net12				= Bii_net12(index_X);
	Bii_dual_net1			= Bii_dual_net1(index_X);
	Bii_dual_net2			= Bii_dual_net2(index_X);
	Bii_dual_net12			= Bii_dual_net12(index_X);
	end
end

if strcmp(type_of_algorithm,'JLL')
	
	%# of draws
	disp('# of Simulated Projections for JLL:')
	scale
		
	%set up the matrices
	tic
	N						= size(X,1);
	NSEL					= size(A_diff_network1,1);
	NBAR					= size(A_levels_network1,1);
	xx_c					= parallel.pool.Constant(xx);
	X						= parallel.pool.Constant(X);
	A_diff_network1			= parallel.pool.Constant(A_diff_network1);
	A_diff_network2			= parallel.pool.Constant(A_diff_network2);
	A_levels_network1		= parallel.pool.Constant(A_levels_network1);
	A_levels_network2		= parallel.pool.Constant(A_levels_network2);
	disp('time to send the matrices to each core')
	toc
	
	%preallocate
	Pii						= zeros(N,1);
	Pii_R2					= zeros(N,1);
	Bii_net1				= zeros(N,1);
	Bii_net2				= zeros(N,1);
	Bii_net12				= zeros(N,1);
	Bii_dual_net1			= zeros(N,1);
	Bii_dual_net2			= zeros(N,1);
	Bii_dual_net12			= zeros(N,1);
	
	parfor i=1:scale
				
				%Random Projection Matrix (Rademacher)
                ons = (rand(1,N) > tolProb);
                ons = ons - not(ons);
                ons = ons./sqrt(scale);
                
           
                %Get me the row
                [Z, flag]= pcg(xx_c.Value,(ons*(X.Value))',tol,numIterations)
                
                %Collect
                Z		 = X.Value*Z;
                Z		 = Z.^2;
                Pii		 = Pii+Z;
				
				%Mean adjusted Pii
				ons		 = ons-mean(ons);
				[Z, flag]= pcg(xx_c.Value,(ons*(X.Value))',tol,numIterations)
				Z		 = X.Value*Z;
                Z		 = Z.^2;
                Pii_R2	 = Pii_R2+Z;
				
				if  variance_decomp == 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
                
                %Random Projection Variance decomposition in differences
                ons = (rand(1,NSEL) > tolProb);
                ons = ons - not(ons);
                ons = ons./sqrt(scale);
                ons = ons-mean(ons);
                
                %Get me the row
                [Z_1, flag]= pcg(xx_c.Value,(ons*(A_diff_network1.Value))',tol,numIterations)
                [Z_2, flag]= pcg(xx_c.Value,(ons*(A_diff_network2.Value))',tol,numIterations)
                
                %Collect
                Z_1	 	 = X.Value*Z_1;
                Z_2 	 = X.Value*Z_2;
                Bii_net1 = Bii_net1+(Z_1.^2)/NSEL;
                Bii_net2 = Bii_net2+(Z_2.^2)/NSEL;
                Bii_net12= Bii_net12+(Z_1.*Z_2)/NSEL;
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
                
                %Random Projection Variance decomposition in differences
                ons = (rand(1,NBAR) > tolProb);
                ons = ons - not(ons);
                ons = ons./sqrt(scale);
                ons = ons-mean(ons);
                
                %Get me the row
                [Z_1, flag]= pcg(xx_c.Value,(ons*(A_levels_network1.Value))',tol,numIterations)
                [Z_2, flag]= pcg(xx_c.Value,(ons*(A_levels_network2.Value))',tol,numIterations)
                
                %Collect
                Z_1	 	 		= X.Value*Z_1;
                Z_2 	 		= X.Value*Z_2;
                Bii_dual_net1 	= Bii_dual_net1+(Z_1.^2)/NBAR;
                Bii_dual_net2 	= Bii_dual_net2+(Z_2.^2)/NBAR;
                Bii_dual_net12	= Bii_dual_net12+(Z_1.*Z_2)/NBAR;
                end
	end
end	

%Censor
	sel						= Pii>=0.99;
	Pii(sel)				= 0.99;		

%KSS Estimated variances
	invP				    = 1-Pii;
	eta_h				    = eta./invP;
	sigma_i					= y.*eta_h;
	
%KSS R2	
	R2_KSS				    = ESS - mean(Pii_R2.*sigma_i);
	R2_KSS					= R2_KSS/TSS;
	
%Now compute theta = [var(FDELTA1*psi1); var(FDELTA2*psi2); cov(FDELTA1*psi1,FDELTA2*psi2); var(FBAR1*beta); var(FBAR2*beta); cov(FBAR1*beta,FBAR2*beta) ]

if variance_decomp == 1

	for parametro = 1:6
		for metodo 	  = 2:2

			if parametro == 1
			Bii			= Bii_net1;
			end

			if parametro == 2
			Bii			= Bii_net2;
			end
			
			if parametro == 3
			Bii			= Bii_net12;
			end
			
			if parametro == 4
			Bii			= Bii_dual_net1;
			end

			if parametro == 5
			Bii			= Bii_dual_net2;
			end
			
			if parametro == 6
			Bii			= Bii_dual_net12;
			end
			
			theta(parametro,metodo)=theta(parametro,1)-sum(Bii.*sigma_i);
			
		end
	end		
end		
if  variance_decomp == 1	
%Check how far we are from full rank
	disp('Sum of Pii --- PVR MODEL')
	sum(Pii)	
	disp('Size of xx --- PVR MODEL')
	K
end	
	
end	
