function theta    = estimate_rakim_sector(y,sigma_i,id,firmid,lagfirmid,firmid_orig,lagfirmid_orig,controls,identified,dual_list,scale,STRINGA);	

%rescale sector
	dual_list(:,6)=dual_list(:,6)+1;

%Auxiliaries
	numIterations	= 10000;
	tol				= 1e-10;
	tolProb			= 0.5;
	Nsectors		= max(dual_list(:,6));
	theta			= zeros(Nsectors,3,2);

	
%Setting up the matrices	
	NT				= size(id,1);
	J				= max(firmid);
	Jlag			= max(lagfirmid);
	N				= max(id);

	F				= sparse((1:NT)',firmid',1,NT,J);
	Flag			= sparse((1:NT)',lagfirmid',1,NT,Jlag);	
	D				= sparse((1:NT)',id',1,NT,N);

%Which model are we estimating?
	X				= [D F Flag controls];
	K				= size(X,2);
			
%Estimate the model
	xy				= X'*y;
	xx				= X'*X;
	L				= ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
	tic
	b				= pcg(xx,xy,tol,numIterations,L,L');	
	toc
	xb				= X*b;		

	
%PI variance components
for metodo 	  = 1:1
   for settore = 1:Nsectors
   
				%Create the Key matrices for a given configuration.
				[DUAL_F DUAL_Flag]	= pick_matrix(dual_list,settore,N,J,Jlag,K);	

   	   				for parametro = 1:3 	
							if parametro == 1 %%Variance of contemporaneous firm effects (dual connected)
										left		= DUAL_F*b;
										right		= DUAL_F*b;
							end
			
							if parametro == 2 %%Variance of lagged firm effects (dual connected)
										left		= DUAL_Flag*b;
										right		= DUAL_Flag*b;
							end
			
							if parametro == 3 %%Covariance of lagged firm and firm effects (dual connected)
										left		= DUAL_Flag*b;
										right		= DUAL_F*b;
							end
						
							COV		 						= cov(left,right);
					 		theta(settore,parametro,metodo)	= COV(1,2);
					
					end				
	end		
end	
theta(:,:,1)
%SAVE OUTPUT
out						= theta(:,:,1);
s						= ['../build/src/RAKIM/SORKIN_GRAPH_PI.csv']
dlmwrite(s, out, 'precision', 16); 
			
%KSS-Correction	
	numIterations			= 1000;
	tol						= 1e-5;
	tolProb					= 0.5;

	%# of draws
	disp('# of Simulated Projections for JLL:')
	scale
		
	%set up the matrices
	tic	
	xx_c					= parallel.pool.Constant(xx);
	X						= parallel.pool.Constant(X);
	L						= parallel.pool.Constant(L);
	DUALE					= parallel.pool.Constant(dual_list);

	disp('time to send the matrices to each core')
	toc
	
	%preallocate
	correzione			    = zeros(Nsectors,3);
	tic
	
	parfor i=1:scale
				correzione_aux=zeros(Nsectors,3);
				for settore = 1:Nsectors
						
					%Build the associated matrices
					[DUAL_F DUAL_Flag]		= pick_matrix(DUALE.Value,settore,N,J,Jlag,K);
					Z						= sparse(NT,2);
						
					%Rademacher
					NDUAL					= size(DUAL_F,1);
					ons 					= (rand(1,NDUAL) > tolProb);
					ons 					= ons - not(ons);
                	ons 					= ons./sqrt(scale);
                	ons 					= ons-mean(ons);
                	
                	%Inversion
                	for parametro = 1:2
                
					if parametro == 1
					design	=  DUAL_F;
					end
				
					if parametro == 2
					design	=  DUAL_Flag;
					end
				
					[aux flag]	= pcg(xx_c.Value,(ons*(design))',tol,numIterations,L.Value,(L.Value)');
					Z(:,parametro) = X.Value*aux;
					design		= [];
					end
					
					%Complete the three-parameters for a given sector
					
					for parametro = 1:3
        		
        			if parametro == 1 %variance of firm effects (dual)
					Z_left	=  Z(:,1);
					Z_right	=  Z(:,1);
					end
				
					if parametro == 2 %variance of lag firm effects (dual)
					Z_left	=  Z(:,2);
					Z_right	=  Z(:,2);
					end
				
					if parametro == 3 %covariance of firm, lag firm effects (dual)
					Z_left	=  Z(:,2);
					Z_right	=  Z(:,1);
					end
					
					Bii					  			 = (Z_left.*Z_right)/NDUAL;
					correzione_aux(settore,parametro)= sum(Bii.*sigma_i);
					end
        		
        		end
					
		correzione = correzione + correzione_aux			
    end			
    toc
    disp('Time to perform JLL')         

%KSS variance components
for metodo 	  = 2:2
   for settore = 1:Nsectors
   		for parametro=1:3
   					theta(settore,parametro,metodo) = theta(settore,parametro,1)-correzione(settore,parametro);
   		end					
   end
end 

%SAVE OUTPUT
theta(:,:,2)
out						= theta(:,:,2);
s						= ['../build/src/RAKIM/SORKIN_GRAPH_KSS.csv']
dlmwrite(s, out, 'precision', 16); 	
				
end
