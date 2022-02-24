function [R2, R2_KSS,Pii_R2,theta,sigma_i,xb]    = estimate_rakim_levels(model,y,id,firmid,lagfirmid,firmid_orig,lagfirmid_orig,controls,identified,dual_list,scale,STRINGA);	

%Auxiliaries
	numIterations			= 10000;
	tol						= 1e-10;
	tolProb					= 0.5;
	theta					= zeros(9,2);

	
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
	if model == 1
	X				= [D Flag controls];
	K				= size(X,2);
	end
	if model == 2
	X				= [D F controls];
	K				= size(X,2);
	end
	
	
%Now the matrices to do the variance decomposition (job-year)
	sel						= identified;
	Niden					= sum(identified);

if model == 3	
	Dvar					= [D(sel,:)  sparse(Niden,K-N)];
	Fvar					= [sparse(Niden,N) 	 F(sel,:)    sparse(Niden,K-N-J)];
	Flagvar					= [sparse(Niden,N+J) Flag(sel,:) sparse(Niden,K-N-J-Jlag)];
	Wvar					= [sparse(Niden,N+J+Jlag) controls(sel,:)];	
end

if model == 2	
	Dvar					= [D(sel,:)  sparse(Niden,K-N)];
	Fvar					= [sparse(Niden,N) 	 F(sel,:)    sparse(Niden,K-N-J)];
	Flagvar					= [sparse(Niden,N+J) controls(sel,:)]; %controls now
	Wvar					= [sparse(Niden,N+J) controls(sel,:)];	
end

if model == 1	
	Dvar					= [D(sel,:)  sparse(Niden,K-N)]; %won't be read
	Fvar					= [D(sel,:)  sparse(Niden,K-N)]; %won't be read
	Flagvar					= [D(sel,:)  sparse(Niden,K-N)]; %won't be read
	Wvar					= [D(sel,:)  sparse(Niden,K-N)]; %won't be read
end
	
%Now the matrices to do the variance decomposition (firms that are dual connected)
if model == 3
	N_dual					= size(dual_list,1)
	sel						= dual_list(:,2)';	
	DUAL_F					= [sparse(N_dual,N) 	 sparse((1:N_dual)',sel,1,N_dual,J)   	  sparse(N_dual,K-N-J)];	
	DUAL_F					= repelem(DUAL_F,dual_list(:,5),1); %span it w.r.t. size of each firm
	sel						= dual_list(:,3)';
	DUAL_Flag				= [sparse(N_dual,N+J) 	 sparse((1:N_dual)',sel,1,N_dual,Jlag)    sparse(N_dual,K-N-J-Jlag)];	
	DUAL_Flag				= repelem(DUAL_Flag,dual_list(:,5),1); %span it w.r.t. size of each firm	
end

if model == 2 | model ==1 %just to make sure code runs
	DUAL_F					= Fvar;
	DUAL_Flag				= Flagvar;
end
	
%Estimate the model
	xy						= X'*y;
	xx						= X'*X;
	L						= ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
	tic
	b						= pcg(xx,xy,tol,numIterations,L,L');	
	toc
	xb						= X*b;		

%Plug-in R2
	ESS					    = var(xb);
	TSS					    = var(y);
	R2					    = ESS/TSS
	adjR2					= 1-((1-R2)*(NT-1))/(NT-size(X,2)-1)
	eta					    = y-xb;
	
%Export the auxiliary file for Stata
if model == 3
	out 					= [y(identified),id(identified),firmid_orig(identified),lagfirmid_orig(identified),Dvar*b,Fvar*b,Flagvar*b,Wvar*b,full(controls(identified,:))];
	s						= ['../build/src/RAKIM/' STRINGA '.csv']
	dlmwrite(s, out, 'precision', 16); 
	out						= dual_list;
	s						= ['../build/src/RAKIM/DUAL_LIST' STRINGA '.csv']
	dlmwrite(s, out, 'precision', 16); 
end

if model == 2
	fe						= b(N+1:N+J);
	fe						= F*fe;
	out 					= [y(identified),id(identified),firmid_orig(identified),lagfirmid_orig(identified),fe];
	s						= ['../build/src/RAKIM/AKM' STRINGA '.csv']
	dlmwrite(s, out, 'precision', 16); 
end	

if model == 1
	fe						= b(N+1:N+Jlag);
	fe						= Flag*fe;
	out 					= [y(identified),id(identified),firmid_orig(identified),lagfirmid_orig(identified),fe];
	s						= ['../build/src/RAKIM/LAKM' STRINGA '.csv']
	dlmwrite(s, out, 'precision', 16); 
end	
	

if model >1	
%PI variance components
   for parametro = 1:9
		for metodo 	  = 1:1
		
			if parametro == 1 %variance of person effects
			left		= Dvar*b;
			right		= Dvar*b;
			end

			if parametro == 2 %variance of firm effects
			left		= Fvar*b;
			right		= Fvar*b;
			end
			
			if parametro == 3 %variance of firm lag effects
			left		= Flagvar*b;
			right		= Flagvar*b;
			end

			if parametro == 4 %covariance of person,firm effects
			left		= Dvar*b;
			right		= Fvar*b;
			end

			if parametro == 5 %%covariance of person,firm lagged effects
			left		= Dvar*b;
			right		= Flagvar*b;
			end
			
			if parametro == 6 %%covariance of firm,firm lagged effects
			left		= Fvar*b;
			right		= Flagvar*b;
			end
			
			if parametro == 7 %%Variance of contemporaneous firm effects (dual connected)
			left		= DUAL_F*b;
			right		= DUAL_F*b;
			end
			
			if parametro == 8 %%Variance of lagged firm effects (dual connected)
			left		= DUAL_Flag*b;
			right		= DUAL_Flag*b;
			end
			
			if parametro == 9 %%Covariance of lagged firm and firm effects (dual connected)
			left		= DUAL_Flag*b;
			right		= DUAL_F*b;
			end

			COV		 	= cov(left,right);
			theta(parametro,metodo)=COV(1,2);
			
		end
	end	
theta(:,1)	
end

	
%KSS-Correction	
	numIterations			= 1000;
	tol						= 1e-5;
	tolProb					= 0.5;

	%# of draws
	disp('# of Simulated Projections for JLL:')
	scale
		
	%set up the matrices
	tic	
	N						= size(X,1);
	NSEL					= size(Dvar,1);
	NDUAL					= size(DUAL_Flag,1);
	xx_c					= parallel.pool.Constant(xx);
	X						= parallel.pool.Constant(X);
	L						= parallel.pool.Constant(L);
	
	Dvar					= parallel.pool.Constant(Dvar);
	Fvar					= parallel.pool.Constant(Fvar);
	Flagvar					= parallel.pool.Constant(Flagvar);
	
	DUAL_Flag				= parallel.pool.Constant(DUAL_Flag);
	DUAL_F					= parallel.pool.Constant(DUAL_F);

	disp('time to send the matrices to each core')
	toc
	
	%preallocate
	Pii						= zeros(N,1);
	Pii_R2					= zeros(N,1);
	Bii						= zeros(N,9);
	tic
	
	parfor i=1:scale
				
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Pii				
				%Random Projection Matrix (Rademacher)
                ons 	= (rand(1,N) > tolProb);
                ons 	= ons - not(ons);
                ons 	= ons./sqrt(scale);
                
           
                %Get me the row
                [Z, flag]= pcg(xx_c.Value,(ons*(X.Value))',tol,numIterations,L.Value,(L.Value)');
                
                %Collect
                Z		 = X.Value*Z;
                Z		 = Z.^2;
                Pii		 = Pii+Z;
                
                %Mean adjusted Pii
				ons		 = ons-mean(ons);
				[Z, flag]= pcg(xx_c.Value,(ons*(X.Value))',tol,numIterations,L.Value,(L.Value)')
				Z		 = X.Value*Z;
                Z		 = Z.^2;
                Pii_R2	 = Pii_R2+Z;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Job-Year Decomposition	
				
				if model > 1
				Bii_aux	= zeros(N,9);               
                ons 	= (rand(1,NSEL) > tolProb);
                ons 	= ons - not(ons);
                ons 	= ons./sqrt(scale);
                ons 	= ons-mean(ons);
                Z		= sparse(N,3);
                
                for pp = 1:3
                
					if pp == 1
					design	=  Dvar.Value;
					end
				
					if pp == 2
					design	=  Fvar.Value;
					end
				
					if pp == 3 %if model =2, then this is capturing the variance components associated with controls.
					design	=  Flagvar.Value;
					end
				
					[aux flag]		= pcg(xx_c.Value,(ons*(design))',tol,numIterations,L.Value,(L.Value)');
					Z(:,pp) = X.Value*aux;
                end
                design		= [];
                 
        		
        		for pp = 1:6
        		
        			if pp == 1 %variance of person effects
					Z_left	=  Z(:,1);
					Z_right	=  Z(:,1);
					end
				
					if pp == 2 %variance of firm effects
					Z_left	=  Z(:,2);
					Z_right	=  Z(:,2);
					end
				
					if pp == 3 %variance of lag firm effects
					Z_left	=  Z(:,3);
					Z_right	=  Z(:,3);
					end
					
					if pp == 4 %covariance of person,firm effects
					Z_left	=  Z(:,1);
					Z_right	=  Z(:,2);
					end
					
					if pp == 5 %covariance of person, lag firm effects
					Z_left	=  Z(:,1);
					Z_right	=  Z(:,3);
					end
					
					if pp == 6 %covariance of firm, lag firm effects
					Z_left	=  Z(:,2);
					Z_right	=  Z(:,3);
					end
					
					Bii_aux(:,pp) = (Z_left.*Z_right)/NSEL;
        		
        		end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Dual Decomposition
 				if 		model == 3       		
        		ons 	= (rand(1,NDUAL) > tolProb);
                ons 	= ons - not(ons);
                ons 	= ons./sqrt(scale);
                ons 	= ons-mean(ons);
                Z		= sparse(N,2);
                
                for pp = 1:2
                
					if pp == 1
					design	=  DUAL_F.Value;
					end
				
					if pp == 2
					design	=  DUAL_Flag.Value;
					end
				
					[aux flag]	= pcg(xx_c.Value,(ons*(design))',tol,numIterations,L.Value,(L.Value)');
					Z(:,pp) = X.Value*aux;
                end
                design		= [];
                 
        		
        		for pp = 7:9
        		
        			if pp == 7 %variance of firm effects (dual)
					Z_left	=  Z(:,1);
					Z_right	=  Z(:,1);
					end
				
					if pp == 8 %variance of lag firm effects (dual)
					Z_left	=  Z(:,2);
					Z_right	=  Z(:,2);
					end
				
					if pp == 9 %covariance of firm, lag firm effects (dual)
					Z_left	=  Z(:,2);
					Z_right	=  Z(:,1);
					end
					
					Bii_aux(:,pp) = (Z_left.*Z_right)/NDUAL;
        		
        		end
        		end
        		
        	Bii= Bii + Bii_aux;
        		
        	end	
                          
    end 
    toc
    disp('Time to perform JLL')         

%Display max(Pii)
	disp('Max-Pii')  
	max(Pii)
	disp('90th,95th,99th quantile of Pii')  
	quantile(Pii,[0.90 0.95 0.99])	

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

if model >1
%KSS variance components
   for parametro = 1:9
		for metodo 	  = 2:2
			theta(parametro,metodo)=theta(parametro,1)-sum(Bii(:,parametro).*sigma_i);	
		end
   end
end
theta		
				
end
