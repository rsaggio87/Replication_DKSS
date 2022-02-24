function [theta EIG_NORM max_x1bar_sq] = kss_quadratic_form(y,id,firmid,leave_out_level,controls,eigen_diagno,eigen_fast,type_of_algorithm,A_right,A_left,deMean,epsilon,Q)

theta=0;
if deMean== 1 & nargin == 13
	error('The user specified to demean the resulting quadratic form but included an inner matrix Q, that is not consistent. Set demean =0 and specify Q is the user wants to user an inner matrix Q ')	
		
	
end

if deMean== 0 & nargin == 12
	error('Need to specify Q. If beta`Aleft*(A_right*beta) is the right quadratic form, then set Q=identity')	
end


%Define key objects
[~,~,n]					= unique(firmid);
firmid					= n;
[~,~,n]					= unique(id);
id						= n; 
NT						= size(y,1);
D						= sparse(1:NT,id',1);
N						= size(D,2);
F						= sparse(1:NT,firmid',1);
J						= size(F,2);
no_controls				= 1.*(size(controls,2)==0)+ 0.*(size(controls,2)>0);


gcs 					= [NaN; id(1:end-1)];
gcs 					= id~=gcs;
lagfirmid				= [NaN; firmid(1:end-1)];
lagfirmid(gcs==1)		= NaN; %%first obs for each worker
stayer					= (firmid==lagfirmid);
stayer(gcs==1)			= 1;
stayer					= accumarray(id,stayer);
T						= accumarray(id,1);
stayer					= T==stayer;
movers					= stayer~=1;
movers					= movers(id);
T						= T(id);

if deMean == 1
Q=speye(N+J);
end

%Step 1: Residualize controls if controls are specified
if no_controls == 0
   			S			= speye(J-1);
   			S			= [S;sparse(-zeros(1,J-1))];  %N+JxN+J-1 restriction matrix 
   			X			= [D,F*S,controls];
   			xx			= X'*X;
   			Lchol		= ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
   			xy			= X'*y;
   			b			= pcg(xx,xy,1e-10,1000,Lchol,Lchol');
   			y			= y-X(:,N+J:end)*b(N+J:end);
   			no_controls	= 1; 
end

K						= 0;
S						= speye(J-1);
S						= [S;sparse(-zeros(1,J-1))];
X						= [D,F*S];
xx						= X'*X;
Lchol					= ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));

A_right_orig			= A_right;
A_left_orig				= A_left;

aux						= [speye(N) sparse(N,J-1); sparse(J,N) S];
A_right					= A_right*aux;
A_left					= A_left*aux;


%Step 2: Diagnostics
if eigen_diagno == 1 & eigen_fast == 0
			tic    
			disp('Calculating Sum of Squared Eigenvalues via Hutchinson...')
			%Compute the sum of squared eigenvalues associated to Atilde via
			%Hutchinson trick to estimate the trace of large matrices.
			NSIM		= 1000;
			trace		= zeros(NSIM,1);
			
			parfor 	   	s =1:NSIM  
			%inversion    
			xsimul=rand(size(X,1),1);
			xsimul=1.*(xsimul>0.5)-1.*(xsimul<=0.5);
			[coeff flag]=pcg(xx,X'*xsimul,1e-10,1000,Lchol,Lchol');


			%Variance of firm effects
			eff=A_right*coeff;
			if deMean == 1
			A_b=A_left'*(eff-mean(eff));			
			end
			if deMean == 0
			A_b=A_left'*Q*eff;
			end
			[aux, flag]=pcg(xx,A_b,1e-10,1000,Lchol,Lchol');
			aux=X*aux;
			trace(s)=aux'*aux;
			end
			
			%Results
			SUM_EIG		= mean(trace);
			disp('done')

			%Auxiliary for spectral decomposition
			if isequal(A_left,A_right) == 0
				a		= [A_right];
				b		= [A_left];
				abar	= sum(a,1)'/(sqrt(size(a,1)));
				bbar	= sum(b,1)'/(sqrt(size(b,1))); 
				ainv	= pcg(xx,abar,1e-5,1000,Lchol,Lchol');
				Sxxdot	= [xx sparse(N+J-1,1) abar+bbar; sparse(1,N+J-1) 1 1; sparse(1,N+J-1) 0 1];
				Adot	= 0.5*[S_ab -S_ab*ainv abar+bbar; bbar' -bbar'*ainv 1; abar' -abar'*ainv 1];
				[v, lambda_eig] = eigs(Adot,Sxxdot,3);
				Q=v(1:end-2,1)-ainv*v(end-1,1); %q=1
			end
			
			if isequal(A_left,A_right) == 1
				Sxxdot	= [xx sparse(N+J-1,1); sparse(1,N+J-1) 1];
				a		= A_right; 		
				abar	= sum(a,1)'/(sqrt(size(a,1)));
				S_aa	= a'*a;
				ainv	= pcg(xx,abar,1e-5,1000,Lchol,Lchol');
				Adot	= [S_aa -S_aa*ainv;abar' -abar'*ainv];
				[v, lambda_eig] = eigs(Adot,Sxxdot,3);
				Q		= v(1:end-1,1)-ainv*v(end,1); %q=1
				clear a abar S_aa ainv Adot v 
			end

			%Finalize spectral decomposition
    		lambda_eig	= diag(lambda_eig);
   			lambda_1	= lambda_eig(1)
    		[EIG_NORM x1bar_all]= eig_x1bar(X,Q,lambda_eig,SUM_EIG); %Note: What we label as x1bar in the code corresponds to w_iq in KSS.
    		EIG_NORM(1:3)
    		max_x1bar_sq=max(x1bar_all.^2)
    		disp('checking, must report 1')
			sum(x1bar_all.^2,1)
			clear v Adot Sxxdot 		          
		
end

if eigen_diagno == 1 & eigen_fast == 1

	 if isequal(A_left,A_right) == 0
				a		= [A_right];
				b		= [A_left];
				abar	= sum(a,1)'/(sqrt(size(a,1)));
				bbar	= sum(b,1)'/(sqrt(size(b,1))); 
				ainv	= pcg(xx,abar,1e-5,1000,Lchol,Lchol');
				Sxxdot	= [xx sparse(N+J-1,1) abar+bbar; sparse(1,N+J-1) 1 1; sparse(1,N+J-1) 0 1];
				Adot	= 0.5*[S_ab -S_ab*ainv abar+bbar; bbar' -bbar'*ainv 1; abar' -abar'*ainv 1];
				[v, lambda_eig] = eigs(Adot,Sxxdot,3);
				Q=v(1:end-2,1)-ainv*v(end-1,1); %q=1
	end
			
	if isequal(A_left,A_right) == 1
				Sxxdot	= [xx sparse(N+J-1,1); sparse(1,N+J-1) 1];
				a		= A_right; 		
				abar	= sum(a,1)'/(sqrt(size(a,1)));
				S_aa	= a'*a;
				ainv	= pcg(xx,abar,1e-5,1000,Lchol,Lchol');
				Adot	= [S_aa -S_aa*ainv;abar' -abar'*ainv];
				[v, lambda_eig] = eigs(Adot,Sxxdot,30);
				Q		= v(1:end-1,1)-ainv*v(end,1); %q=1
				clear a abar S_aa ainv Adot v 
	end

	%Finalize spectral decomposition
	lambda_eig			= diag(lambda_eig);
	SUM_EIG				= sum(lambda_eig.^2);
	[EIG_NORM x1bar_all]= eig_x1bar(X,Q,lambda_eig,SUM_EIG); %Note: What we label as x1bar in the code corresponds to w_iq in KSS.
	EIG_NORM(1:3)
	max_x1bar_sq=max(x1bar_all.^2)
	disp('checking, must report 1')
	sum(x1bar_all.^2,1)
	clear v Adot Sxxdot 	

end

if 1 == 1

%Step 3: Obtain (Bii,Pii)
X						= [D,-F];
xx						= X'*X;
disp('Building preconditioner for Laplacian Matrix...')
Lchol 					= cmg_sdd(xx); %preconditioner for Laplacian matrices.
A_right					= A_right_orig;
A_left					= A_left_orig;

[~,~,match_id]			= unique([id firmid],'rows','stable');
	
if strcmp(leave_out_level,'obs')
		index   		= (1:NT)';    
		clustering_var  = index;
end
	
if strcmp(leave_out_level,'matches')
		clustering_var	= match_id;
end
	
elist					= index_constr(clustering_var,id,match_id);

[Lambda_P, Lambda_B] 	= eff_res_inputA(A_right,A_left,deMean,Q,X,xx,Lchol,N,J,K,elist,leave_out_level,movers,T,type_of_algorithm,id,firmid,match_id,epsilon);

%Step 4: Get theta
X						= [D,F*S];
xx						= X'*X;
Lchol					= ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
aux						= [speye(N) sparse(N,J-1); sparse(J,N) S];
A_right					= A_right*aux;
A_left					= A_left*aux;
xy						= X'*y;
b						= pcg(xx,xy,1e-10,1000,Lchol,Lchol');
xb						= X*b;
eta						= y-xb;

if deMean == 1
	theta_bias			= A_right*b ;
	theta_bias			= theta_bias-mean(theta_bias);
	theta_bias		    = (A_left*b)'*theta_bias;
end

if deMean == 0
	theta_bias			= (A_left*b)'*Q*A_right*b ;
end

I_Lambda_P=(speye(NT,NT)-Lambda_P);
eta_h=I_Lambda_P\eta; %Leave one out residual
msg = lastwarn ; 
if contains(msg, 'singular') 
            	s=['******************************************'];
                disp(s);
                disp(s); 
				disp('Warning: OLS coefficient not always identified when leaving a particular set of observation out as specified by "leave_out_level"')
				disp('One example where this occurs is when the user asks to run leave out on matches without restricting the analysis to movers only.')
				s=['******************************************'];
                disp(s);
                disp(s); 
				
end

%KSS formula
theta 					= theta_bias-y'*Lambda_B*eta_h;


%Step 5: Output results
%% STEP 6: REPORTING
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s);
disp(s);
s=['Maximum Leverage: ' num2str(max(diag(Lambda_P)))];
disp(s)
s=['Theta - Plug in: ' num2str(theta_bias)];
disp(s)
s=['Theta - KSS: ' num2str(theta)];
disp(s)
if eigen_diagno==1  
        s=['*********************Diagnostics on the Quadratic Form*********************'];
        disp(s);
        s=['ratio of eigenvalues: '];
        disp(s)
        EIG_NORM(1:3)
        s=['Weak Lindeberg Condition: ' num2str(max_x1bar_sq)];
        disp(s)
        s=['Sum of squared eigenvalues: ' num2str(SUM_EIG)];
        disp(s)
       
end

end
end

