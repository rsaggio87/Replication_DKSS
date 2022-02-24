function [EIG1 EIG2 EIG3 sigma2_true theta stdTheta V_sym_SIM V_nosym_SIM sdV_sym_SIM sdV_nosym_SIM oracle_coverage coverage_SYM coverage_NOSYM]  = leave_out_FD_task15(y,id,firmid,leave_out_level,leave_2_out,controls,type_algorithm,eigen_diagno,eigen_fast,do_montecarlo,restrict_movers,filename)
s =[filename 'TUTTO_FD_NEW.mat'];
load(s)
disp('time to brute force C')
toc

%Auxiliary
movers_sel																  = firmid_delta~=firmid_delta_f; %focus on movers only
vec 																	  = C_build_FD(ydelta,Fdelta,L,pfun_,Lambda_B,I_Lambda_P,L_P,F);

%Standard errors
[V_sym V_nosym Q_share C_share_1 C_share_2 C_share_3 C_share_4 C_share_5 Chamard Dsym vec index_movers index_movers_i index_movers_2 index_movers_i_2 Ppath_1 Ppath_2 no_A2 list_problem taxonomy magical_number TAXONOMY_I TAXONOMY_J MOVERS_i MOVERS_j FOCUS_I FOCUS_J Ppath_PROBLEM_i Ppath_PROBLEM_j taxonomyALT TAXONOMY_ALT_I TAXONOMY_ALT_J MOVERS_ALT_i MOVERS_ALT_j FOCUS_ALT_I FOCUS_ALT_J Ppath_ALT_PROBLEM_i Ppath_ALT_PROBLEM_j]  = paths_network(id_movers,firmid_delta,firmid_delta_f,ydelta,C);

V_nosym 											 					  = sqrt(V_nosym/NT^2)
V_sym 																      = sqrt(V_sym/NT^2)


%[V_sym V_nosym Q_share C_share_1 C_share_2 C_share_3 C_share_4 C_share_5 Dsym vec index_movers index_movers_i index_movers_2 index_movers_i_2 Ppath_1 Ppath_2 no_A2 list_problem taxonomy magical_number TAXONOMY_I TAXONOMY_J MOVERS_i MOVERS_j FOCUS_I FOCUS_J Ppath_PROBLEM_i Ppath_PROBLEM_j taxonomyALT TAXONOMY_ALT_I TAXONOMY_ALT_J MOVERS_ALT_i MOVERS_ALT_j FOCUS_ALT_I FOCUS_ALT_J Ppath_ALT_PROBLEM_i Ppath_ALT_PROBLEM_j]  = paths_network_v2(id_movers,firmid_delta,firmid_delta_f,ydelta,vec,Fdelta,L,F,tol,epsilon,type_algorithm,pfun_,Lambda_P,Lambda_B);

%V_nosym 											 					  = sqrt(V_nosym/NT^2)
%V_sym 																      = sqrt(V_sym/NT^2)

%[CHECK1 CHECK2]	 													  = MC_fixed_FD_v2(ydelta,Fdelta,L,F,tol,epsilon,type_algorithm,pfun_,Lambda_P,Lambda_B,Dsym,vec,id_movers(movers_sel),movers_sel,index_movers,index_movers_i,index_movers_2,index_movers_i_2,Ppath_1,Ppath_2,no_A2,list_problem,taxonomy,magical_number,TAXONOMY_I,TAXONOMY_J,MOVERS_i,MOVERS_j,FOCUS_I, FOCUS_J,Ppath_PROBLEM_i,Ppath_PROBLEM_j,taxonomyALT,TAXONOMY_ALT_I,TAXONOMY_ALT_J,MOVERS_ALT_i,MOVERS_ALT_j,FOCUS_ALT_I, FOCUS_ALT_J,Ppath_ALT_PROBLEM_i,Ppath_ALT_PROBLEM_j);
%CHECK1 											 					  = sqrt(CHECK1/NT^2)
%CHECK2 																  = sqrt(CHECK2/NT^2)



if do_montecarlo == 1

%%%%%%%%%%%%%%%%%%%%%%Montecarlo
%fit the sigmas first
sigma_i 	= ydelta.*eta_h; 
SIZE		= diag(F'*F);
SIZE		= log(SIZE);
F1			= F(gcs==1,:);
F2			= F(gcs==2,:);
SIZE1 		= F1*SIZE;
SIZE2		= F2*SIZE;
xdata		= [ones(Ndelta,1) diag(Lambda_P) diag(Lambda_B) SIZE1 SIZE2];
fun 		= @(x,xdata)exp(x(1)*xdata(:,1)+x(2)*xdata(:,2)+x(3)*xdata(:,3)+x(4)*xdata(:,4) + x(5)*xdata(:,5));
x0 			= [0,0,0,0,0];
[b] 		= lsqcurvefit(fun,x0,xdata,sigma_i)
sigma_i		= fun(b,xdata);
figure
histogram(sigma_i(xdata(:,2)>0))
saveas(gcf,[ s 'HISTOGRAM.png'])

%Fill zeros
NSIM			= 1000;
theta			= zeros(NSIM,1);
V_sym_SIM		= zeros(NSIM,1);
V_nosym_SIM		= zeros(NSIM,1);
coverage_SYM	= zeros(NSIM,1);
coverage_NOSYM	= zeros(NSIM,1);
V_sym_rand_SIM	= zeros(NSIM,1);
V_nosym_rand_SIM= zeros(NSIM,1);
coverage_rand_SYM= zeros(NSIM,1);
coverage_rand_NOSYM= zeros(NSIM,1);

%Create the DGP
psi_true		= psi_hat_norm*sqrt(sigma2_psi/sigma_2_psi_AKM);
sigma2_true 	= var(F*psi_true);


%This is to smooth the parallel process
tic
C				= parallel.pool.Constant(C);
Chamard			= parallel.pool.Constant(Chamard);
Dsym			= parallel.pool.Constant(Dsym);
index_movers    = parallel.pool.Constant(index_movers);
index_movers_i  = parallel.pool.Constant(index_movers_i);
index_movers_2  = parallel.pool.Constant(index_movers_2);
index_movers_i_2= parallel.pool.Constant(index_movers_i_2);

Ppath_1			= parallel.pool.Constant(Ppath_1);
Ppath_2			= parallel.pool.Constant(Ppath_2);
no_A2			= parallel.pool.Constant(no_A2);
list_problem    = parallel.pool.Constant(list_problem);
taxonomy  		= parallel.pool.Constant(taxonomy);
magical_number  = parallel.pool.Constant(magical_number);


TAXONOMY_I		= parallel.pool.Constant(TAXONOMY_I);
TAXONOMY_J		= parallel.pool.Constant(TAXONOMY_J);
MOVERS_i		= parallel.pool.Constant(MOVERS_i);
MOVERS_j		= parallel.pool.Constant(MOVERS_j);
FOCUS_I			= parallel.pool.Constant(FOCUS_I);
FOCUS_J			= parallel.pool.Constant(FOCUS_J);
Ppath_PROBLEM_i	= parallel.pool.Constant(Ppath_PROBLEM_i);
Ppath_PROBLEM_j	= parallel.pool.Constant(Ppath_PROBLEM_j);


taxonomyALT		= parallel.pool.Constant(taxonomyALT);
TAXONOMY_ALT_I	= parallel.pool.Constant(TAXONOMY_ALT_I);
TAXONOMY_ALT_J	= parallel.pool.Constant(TAXONOMY_ALT_J);
MOVERS_ALT_i	= parallel.pool.Constant(MOVERS_ALT_i);
MOVERS_ALT_j	= parallel.pool.Constant(MOVERS_ALT_j);
FOCUS_ALT_I		= parallel.pool.Constant(FOCUS_ALT_I);
FOCUS_ALT_J		= parallel.pool.Constant(FOCUS_ALT_J);
Ppath_ALT_PROBLEM_i= parallel.pool.Constant(Ppath_ALT_PROBLEM_i);
Ppath_ALT_PROBLEM_j= parallel.pool.Constant(Ppath_ALT_PROBLEM_j);
toc

MSE             = 0.10*var(ydelta);

parfor s=1:NSIM
	ydelta						  = Fdelta*psi_true+sqrt(sigma_i).*randn(Ndelta,1);
	%ydelta						  = Fdelta*psi_true+sqrt(MSE).*randn(Ndelta,1);
	
	%vec						  = C.Value*ydelta;
	%vec						  = C*ydelta;
	vec 						  = C_build_FD(ydelta,Fdelta,L,pfun_,Lambda_B,I_Lambda_P,L_P,F);
	
	%Run KSS
	xy                        	  = Fdelta'*ydelta;
	[b flag]					  = pcg(L,xy,1e-10,1000,pfun_);
	r						      = ydelta-Fdelta*b;
	[eta_h, flag]			      = pcg(I_Lambda_P,r,1e-5,1000,L_P,L_P');
	fe						      = F*b;
	theta(s)				      = var(fe)-(ydelta'*Lambda_B*eta_h)/(NT-1);
	
	%Get the SE
	%[V_sym V_nosym]			  = MC_fixed_FD(ydelta,Chamard,Dsym,vec,id_movers(movers_sel),movers_sel,index_movers,index_movers_i,index_movers_2,index_movers_i_2,Ppath_1,Ppath_2,no_A2,list_problem,taxonomy,magical_number,TAXONOMY_I,TAXONOMY_J,MOVERS_i,MOVERS_j,FOCUS_I, FOCUS_J,Ppath_PROBLEM_i,Ppath_PROBLEM_j,taxonomyALT,TAXONOMY_ALT_I,TAXONOMY_ALT_J,MOVERS_ALT_i,MOVERS_ALT_j,FOCUS_ALT_I, FOCUS_ALT_J,Ppath_ALT_PROBLEM_i,Ppath_ALT_PROBLEM_j);
	[V_sym V_nosym]				  = MC_fixed_FD(ydelta,Chamard.Value,Dsym.Value,vec,id_movers(movers_sel),movers_sel,index_movers.Value,index_movers_i.Value,index_movers_2.Value,index_movers_i_2.Value,Ppath_1.Value,Ppath_2.Value,no_A2.Value,list_problem.Value,taxonomy.Value,magical_number.Value,TAXONOMY_I.Value,TAXONOMY_J.Value,MOVERS_i.Value,MOVERS_j.Value,FOCUS_I.Value, FOCUS_J.Value,Ppath_PROBLEM_i.Value,Ppath_PROBLEM_j.Value,taxonomyALT.Value,TAXONOMY_ALT_I.Value,TAXONOMY_ALT_J.Value,MOVERS_ALT_i.Value,MOVERS_ALT_j.Value,FOCUS_ALT_I.Value, FOCUS_ALT_J.Value,Ppath_ALT_PROBLEM_i.Value,Ppath_ALT_PROBLEM_j.Value);
	%[V_sym V_nosym]			  = MC_fixed_FD_v2(ydelta,Fdelta,L,F,tol,epsilon,type_algorithm,pfun_,Lambda_P,Lambda_B,Dsym.Value,vec,id_movers(movers_sel),movers_sel,index_movers.Value,index_movers_i.Value,index_movers_2.Value,index_movers_i_2.Value,Ppath_1.Value,Ppath_2.Value,no_A2.Value,list_problem.Value,taxonomy.Value,magical_number.Value,TAXONOMY_I.Value,TAXONOMY_J.Value,MOVERS_i.Value,MOVERS_j.Value,FOCUS_I.Value, FOCUS_J.Value,Ppath_PROBLEM_i.Value,Ppath_PROBLEM_j.Value,taxonomyALT.Value,TAXONOMY_ALT_I.Value,TAXONOMY_ALT_J.Value,MOVERS_ALT_i.Value,MOVERS_ALT_j.Value,FOCUS_ALT_I.Value, FOCUS_ALT_J.Value,Ppath_ALT_PROBLEM_i.Value,Ppath_ALT_PROBLEM_j.Value);
	%[V_sym_rand V_nosym_rand] 	  = paths_network_old(id_movers,firmid_delta,firmid_delta_f,ydelta,C);
	
	
	%Normalize
	V_nosym_SIM(s) 				  = sqrt(V_nosym/NT^2);
	V_sym_SIM(s) 				  = sqrt(V_sym/NT^2);
	%V_nosym_rand_SIM(s) 		  = sqrt(V_nosym_rand/NT^2);
	%V_sym_rand_SIM(s) 			  = sqrt(V_sym_rand/NT^2);
	
	%COVERAGE
	UB							  = theta(s)+1.96*V_sym_SIM(s);
	LB							  = theta(s)-1.96*V_sym_SIM(s);
	coverage_SYM(s)				  = (sigma2_true>=LB && sigma2_true<= UB);
	
	UB							  = theta(s)+1.96*V_nosym_SIM(s);
	LB							  = theta(s)-1.96*V_nosym_SIM(s);
	coverage_NOSYM(s)			  = (sigma2_true>=LB && sigma2_true<= UB);
	
	%UB							  = theta(s)+1.96*V_sym_rand_SIM(s);
	%LB							  = theta(s)-1.96*V_sym_rand_SIM(s);
	%coverage_rand_SYM(s)		  = (sigma2_true>=LB && sigma2_true<= UB);
	
	%UB							  = theta(s)+1.96*V_nosym_rand_SIM(s);
	%LB							  = theta(s)-1.96*V_nosym_rand_SIM(s);
	%coverage_rand_NOSYM(s)		  = (sigma2_true>=LB && sigma2_true<= UB);
	
	disp('done')
	
end
s=[filename 'MCCC_FD_NEW'];
save(s)

sel  							  = (imag(V_sym_SIM)>0) + (imag(V_nosym_SIM)>0);
sel								  = sel>0;
s=['Share of draws where I got imaginary SEs'];
disp(s);
mean(sel)


theta			= theta(~sel);
V_sym_SIM		= V_sym_SIM(~sel);
V_nosym_SIM		= V_nosym_SIM(~sel);
coverage_SYM	= coverage_SYM(~sel);
coverage_NOSYM	= coverage_NOSYM(~sel);

oracle_UB=theta+1.96*std(theta);
oracle_LB=theta-1.96*std(theta);
oracle_coverage=(sigma2_true>=oracle_LB).*(sigma2_true<=oracle_UB);
 

s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['True Value of Variance of Firm Effects: ' num2str(sigma2_true)];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['Expected Value of Leave one Out Variance Estimator: ' num2str(mean(theta))];
disp(s);
s=['Std of Estimated Variance of Firm Effects: ' num2str(std(theta))];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['Expected Value of standard error estimator (Symmetry): '    num2str(mean(V_sym_SIM))];
disp(s);
s=['Expected Value of standard error estimator (No-Symmetry): ' num2str(mean(V_nosym_SIM))];
disp(s);
%s=['Expected Value of Rand standard error estimator (Symmetry): '    num2str(mean(V_sym_rand_SIM))];
%disp(s);
%s=['Expected Value of Rand standard error estimator (No-Symmetry): ' num2str(mean(V_nosym_rand_SIM))];
%disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['Std Deviation of standard error estimator (Symmetry): '    num2str(std(V_sym_SIM))];
disp(s);
s=['Std Deviation of standard error estimator (No-Symmetry): ' num2str(std(V_nosym_SIM))];
disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['p1 p10 p25 p50 p75 p90 p99 of mathcal{V}%:']
disp(s)
quantile(V_sym_SIM,[0.01 0.10 0.25 0.50 0.75 0.90 0.99])
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
s=['Coverage Rate of oracle estimator 95%: ' num2str(mean(oracle_coverage))];
disp(s);
s=['Coverage Rate of Proposed Bounds for inference (SYMMETRIC) 95%: ' num2str(mean(coverage_SYM))];
disp(s);
s=['Coverage Rate of Proposed Bounds for inference (NON-SYMMETRIC) 95%: ' num2str(mean(coverage_NOSYM))];
disp(s);
%s=['Coverage Rate of Proposed Bounds for inference (Rand SYMMETRIC) 95%: ' num2str(mean(coverage_rand_SYM))];
%disp(s);
%s=['Coverage Rate of Proposed Bounds for inference (Rand NON-SYMMETRIC) 95%: ' num2str(mean(coverage_rand_NOSYM))];
%disp(s);
s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
disp(s)
stdTheta=std(theta);
theta=mean(theta);
sdV_sym_SIM=std(V_sym_SIM);
sdV_nosym_SIM=std(V_nosym_SIM);
V_sym_SIM=mean(V_sym_SIM);
V_nosym_SIM=mean(V_nosym_SIM);
oracle_coverage=mean(oracle_coverage);
coverage_SYM=mean(coverage_SYM);
coverage_NOSYM=mean(coverage_NOSYM);

end
end

