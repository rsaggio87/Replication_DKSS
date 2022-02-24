	function [coverage, coverage_AM]  = mc_read(filename,NSIM)
	
	%Load MC results 
	theta			= zeros(NSIM,1);
	theta_homo		= zeros(NSIM,1);
	theta_AKM		= zeros(NSIM,1);
	V_SIM			= zeros(NSIM,1);
	coverage		= zeros(NSIM,1);
	coverage_AM		= zeros(NSIM,1);
	coverage_oracle	= zeros(NSIM,1);
	COV_R1_sim		= zeros(2,2,NSIM);
	bar_Beta_1_simul= zeros(NSIM,1);
	theta_2_simul	= zeros(NSIM,1);
	Curvature 		= zeros(NSIM,1);
	taxnm			= zeros(NSIM,1);
	
	parfor s=1:NSIM
	name_save		= ['results/SIMULATIONS/result_sim' num2str(s)]
	load('name_save')
	theta(s)		= theta_s;
	theta_homo(s)	= theta_homo_s;
	theta_AKM(s)	= theta_AKM_s
	V_SIM(s)		= V_s;
	coverage(s)		= coverage_s;
	COV_R1_sim(:,:,s)=COV_R1_s;
	bar_Beta_1_simul(s)=bar_Beta_1_s;
	theta_2_simul(s)=theta_2_s;
	taxnm(s)		= taxnm_s;
	end
	


%Calculate oracle
	COV_R1_oracle					  = cov(bar_Beta_1_simul,theta_2_simul);
	correlation_oracle				  = corr(bar_Beta_1_simul,theta_2_simul);
	gamma_sq_oracle					  = ((lambda_1^2/NT^2)*(COV_R1_oracle(1,1)^2))/(COV_R1_oracle(2,2));
	oracle_coverage_fancy			  = zeros(NSIM,1);
	Curvature_true					  = zeros(NSIM,1);
	parfor s=1:NSIM
		[UB,LB,Curvature_true(s)] 	  = AM_CI(NT,lambda_1,COV_R1_oracle,bar_Beta_1_simul(s),theta_2_simul(s));
		oracle_coverage_fancy(s) 	  = (sigma2_true>=LB).*(sigma2_true<=UB);
		lunghezza_oracle(s)			  = UB-LB;
	
		[UB,LB,Curvature(s)]	  	  = AM_CI(NT,lambda_1,COV_R1_sim(:,:,s),bar_Beta_1_simul(s),theta_2_simul(s));
		coverage_AM(s)			      = (sigma2_true>=LB && sigma2_true<= UB);
		lunghezza_q1(s)				  = UB-LB;
	end

	save([filename '_RESULTS_MC'])

%Imaginary problems
		s=['Tabulate cases'];
		tabulate(taxnm)

		oracle_UB=theta+1.96*std(theta);
		oracle_LB=theta-1.96*std(theta);
		oracle_coverage=(sigma2_true>=oracle_LB).*(sigma2_true<=oracle_UB);

		s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
		disp(s)
		s=['True Value of Variance of Firm Effects: ' num2str(sigma2_true)];
		disp(s);
		s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
		disp(s)
		s=['Expected Value of KSS: ' num2str(mean(theta))];
		disp(s);
		s=['Std of KSS: ' num2str(std(theta))];
		disp(s);
		s=['Expected Value of Andrews: ' num2str(mean(theta_homo))];
		disp(s);
		s=['Std of Andrews: ' num2str(std(theta_homo))];
		disp(s);
		s=['Expected Value of AKM: ' num2str(mean(theta_AKM))];
		disp(s);
		s=['Std of AKM: ' num2str(std(theta_AKM))];
		disp(s);
		s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
		disp(s);
		s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
		disp(s)
		s=['Expected Value of standard error estimator: '    num2str(mean(V_SIM))];
		disp(s);
		s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
		disp(s)
		s=['Std Deviation of standard error estimator: '    num2str(std(V_SIM))];
		disp(s);
		s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
		disp(s)
		s=['p1 p10 p25 p50 p75 p90 p99 of mathcal{V}%:']
		disp(s)
		quantile(V_SIM,[0.01 0.10 0.25 0.50 0.75 0.90 0.99])
		s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
		disp(s)
		s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
		disp(s)
		s=['Coverage Rate of oracle estimator 95%: ' num2str(mean(oracle_coverage))];
		disp(s);
		s=['Coverage Rate of Proposed Bounds for inference 95%: ' num2str(mean(coverage))];
		disp(s);
		s=['Coverage Rate of Proposed Bounds for inference (AM) 95%: ' num2str(mean(coverage_AM))];
		disp(s);
		s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
		disp(s)
		s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
		disp(s)
		s=['-*-*-*-*-*-*-*-*-*-*-*-*'];
		disp(s)
		s=['Oracle SIGMA:'];
		disp(s);
		num2str(COV_R1_oracle)

		s=['Mean Estimated SIGMA:'];
		disp(s);
		num2str(mean(COV_R1_sim,3))

		s=['-*-*-*-*-*-*-*-*-*-*-*-*'];

end