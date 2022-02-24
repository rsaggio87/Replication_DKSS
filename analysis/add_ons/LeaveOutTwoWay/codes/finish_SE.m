function [V, share_bp, COV_R1]=finish_SE(ydelta,taxonomy,A_i2_in_Aj1,A_i1_in_Aj2,A_i2_in_Aj2,vec2,pfun_,prob_pairs,j_in_Ai1,j_in_Ai2,i_in_Aj1,i_in_Aj2,index_movers_i,index_movers,index_movers_i_2,index_movers_2,no_A2_ii,no_A2_jj,Q_share,y_ii,y_jj,sigma_1_i,sigma_1_j,sigma_2_i,sigma_2_j,yhat_1_ii,yhat_2_ii,yhat_1_jj,yhat_2_jj,vec,sigma_2,Ndelta,Lambda_P,sigma_1,Fdelta,L,Lambda_B,do_weak,Lambda_B2,F,x1bar,lambda_1,D,magical_number,magical_number_q,Dq,sigma_i,GRANDE_MEDIA);

sigma_3_i				= y_ii.*(y_ii-yhat_2_ii);
sigma_3_j				= y_jj.*(y_jj-yhat_2_jj);

sigma_4_i				= (y_ii-GRANDE_MEDIA).^2;
sigma_4_j				= (y_jj-GRANDE_MEDIA).^2;

%Ready to fill the product of sigmas for the problematic pairs now,
[sigma_product bad_pairs tag] = compute_sigma_product(taxonomy,no_A2_jj,no_A2_ii,j_in_Ai1,j_in_Ai2,i_in_Aj1,i_in_Aj2,A_i2_in_Aj1,A_i1_in_Aj2,A_i2_in_Aj2,sigma_1_i,sigma_2_i,sigma_3_i,sigma_4_i,sigma_1_j,sigma_2_j,sigma_3_j,sigma_4_j,magical_number);
sigma_product_q			 = compute_sigma_product(taxonomy,no_A2_jj,no_A2_ii,j_in_Ai1,j_in_Ai2,i_in_Aj1,i_in_Aj2,A_i2_in_Aj1,A_i1_in_Aj2,A_i2_in_Aj2,sigma_1_i,sigma_2_i,sigma_3_i,sigma_4_i,sigma_1_j,sigma_2_j,sigma_3_j,sigma_4_j,magical_number_q);

%Now building the four pieces that compose the sampling variance estimator of the firm effects.
first_piece 			= 4*sum((vec.^2).*sigma_2);

%Use Hutchison trick and this (https://en.wikipedia.org/wiki/Hadamard_product_(matrices)) to get the second piece
I_Lambda_P				= (speye(Ndelta,Ndelta)-Lambda_P);
L_P						= ichol(I_Lambda_P,struct('type','ict','droptol',1e-2,'diagcomp',2));
VCM 					= parallel.pool.Constant(diag(sigma_1));
Fdelta  				= parallel.pool.Constant(Fdelta);
L  						= parallel.pool.Constant(L);
Lambda_B				= parallel.pool.Constant(Lambda_B);

if do_weak == 1
 	Lambda_B2			= parallel.pool.Constant(Lambda_B2);
 	Q_use				= size(x1bar,2);
end

I_Lambda_P				= parallel.pool.Constant(I_Lambda_P);
L_P  					= parallel.pool.Constant(L_P);
F 						= parallel.pool.Constant(F);   
%NSIM					= floor(20*0.05^(-2)*log(2/0.05)); %w.p 95% the absolute difference b/w the trace estimator and the actual trace is going to be less than 5% - see table 1 of http://www.cs.tau.ac.il/~stoledo/Bib/Pubs/trace3.pdf 
NSIM					= 1000;
aux_SIM					= zeros(NSIM,1);
aux_SIM_2				= zeros(NSIM,1);
%rng(3815)
%XSIMUL					= randn(Ndelta,NSIM);

    parfor s=1:NSIM    
        xsimul			= randn(Ndelta,1);
        %xsimul			= XSIMUL(:,s);       
	 	left			= C_build_FD(xsimul,Fdelta.Value,L.Value,pfun_,Lambda_B.Value,I_Lambda_P.Value,L_P.Value,F.Value);
	 	
	 	if do_weak == 1
	 		left_2		= C_build_FD(xsimul,Fdelta.Value,L.Value,pfun_,Lambda_B2.Value,I_Lambda_P.Value,L_P.Value,F.Value);
	 		for qqq=1:Q_use
	 		norm_fact	= (lambda_1(qqq)*x1bar(qqq))*sum(x1bar(qqq).*xsimul);
	 		left_2 		= left_2 - norm_fact;
	 		end
	 	end
	 	
		xsimul			= VCM.Value*xsimul;
		right			= C_build_FD(xsimul,Fdelta.Value,L.Value,pfun_,Lambda_B.Value,I_Lambda_P.Value,L_P.Value,F.Value);
		
		if do_weak == 1
	 		right_2		= C_build_FD(xsimul,Fdelta.Value,L.Value,pfun_,Lambda_B2.Value,I_Lambda_P.Value,L_P.Value,F.Value);
	 		for qqq=1:Q_use
	 		norm_fact	= (lambda_1(qqq)*x1bar(qqq))*sum(x1bar(qqq).*xsimul);
	 		right_2 	= right_2 - norm_fact;
	 		end
	 	end
		
		aux_SIM(s)		= left'*VCM.Value*right;
		
		if do_weak == 1
			aux_SIM_2(s)= left_2'*VCM.Value*right_2;
		end
    end  

second_piece			= 2*mean(aux_SIM);

%The last two pieces
Dsym					= 2*D + 2*D';
third_piece				= 2*sigma_1'*Dsym*sigma_1;
fourth_piece			= 2*sum((magical_number).*(sigma_product-(sigma_1_i.*sigma_1_j))); %multiply by 4 instead of 2 because I only considered upper triangular part of the that identifies problematic pairs

%Finished
V						= first_piece - second_piece - third_piece - fourth_piece;
taxonomy_V				= 0;
taxonomy_SIGMA_1		= 0;
taxonomy_SIGMA_2		= 0;

if V < 0						
	V					= first_piece;
	taxonomy_V			= 1;
end

if V < 0						
	V 					= 4*sum((vec.^2).*(ydelta-GRANDE_MEDIA).^2);
	taxonomy_V			= 2;
end



%VCM for weak id
if do_weak == 1

	%SIGMA(2,2)
	first_piece 								= 4*sum((vec2.^2).*sigma_2);
	second_piece								= 2*mean(aux_SIM_2);
	Dsym										= 2*Dq + 2*Dq';
	third_piece									= 2*sigma_1'*Dsym*sigma_1;
	fourth_piece								= 2*sum((magical_number_q).*(sigma_product_q-(sigma_1_i.*sigma_1_j)));
	COV_R1(Q_use+1,Q_use+1)						= first_piece - second_piece - third_piece - fourth_piece;
	
	if COV_R1(Q_use+1,Q_use+1)<0
				COV_R1(Q_use+1,Q_use+1) 		= first_piece;
				taxonomy_SIGMA_2				= 1;			
	end
	
	if COV_R1(Q_use+1,Q_use+1)<0
				COV_R1(Q_use+1,Q_use+1) 		= 4*sum((vec2.^2).*(ydelta-GRANDE_MEDIA).^2);;
				taxonomy_SIGMA_2				= 2;			
	end
	
	%SIGMA(1,1)
	VCM					    					= sparse((1:size(sigma_i,1))',1:size(sigma_i,1),sigma_i,size(sigma_i,1),size(sigma_i,1));
	COV_R1(1:Q_use,1:Q_use) 					= (x1bar'*VCM*x1bar);
	
	if COV_R1(1:Q_use,1:Q_use)	< 0
				VCM								= sparse((1:size(sigma_i,1))',1:size(sigma_i,1),(ydelta-GRANDE_MEDIA).^2,size(sigma_i,1),size(sigma_i,1));
				COV_R1(1:Q_use,1:Q_use)			= (x1bar'*VCM*x1bar);
				taxonomy_SIGMA_1				= 1;
	end
	
	%SIGMA(1,2)
	COV_R1(1:Q_use,Q_use+1)						= 2*x1bar'*(sigma_2.*vec2);
	COV_R1(Q_use+1,1:Q_use)		    			= COV_R1(1:Q_use,Q_use+1)';

end


%Taxonomy of cases
if taxonomy_V == 0 & taxonomy_SIGMA_1 == 0 & taxonomy_SIGMA_2 == 0 
share_bp				= 1;
end

if taxonomy_V == 0 & taxonomy_SIGMA_1 == 0 & taxonomy_SIGMA_2 == 1 
share_bp				= 2;
end

if taxonomy_V == 0 & taxonomy_SIGMA_1 == 1 & taxonomy_SIGMA_2 == 0
share_bp				= 3;
end

if taxonomy_V == 0 & taxonomy_SIGMA_1 == 1 & taxonomy_SIGMA_2 == 1
share_bp				= 4;
end

if taxonomy_V == 1 & taxonomy_SIGMA_1 == 0 & taxonomy_SIGMA_2 == 0 
share_bp				= 5;
end

if taxonomy_V == 1 & taxonomy_SIGMA_1 == 0 & taxonomy_SIGMA_2 == 1 
share_bp				= 6;
end

if taxonomy_V == 1 & taxonomy_SIGMA_1 == 1 & taxonomy_SIGMA_2 == 0
share_bp				= 7;
end

if taxonomy_V == 1 & taxonomy_SIGMA_1 == 1 & taxonomy_SIGMA_2 == 1
share_bp				= 8;
end

if taxonomy_V == 2 & taxonomy_SIGMA_1 == 0 & taxonomy_SIGMA_2 == 0 
share_bp				= 9;
end

if taxonomy_V == 2 & taxonomy_SIGMA_1 == 0 & taxonomy_SIGMA_2 == 1 
share_bp				= 10;
end

if taxonomy_V == 2 & taxonomy_SIGMA_1 == 1 & taxonomy_SIGMA_2 == 0
share_bp				= 11;
end

if taxonomy_V == 2 & taxonomy_SIGMA_1 == 1 & taxonomy_SIGMA_2 == 1
share_bp				= 12;
end

if taxonomy_V == 0 & taxonomy_SIGMA_2 == 2 & taxonomy_SIGMA_1 == 0 
share_bp				= 13;
end

if taxonomy_V == 0 & taxonomy_SIGMA_2 == 2 & taxonomy_SIGMA_1 == 1 
share_bp				= 14;
end

if taxonomy_V == 1 & taxonomy_SIGMA_2 == 2 & taxonomy_SIGMA_1 == 0
share_bp				= 15;
end

if taxonomy_V == 1 & taxonomy_SIGMA_2 == 2 & taxonomy_SIGMA_1 == 1
share_bp				= 16;
end

if taxonomy_V == 2 & taxonomy_SIGMA_2 == 2 & taxonomy_SIGMA_1 == 0
share_bp				= 17;
end

if taxonomy_V == 2 & taxonomy_SIGMA_2 == 2 & taxonomy_SIGMA_1 == 1
share_bp				= 18;
end


end

    
    
    