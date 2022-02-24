function [B_u B_l] = quality_JLL(Lambda_P,Lambda_B_fe,y,xx,X,Lchol,scale);
	
	%First step is to calculate b
	NT			 = size(y,1);
	xy			 = X'*y;
	[beta, flag] = pcg(xx,xy,1e-10,1000,Lchol);
	eta			 = y-X*beta;
	
	%Now Leave out resids.
	NT 			 = size(y,1);
	I_Lambda_P   = (speye(NT,NT)-Lambda_P);
	eta_h	     = I_Lambda_P\eta; %Leave one out residual
	
	%Now Bii,Pii
	Bii			= diag(Lambda_B_fe);
	Pii			= diag(Lambda_P);
	
	disp('Average Leverage')
	mean(Pii)
	
	disp('Std. Leverage')
	std(Pii)
	
	disp('Maximum Leverage')
	max(Pii)
	
	%Now construct the bound for the approximation
	sigma_i 	= y.*eta_h;
	den			= 1./(1-Pii);
	B_l			= (1/NT)*(2/scale)*sum(Bii.*sigma_i.*(2*Pii.^3).*den)
	B_u			= (1/NT)*(2/scale)*sum(Bii.*sigma_i.*(Pii.^2).*den)
	
end 
	