function [R2_interacted R2_interacted_KSS]	= interacted_model_levels(model,y,Pii,sigma_i,id,firmid,lagfirmid,controls,cara);

numIterations			= 10000;
tol						= 1e-10;

%create the new firmids
[~,~,firmid]   			= unique([firmid cara],'rows','stable');
[~,~,lagfirmid]   		= unique([lagfirmid cara],'rows','stable');

%create the matrices now
NT						= size(id,1);
J						= max(firmid);
Jlag					= max(lagfirmid);
N						= max(id);

F						= sparse((1:NT)',firmid',1,NT,J);
Flag					= sparse((1:NT)',lagfirmid',1,NT,Jlag);	
D						= sparse((1:NT)',id',1,NT,N);

if model == 3
X						= [D F Flag controls];
end

if model == 1
X						= [D Flag controls];
end

if model == 2
X						= [D F controls];
end


%Estimate the model
xy						= X'*y;
xx						= X'*X;
b						= pcg(xx,xy,tol,numIterations);	
xb						= X*b;	

%R2
ESS					    = var(xb);
TSS					    = var(y);
R2					    = ESS/TSS
R2_interacted			= 1-((1-R2)*(NT-1))/(NT-size(X,2)-1)

%Correct the R2 now using the model under the null
R2_interacted_KSS		= ESS - mean(Pii.*sigma_i);
R2_interacted_KSS		= R2_interacted_KSS/TSS;

end	
