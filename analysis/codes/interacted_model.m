function [R2_interacted R2_interacted_KSS]	= interacted_model(y,Pii,sigma_i,id,firmid,controls,cara);

%key auxiliary
NT						= length(y);
count					= ones(NT,1);
gcs 					= cell2mat(accumarray(id,count,[],@(x){cumsum(x)}));

%create the new firmids
[~,~,firmid]   			= unique([firmid cara],'rows','stable');

%now build the estimation sample for the interacted model
sel						= (gcs == 1 | gcs ==2);
firmid_orig1			= firmid(sel); %T=2.
[~,~,firmid_use]		= unique(firmid_orig1);
gcs_use					= gcs(sel);
firmid_lag_1			= firmid_use(gcs_use==1);
firmid_con_1			= firmid_use(gcs_use==2);
F1						= sparse((1:max(id))',firmid_lag_1',1,max(id),max(firmid_use));
F2						= sparse((1:max(id))',firmid_con_1',1,max(id),max(firmid_use));
FDELTA1					= F2-F1;		
sel						= (gcs == 2 | gcs ==3);
firmid_orig2			= firmid(sel);
[~,~,firmid_use]		= unique(firmid_orig2);
gcs_use					= gcs(sel);
firmid_lag_2			= firmid_use(gcs_use==2);
firmid_con_2			= firmid_use(gcs_use==3);
F1						= sparse((1:max(id))',firmid_lag_2',1,max(id),max(firmid_use));
F2						= sparse((1:max(id))',firmid_con_2',1,max(id),max(firmid_use));
FDELTA2					= F2-F1; 

%Controls
controls_1				= controls(gcs==2,:);
controls_2				= controls(gcs==3,:);
controls_DELTA  		= controls_2-controls_1; 

%Estimate the augmented PVR model, find the PI R2
ydelta					= y(gcs==3)-y(gcs==2);
X						= [FDELTA1 FDELTA2 controls_DELTA ones(size(FDELTA1,1),1)];
xx						= X'*X;
xy						= X'*ydelta;
b						= pcg(xx,xy,1e-10,10000);	
xb						= X*b;
ESS					    = var(xb);
TSS					    = var(ydelta);
R2_interacted			= ESS/TSS;
NT						= size(y,1);
dof						= size(y,1)-size(X,2)-1;
R2						= 1-(1-R2_interacted)*((NT-1)/dof);

%Correct the R2 now using the model under the null
R2_interacted_KSS		= ESS - mean(Pii.*sigma_i);
R2_interacted_KSS		= R2_interacted_KSS/TSS;



end	
