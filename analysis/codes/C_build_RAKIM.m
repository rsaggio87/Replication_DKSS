function vec = C_build_RAKIM(v,X_unrestrict,X_restrict,Pii_unrestrict,Pii_restrict);
%This computes C*v where v is a general Ndelta x 1 vector and the C matrix is the one from KSS and
%corresponds to the case where we want to test for equality of firm effects across two the PVR and the AKM model.

%Auxiliary
Ndelta							= size(X,1);
	
%Compute v_pred_restricted
X								= X_restrict;
xx								= (X'*X);
xy  							= X'*v;
[b flag]						= pcg(xx,xy,1e-10,1000);
v_pred_restricted				= X*b;

%Compute v_pred_unrestricted
X								= X_unrestrict;
xx								= (X'*X);
xy  							= X'*v;
[b flag]						= pcg(xx,xy,1e-10,1000);
v_pred_unrestricted				= X*b;

%Next line is essentially B*v
my_first_part					= (v_pred_unrestricted-v_pred_restricted);

%Now construct Lambda_B, Lambda_P
Lambda_B						= sparse((1:Ndelta)',(1:Ndelta), Pii_unrestrict-Pii_restrict,Ndelta,Ndelta);
Lambda_P						= sparse((1:Ndelta)',(1:Ndelta), Pii_unrestrict,Ndelta,Ndelta);
I_Lambda_P						= (speye(Ndelta,Ndelta)-Lambda_P);

%Now compute xi=M*(I-Lambda_P)^(-1)*Lambda_B*v
xy								= Lambda_B*v;
L_P								= ichol(I_Lambda_P,struct('type','ict','droptol',1e-2,'diagcomp',2));
[ydelta_xi, flag]				= pcg(I_Lambda_P,xy,1e-10,1000,L_P,L_P');
xy  							= X'*ydelta_xi;
[b flag]						= pcg(xx,xy,1e-10,1000);
pred							= X*b;
xi								= ydelta_xi-pred;

%Now compute Lambda_B*(I-Lambda_P)^(-1)*M*v;
res								= v-v_pred_unrestricted;
[res, flag]						= pcg(I_Lambda_P,res,1e-10,1000,L_P,L_P');
my_second_part					= 0.5*(Lambda_B*res+xi);


%Finish up
vec								= my_first_part-my_second_part;

end

