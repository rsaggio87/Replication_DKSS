function vec = C_build_cck(v,GROUP_delta,group_stack,firmid_stack,id_stack,firmid_delta_f,firmid_delta,Pii_unrestrict,Pii_restrict);
%This computes C*v where v is a general Ndelta x 1 vector and the C matrix is the one from KSS and
%corresponds to the case where we want to test for equality of firm effects across two groups.

Ndelta							= size(v,1);
J								= max([firmid_delta_f;firmid_delta]);

%Build X'X for restricted model
Fdelta1							= sparse((1:Ndelta)',firmid_delta_f,1,Ndelta,J);
Fdelta2							= sparse((1:Ndelta)',firmid_delta,1,Ndelta,J);
Fdelta							= Fdelta1-Fdelta2;
clear Fdelta1 Fdelta2	

%First thing is to compute B*v=P_Z*v
L								= (Fdelta'*Fdelta); %Laplacian matrix
evalc('pfun_= cmg_sdd(L)'); %preconditioner for Laplacian matrices.
xy  							= Fdelta'*v;
[b flag]						= pcg(L,xy,1e-10,1000,pfun_);
v_pred_restricted				= Fdelta*b;

%Now use a simple routine that projects a vector onto the Fdelta associated with a given group (either 0/1).
v_pred_0						= project_group(0,v(GROUP_delta==0),group_stack,firmid_stack,id_stack);
v_pred_1						= project_group(1,v(GROUP_delta==1),group_stack,firmid_stack,id_stack);
v_pred							= [v_pred_0;v_pred_1];
clear v_pred_0 v_pred_1

%Next line is B*v
my_first_part					= (v_pred-v_pred_restricted);

%Now construct Lambda_B, Lambda_P
Lambda_B						= sparse((1:Ndelta)',(1:Ndelta), Pii_unrestrict-Pii_restrict,Ndelta,Ndelta);
Lambda_P						= sparse((1:Ndelta)',(1:Ndelta), Pii_unrestrict,Ndelta,Ndelta);
I_Lambda_P						= (speye(Ndelta,Ndelta)-Lambda_P);


%Now compute xi=M*(I-Lambda_P)^(-1)*Lambda_B*v
xy								= Lambda_B*v;
L_P								= ichol(I_Lambda_P,struct('type','ict','droptol',1e-2,'diagcomp',2));
[ydelta_xi, flag]				= pcg(I_Lambda_P,xy,1e-10,1000,L_P,L_P');
pred_0							= project_group(0,ydelta_xi(GROUP_delta==0),group_stack,firmid_stack,id_stack);
pred_1							= project_group(1,ydelta_xi(GROUP_delta==1),group_stack,firmid_stack,id_stack);
pred							= [pred_0;pred_1];
clear pred_0 pred_1
xi								= ydelta_xi-pred;

%Now compute Lambda_B*(I-Lambda_P)^(-1)*M*v;
res								= v-v_pred;
[res, flag]						= pcg(I_Lambda_P,res,1e-10,1000,L_P,L_P');
my_second_part					= 0.5*(Lambda_B*res+xi);


%Finish up
vec								= my_first_part-my_second_part;




end

