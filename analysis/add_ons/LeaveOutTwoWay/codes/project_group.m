function prrred					= project_group(group,v,group_stack,firmid_stack,id_stack);

sel								= group_stack==group;

firmid							= firmid_stack(sel);
[~,~,firmid]					= unique(firmid);
J								= max(firmid);
id								= id_stack(sel);
[~,~,id]						= unique(id);

gcs 							= [NaN; id(1:end-1)];
gcs 							= id~=gcs;
firmid_delta					= [NaN; firmid(1:end-1)];
firmid_delta(gcs==1)			= NaN;
sel								= gcs~=1;
firmid_delta					= firmid_delta(sel);
firmid_delta_f					= firmid;
firmid_delta_f					= firmid_delta_f(sel);
Ndelta							= size(firmid_delta_f,1);


Fdelta1							= sparse((1:Ndelta)',firmid_delta_f,1,Ndelta,J);
Fdelta2							= sparse((1:Ndelta)',firmid_delta,1,Ndelta,J);
Fdelta							= Fdelta1-Fdelta2;

L								= (Fdelta'*Fdelta); %Laplacian matrix
evalc('pfun_= cmg_sdd(L)'); %preconditioner for Laplacian matrices.
xy  							= Fdelta'*v;
[b flag]						= pcg(L,xy,1e-10,1000,pfun_);
prrred							= Fdelta*b;

end