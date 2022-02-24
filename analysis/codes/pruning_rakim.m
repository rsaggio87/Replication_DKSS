function [y,id,firmid,gcs,controls,ydelta,X,FDELTA1,FDELTA2,firmid_orig1,firmid_orig2,firmid_lag_1,firmid_con_1,firmid_lag_2,firmid_con_2,firmid_DICTIO1,firmid_DICTIO2,cara]=pruning_rakim(y,id,firmid,controls,cara)

bad_obs					= 1;
firmid_orig				= firmid;
id_orig					= id;

while bad_obs>0			
			%Calculate the key auxiliary
			NT				= length(y);
			count			= ones(NT,1);
			gcs 			= cell2mat(accumarray(id,count,[],@(x){cumsum(x)}));
			tabulate(gcs)
			obs_count		= (1:NT)';
			
			if 0 == 1 
			%Leave one out connected set
			[~,~,~,~,~,obs_sel] = pruning_unbal_v3(y,firmid,id,id,firmid,obs_count);
			id_sel			= id(obs_sel);
			sel				= ismember(id,id_sel);
			id				= id(sel)	;
			firmid			= firmid(sel);
			y				= y(sel);
			controls		= controls(sel,:);
			[~,~,id]		= unique(id);
			[~,~,firmid]	= unique(firmid);
			
			s=['After we found the leave one out connected set in the AKM model'];
			disp(s)
			s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
			s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
			disp(s)
			disp(s)
			s=['# of workers: ' int2str(max(id))];
			disp(s);
			s=['# of firms: ' int2str(max(firmid))];
			disp(s);
			s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
			s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
			disp(s);
			
			%Redefine auxiliary
			NT				= length(y);
			count			= ones(NT,1);
			gcs 			= cell2mat(accumarray(id,count,[],@(x){cumsum(x)}));
			tabulate(gcs)
			end
									
			%Define Fdelta for Network 1
			sel				= (gcs == 1 | gcs ==2);
			firmid_orig1	= firmid(sel); %T=2.
			[~,~,firmid_use]= unique(firmid_orig1);
			firmid_DICTIO1	= accumarray(firmid_use,firmid_orig(sel),[],@(x)mean(x));
			gcs_use			= gcs(sel);
			firmid_lag_1	= firmid_use(gcs_use==1);
			firmid_con_1	= firmid_use(gcs_use==2);
			firmid_orig_lag_1= firmid_orig1(gcs_use==1);
			firmid_orig_con_1= firmid_orig1(gcs_use==2);
			
			F1= sparse((1:max(id))',firmid_lag_1',1,max(id),max(firmid_use));
			F2= sparse((1:max(id))',firmid_con_1',1,max(id),max(firmid_use));
			FDELTA1			= F2-F1; 
			
			%Define Fdelta for Network 2			
			sel				= (gcs == 2 | gcs ==3);
			firmid_orig2	= firmid(sel);
			[~,~,firmid_use]= unique(firmid_orig2);
			firmid_DICTIO2	= accumarray(firmid_use,firmid_orig(sel),[],@(x)mean(x));
			gcs_use			= gcs(sel);
			firmid_lag_2	= firmid_use(gcs_use==2);
			firmid_con_2	= firmid_use(gcs_use==3);
			firmid_orig_lag_2= firmid_orig2(gcs_use==2);
			firmid_orig_con_2= firmid_orig2(gcs_use==3);
			F1				= sparse((1:max(id))',firmid_lag_2',1,max(id),max(firmid_use));
			F2				= sparse((1:max(id))',firmid_con_2',1,max(id),max(firmid_use));
			FDELTA2			= F2-F1; 
			clear F1 F2 gcs_use sel 
			
			%Controls
			controls_1		= controls(gcs==2,:);
			controls_2		= controls(gcs==3,:);
			controls_DELTA  = controls_2-controls_1; 
			
			%Estimate the PVR model
			ydelta			= y(gcs==3)-y(gcs==2);
			X				= [FDELTA1 FDELTA2 controls_DELTA ones(size(FDELTA1,1),1)];
			xx				= X'*X;
			xy				= X'*ydelta;
			b				= pcg(xx,xy,1e-10,10000);	
			xb				= X*b;
			
			%Find cases where we have a perfect prediction
			dist			= abs(ydelta-xb);
			good_ids		= dist>1e-4;
			good_ids		= find(good_ids);
			sel				= ismember(id,good_ids);
			bad_obs			= ~ismember(id,good_ids);
			id_orig			= id_orig(sel);
			firmid_orig		= firmid_orig(sel);
			id				= id(sel);
			firmid			= firmid(sel);
			y				= y(sel);
			[~,~,id]		= unique(id);
			[~,~,firmid]	= unique(firmid);
			controls		= controls(sel,:);
			cara			= cara(sel,:);
			
			%Drop these dudes and start over
			bad_obs			= sum(bad_obs)
			
			s=['After pruning obs with PVR Pii ==1'];
			disp(s)
			s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
			s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
			disp(s)
			disp(s)
			s=['# of workers: ' int2str(max(id))];
			disp(s);
			s=['# of firms: ' int2str(max(firmid))];
			disp(s);
			s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
			s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
			disp(s);
			save('pruned_data')		
end	
end	
