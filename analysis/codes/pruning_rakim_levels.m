function [y,id,firmid,lagfirmid,controls,id_orig,firmid_orig,lagfirmid_orig,cara]=pruning_rakim_levels(y,id,firmid,lagfirmid,controls,cara,tolleranza);	

if nargin < 7 
tolleranza = 1e-3
end 

%crucial: all ids must already be normalized.

bad_obs					= 1;
firmid_orig				= firmid;
id_orig					= id;
lagfirmid_orig			= lagfirmid;

%Normalize ids'
[~,~,id]				= unique(id);
[~,~,firmid]			= unique(firmid);
[~,~,lagfirmid]			= unique(lagfirmid);
NT=size(id,1);

%Summarize starting point 
	s=['Summarize Starting Point:'];
	disp(s)
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	disp(s)
	s=['Mean of y: ' num2str(mean(y))];
	disp(s);
	s=['Variance of y: ' num2str(var(y))];
	disp(s);
	s=['# of Jobs x Person observations: ' int2str(NT)];
	disp(s);
	s=['# of workers: ' int2str(max(id))];
	disp(s);
	s=['# of firms: ' int2str(max(firmid))];
	disp(s);
	s=['# of lagged firms: ' int2str(max(lagfirmid))];
	disp(s);
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s)
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s);
	s=['STARTING PRUNING.........'];
	disp(s)
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
	disp(s);
	

while bad_obs>0			
			
			%Sparsify
			NT			= size(y,1);
			J			= max(firmid);
			Jlag		= max(lagfirmid);
			N			= max(id);
			
			F			= sparse((1:NT)',firmid',1,NT,J);
			Flag		= sparse((1:NT)',lagfirmid',1,NT,Jlag);	
			D			= sparse((1:NT)',id',1,NT,N);
			X			= [D F Flag controls];
			K			= size(X,2);
			
			xx			= X'*X;
			xy			= X'*y;
			tic
			b			= pcg(xx,xy,1e-10,1e8);	
			toc
			xb			= X*b;
			
			%Find cases where we have a perfect prediction
			dist		= abs(y-xb);
			sel			= dist>tolleranza;
			
			id_orig		= id_orig(sel);
			firmid_orig	= firmid_orig(sel);
			lagfirmid_orig=lagfirmid_orig(sel);
			id			= id(sel);
			firmid		= firmid(sel);
			lagfirmid	= lagfirmid(sel);
			y			= y(sel);
			controls	= controls(sel,:);
			cara		= cara(sel,:);
			
			[~,~,id]	= unique(id);
			[~,~,firmid]= unique(firmid);
			[~,~,lagfirmid]= unique(lagfirmid);
			
			%Drop these dudes and start over
			bad_obs		= sum(~sel)
			
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
			s=['# of lagged firms: ' int2str(max(lagfirmid))];
			disp(s);
			s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
			s=['-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*'];
			disp(s);	
end	
end	