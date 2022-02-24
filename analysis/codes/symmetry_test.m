function [R2_AKM R2_SEQ R2_PAIR] = symmetry_test(id,firmid,firmid_orig,controls,y,scale);	
R2_AKM=0;
R2_SEQ=0;
R2_PAIR=0;
%Data must be sort by id-job#. All ids must be normalized. T=2, balanced.

%Auxiliaries
	numIterations			= 10000;
	tol						= 1e-10;
	tolProb					= 0.5;
	N						= max(id);
	J						= max(firmid);

%Count Job #
	count					= ones(size(y,1),1);
	gcs 					= cell2mat(accumarray(id,count,[],@(x){cumsum(x)}));
	
%Origin-Destination Firm
	firmid_origin			= firmid(gcs==1);
	firmid_destination		= firmid(gcs==2);
	
	firmid_origin_true		= firmid_orig(gcs==1);
	firmid_destination_true	= firmid_orig(gcs==2);
	 		
%Wage Change
	y						= y(gcs==2)-y(gcs==1);
	y						= y-mean(y);
	disp('average of wage change')
	mean(y)
	disp('variance of wage changes')
	var(y)	
	
%Controls change
	controls				= controls(gcs==2,:)-controls(gcs==1,:);

%AKM matrix
	F1						= sparse((1:N)',firmid_origin',1,N,J);	
	F2						= sparse((1:N)',firmid_destination',1,N,J);
	F_AKM					= F2-F1;	
			
%SEQUENCE matrix
	list					= [firmid_origin,firmid_destination];
	sel						= firmid_origin<firmid_destination;
	list(sel,:)				= [firmid_destination(sel) firmid_origin(sel)];
	[~,~,match_PAIR]		= unique([list(:,1) list(:,2)],'stable','rows');
	uno						= ones(size(firmid_origin,1),1);
	uno(sel)				= -1;
	PAIRS					= max(match_PAIR)
	F_PAIRS					= sparse((1:N)',match_PAIR',uno',N,PAIRS);
	
%PAIR matrix
	[~,~,match_OD]			= unique([firmid_origin,firmid_destination],'rows');
	OD_cells				= max(match_OD)
	F_OD					= sparse((1:N)',match_OD',1,N,OD_cells);	
	
%Tell me number of uniquely populated cells
	peso_edge				= accumarray(match_OD,1);
	peso_edge				= peso_edge(match_OD);
	disp('share of unique (o,d) cells travelled by only one worker')
	mean(peso_edge==1)	
	
%Ready to run estimation for all the three models:
	ESS						= zeros(3,1);
	eta						= zeros(N,3);
	Pii_final				= zeros(N,3);
	Pii_final_R2			= zeros(N,3);
	fe						= zeros(N,3);
	
	for model=1:3
	
	
		if model == 1
			X 				= [F_AKM];
			chiamami		= 'AKM Model'	
		end
		
		if model == 2
			X 				= [F_PAIRS];
			chiamami		= 'Sequence Model'
		end
		
		if model == 3
			X 				= [F_OD];
			chiamami		= 'Pairs Model'
		end
	
		
%Estimate the corresponding model
		xy						= X'*y;
		xx						= X'*X;
		L						= ichol(xx,struct('type','ict','droptol',1e-2,'diagcomp',.1));
		tic
		b						= pcg(xx,xy,tol,numIterations,L,L');	
		toc
		xb						= X*b;
		fe(:,model)				= xb;
		eta(:,model)			= y-xb;
			
%Report Adj R2		
		ESS(model)				= var(xb);
		TSS					    = var(y);
		R2						= ESS(model)/TSS;
		adjR2					= 1-((1-R2)*(N-1))/(N-size(X,2)-1);
		s=['# of parameters: '  chiamami  '  ' num2str(size(X,2))];
		disp(s)
		s=['R2: '  chiamami  '  ' num2str(R2)];
		disp(s)
		s=['Adj R2: '  chiamami  '  ' num2str(adjR2)];
		disp(s)


%Getting stuff ready for JLA
		%# of draws
		disp('# of Simulated Projections for JLL:')
		scale
		
%set up the matrices	
		N						= size(X,1);
		xx_c					= parallel.pool.Constant(xx);
		X						= parallel.pool.Constant(X);
		L						= parallel.pool.Constant(L);
	
%preallocate
		Pii						= zeros(N,1);
		Pii_R2					= zeros(N,1);


%RUNNING JLA for AKM Model
if 1 == 1
tic
 	parfor i=1:scale
                ons 	= (rand(1,N) > tolProb);
                ons 	= ons - not(ons);
                ons 	= ons./sqrt(scale);
                  
                %Get me the row
                [Z, flag]= pcg(xx_c.Value,(ons*(X.Value))',tol,numIterations,L.Value,(L.Value)');
                
                %Collect
                Z		 = X.Value*Z;
                Z		 = Z.^2;
                Pii		 = Pii+Z;
                
                %Mean adjusted Pii
				ons		 = ons-mean(ons);
				[Z, flag]= pcg(xx_c.Value,(ons*(X.Value))',tol,numIterations,L.Value,(L.Value)')
				Z		 = X.Value*Z;
                Z		 = Z.^2;
                Pii_R2	 = Pii_R2+Z;

	end
		Pii_final(:,model)	 = Pii;
		Pii_final_R2(:,model)= Pii_R2;
		s=['Min(Pii): '  chiamami  '  ' num2str(min(Pii))];
		disp(s)
		s=['Mean(Pii): '  chiamami  '  ' num2str(mean(Pii))];
		disp(s)
		s=['Max(Pii): '  chiamami  '  ' num2str(max(Pii))];
		disp(s)
		mean(Pii_final(:,model))
		max(Pii_final(:,model))
		min(Pii_final(:,model))
	end	
toc	
end	

%Get the sigma_i of the restricted model
	invP				    = 1-Pii_final(:,1);
	eta_h				    = eta(:,1)./invP;
	sigma_i					= y.*eta_h;
	
	s=['Min(sigma_i): '    '  ' num2str(min(sigma_i))];
	disp(s)
	s=['Mean(sigma_i): '    '  ' num2str(mean(sigma_i))];
	disp(s)
	s=['Max(sigma_i): '    '  ' num2str(max(sigma_i))];
	disp(s)

%Compute the KSS R2 now imposing the null
	R2_KSS					= zeros(3,1);
	for model=1:3
	R2_KSS(model)			= ESS(model)-mean(Pii_final_R2(:,model).*sigma_i);
	R2_KSS(model)			= R2_KSS(model)/var(y);	
	end

%Done
	R2_AKM 					= R2_KSS(1)
	R2_SEQ 					= R2_KSS(2)
	R2_PAIR					= R2_KSS(3)
	
%Output csv
	out						= [y firmid_origin firmid_destination firmid_origin_true firmid_destination_true fe(:,1) fe(:,2) fe(:,3)];
	s						= ['/scratch/public/leave_out/symmetry.csv'];
	dlmwrite(s, out, 'delimiter', '\t', 'precision', 16); 	
	
%Save
	save('aux')					
end
