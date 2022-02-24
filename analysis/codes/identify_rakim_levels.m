function  identified = identify_rakim_levels(id,firmid,lagfirmid,controls,network_to_look_at);

	%Set-up the matrices
	NT			= size(id,1);
	J			= max(firmid);
	Jlag		= max(lagfirmid);
	N			= max(id);

	F			= sparse((1:NT)',firmid',1,NT,J);
	Flag		= sparse((1:NT)',lagfirmid',1,NT,Jlag);	
	D			= sparse((1:NT)',id',1,NT,N);
	X			= [D F Flag controls];
	K			= size(X,2);

	%Run the simulated pcg
	diff		= simulated_pcg_levels(X,N,J,Jlag,network_to_look_at);
		
	%From this first run get an upper bound by inspecting the kink
	qpp						= (1:1:100);
	quantilI				= quantile(diff,qpp/100); 
	disp('Quantiles iteration 1')
	aiuto					= ([qpp;quantilI])
	quantilI				= quantilI';
	shifted					= [quantilI(1); quantilI(1:end-1)];
	shifted					= quantilI-shifted;
	disp('Quantile cut-off')
	kink_value				= 0.01;
	identified				= diff<=kink_value;
	disp('# of Identified effects in iteration 1')
	sum(identified)

	%Now keep iterating 
	N_identified_old 		= sum(identified);
	N_identified_new		= 0;
	convergenza				= abs(N_identified_old-N_identified_new);
	identified_old			= identified;
	ss						= 2;
	
	while convergenza > 0	
		%new draw
		diff				= simulated_pcg_levels(X,N,J,Jlag,network_to_look_at);
	
		%evaluate who is identified now
		identified			= (identified_old == 1 & diff <= kink_value);
		N_identified_new	= sum(identified);
		convergenza			= abs(N_identified_old-N_identified_new);
		disp(['# of Identified effects in iteration  ' num2str(ss)])
		N_identified_new
		
		%Update
		N_identified_old 	= N_identified_new;
		ss					= ss + 1;
		identified_old		= identified;
	end

end	
