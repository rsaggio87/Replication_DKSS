function  diff = simulated_pcg_levels(X,N,J,Jlag,network_to_look_at);


%Auxiliaries
	numIterations			= 10000;
	tol						= 1e-6;
	K						= size(X,2);
	
%Set up the pseudo parameters
	beta_sim			    = rand(K,1);	
	y						= X*beta_sim;

%Estimate
	xx					    = X'*X;
	xy						= X'*y;
	b						= pcg(xx,xy,tol,numIterations);
	
%Give me a reference group
if network_to_look_at == 2
	L						= X(:,N+J+1:N+J+Jlag)'*X(:,N+J+1:N+J+Jlag);
	L						= diag(L);
	index_firm				= max(find(L==max(L))); %maximum degree firm in corresponding network
	b_normalized			= beta_sim(N+J+index_firm);	
	beta_sim				= beta_sim-b_normalized;
	beta_sim				= beta_sim(N+J+1:N+J+Jlag);
	b						= b-b(N+J+index_firm);
	b						= b(N+J+1:N+J+Jlag);	
end

if network_to_look_at == 1
	L						= X(:,N+1:N+J)'*X(:,N+1:N+J);
	L						= diag(L);
	index_firm				= max(find(L==max(L))); %maximum degree firm in corresponding network
	b_normalized			= beta_sim(N+index_firm);	
	beta_sim				= beta_sim-b_normalized;
	beta_sim				= beta_sim(N+1:N+J);
	b						= b-b(N+index_firm);
	b						= b(N+1:N+J);	
end		

%Compare and contrast
	diff					= abs(b-beta_sim);

end	
