function  diff = simulated_pcg(X,FDELTA1,FDELTA2,network_to_look_at);


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
	L						= FDELTA2'*FDELTA2;
	L						= diag(L);
	index_firm				= find(L==max(L)); %maximum degree firm in contemporaneous network
	index_firm				= size(FDELTA1,2)+index_firm;
	b_normalized			= beta_sim(index_firm);
	beta_sim				= beta_sim-b_normalized;
	beta_sim				= beta_sim(size(FDELTA1,2)+1:size(FDELTA1,2)+1+size(FDELTA2,2)-1,1);
	b						= b-b(index_firm);
	b						= b(size(FDELTA1,2)+1:size(FDELTA1,2)+1+size(FDELTA2,2)-1,1);	
end

if network_to_look_at == 1
	L						= FDELTA1'*FDELTA1;
	L						= diag(L);
	index_firm				= find(L==max(L)); %maximum degree firm in contemporaneous network
	index_firm				= index_firm;
	b_normalized			= beta_sim(index_firm);	
	beta_sim				= beta_sim-b_normalized;
	beta_sim				= beta_sim((1:size(FDELTA1,2))');
	b						= b-b(index_firm);
	b						= b((1:size(FDELTA1,2))');	
end		

%Compare and contrast
	diff					= abs(b-beta_sim);

end	
