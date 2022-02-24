function sel=largest_connected_bipartite(id,firmid);
%find movers
    J=max(firmid);
    gcs = [NaN; id(1:end-1)];
    gcs = id~=gcs;
    lagfirmid=[NaN; firmid(1:end-1)];
    lagfirmid(gcs==1)=NaN; %%first obs for each worker
    move=(firmid~=lagfirmid);
    move(gcs==1)=0;
    move=accumarray(id,move);
    move=(move>0);
    move=move(id);
    
%id and firmids associated with movers only
    id_movers=id(move);
    firmid_mover=firmid(move);

%need to normalize id_movers and keep a dictionary. 
    id_movers_orig=id_movers;
    [~,~,id_movers]=unique(id_movers);
    n_movers=max(id_movers);
    
%unique pairs    
    list=[id_movers firmid_mover];
    list=unique(list,'rows');
    
%Construct Adj. Matrix of Bipartite Graph
    B=sparse(list(:,1),list(:,2),1,n_movers,J);
    G = [ sparse( n_movers, n_movers ), B; B.', sparse(J, J)];
    
 %Find Largest Connected set   
 	[sindex, sz]= components(G); %get connected sets
	idx			= find(sz==max(sz)); %find largest set
	id_sel		= find(sindex==idx); %firms in connected set
    too			= id_sel>max(id);
    id_sel(too) = [];
    sel			= ismember(id,id_sel);
    
end

