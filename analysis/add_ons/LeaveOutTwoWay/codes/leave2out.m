function [n_of_bad_workers,bad_workers,artic_points] = leave2out(i,list)

%Leave out i
	id_movers 		= list(:,1);
	firmid			= list(:,2);
	id_movers_orig	= list(:,3);
	n_movers		= max(id_movers);
	J 				= max(firmid); 
    
%Construct Adj. Matrix of Bipartite Graph leaving i out
	sel = ~ismember(id_movers,i);
    B=sparse(list(sel,1),list(sel,2),1,n_movers,J);
    G = [ sparse( n_movers, n_movers ), B; B.', sparse(J, J)];
    
%Inspect the resulting Graph: find the workers that constitute an
%articulation point
    [artic_points, CC] = biconnected_components(G);
    bad_workers=artic_points(artic_points<=n_movers);
    
%Now get the right index for these workers
    sel=ismember(id_movers,bad_workers);
    bad_workers=id_movers_orig(sel);
    artic_points=unique(bad_workers);
    n_of_bad_workers=size(bad_workers,1);
    
    if n_of_bad_workers>0
    	bad_workers = mean(list(i,3));
    end	
    
    if n_of_bad_workers==0
      	bad_workers = 0;
    end 
    
end    