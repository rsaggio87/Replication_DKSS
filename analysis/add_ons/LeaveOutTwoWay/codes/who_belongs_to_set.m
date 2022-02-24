function in_set = who_belongs_to_set(index_movers,index_movers_i,prob_pairs);


%whether we want A_i1 or A_i2  depends on the inputs index_movers and index_movers_i

%whether we want A_i1 or A_ell1 depends on the first column of prob_pairs.

sel						= ismember(index_movers_i,prob_pairs(:,1)); 
index_movers_i			= index_movers_i(sel);
index_movers			= index_movers(sel);
c2			    		= (1:size(index_movers,1))';
c1		   				= (1:size(prob_pairs,1))';
[c1,c2] 				= meshgrid(c1,c2);
giga_lista				= cat(2,c1',c2');
big_INDEX				= reshape(giga_lista,[],2);
big_INDEX_prob_pairs	= big_INDEX(:,1);
big_INDEX_A				= big_INDEX(:,2);
clear big_INDEX
giga_lista				= [prob_pairs(big_INDEX_prob_pairs,1) prob_pairs(big_INDEX_prob_pairs,2) index_movers_i(big_INDEX_A) index_movers(big_INDEX_A)];
sel						= (giga_lista(:,1)==giga_lista(:,3));
giga_lista				= giga_lista(sel,:);
big_INDEX_prob_pairs	= big_INDEX_prob_pairs(sel);
big_INDEX_A				= big_INDEX_A(sel);
in_set					= giga_lista(:,2)==giga_lista(:,4);
in_set					= splitapply(@max,in_set,big_INDEX_prob_pairs);


end

    
    
    