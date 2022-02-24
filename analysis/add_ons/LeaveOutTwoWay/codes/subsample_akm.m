function [yhat P] = subsample_akm(i,list,list_to_use);
				
				%what i have to predict
				j				= list(i,1);
				jprime			= list(i,2);
				
				%Build the predictor
				sel				= ismember(list(:,3), list_to_use);
				actual_path	    = [list(sel,1) list(sel,2)];
				NPATHS			= size(actual_path,1);
				[aa,bb,cc]		= unique(actual_path);
				cc				= reshape(cc,NPATHS,2);
				Nlist			= size(cc,1);
				Nodes			= max(max(cc));
				F				= sparse((1:Nlist)',cc(:,1),1,Nlist,Nodes) + sparse((1:Nlist)',cc(:,2),-1,Nlist,Nodes) ;
				outcomes		= list(sel,4);
				psi				= pinv(full(F'*F))*F'*outcomes;
				indexJ			= find(aa==j);
				indexJprime		= find(aa==jprime);
				yhat	 		= psi(indexJ) - psi(indexJprime);
				
				
				%Give me back the P associated with the mini regression
				aux 			= sparse(1,indexJ,1,1,Nodes) + sparse(1,indexJprime,-1,1,Nodes) ;
				P				= aux*(pinv(full(F'*F))*F');
			
end


    