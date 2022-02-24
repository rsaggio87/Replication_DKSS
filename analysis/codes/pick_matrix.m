function [DUAL_F DUAL_Flag]				= pick_matrix(dual_list,settore,N,J,Jlag,K);			
				
				sel						= (dual_list(:,6)==settore); %cut the list to only firms that pertain to the sector that you care about
				dual_list_usami			= dual_list(sel,:);
				N_dual					= size(dual_list_usami,1);
				
				sel						= dual_list_usami(:,2)';
				DUAL_F					= [sparse(N_dual,N) 	 sparse((1:N_dual)',sel,1,N_dual,J)   	  sparse(N_dual,K-N-J)];	
				DUAL_F					= repelem(DUAL_F,dual_list_usami(:,5),1); %span it w.r.t. size of each firm
				
				sel						= dual_list_usami(:,3)';	
				DUAL_Flag				= [sparse(N_dual,N+J) 	 sparse((1:N_dual)',sel,1,N_dual,Jlag)    sparse(N_dual,K-N-J-Jlag)];	
				DUAL_Flag				= repelem(DUAL_Flag,dual_list_usami(:,5),1); %span it w.r.t. size of each firm	
				
end				