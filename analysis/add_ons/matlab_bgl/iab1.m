%This program reads in the set of firm to firm transitions
%and computes the largest connected set of firms.
%It then saves the list of connected firmids to disk

path(path,'~/matlabBGL/');


for i=1:2
s=['transitions' int2str(i) '.csv'];
data=importdata(s);
data=data.data;

firmid=data(:,1);
lagfirmid=data(:,2);

A=sparse(lagfirmid',firmid',1); %adjacency matrix
[sindex, sz]=components(A); %get connected sets
idx=find(sz==max(sz)); %find largest set
firmlst=find(sindex==idx); %firms in connected set

s=['connectedfirms' int2str(i) '.txt'];
dlmwrite(s,firmlst);

end


