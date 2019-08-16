function [nod_adjele]=get_nod_adjele(NNODES,NELM,NCA)
%------------------------------------------------------------------------
%  Purpose:
%     find adjacent elements of each node
%
%  Synopsis:
%     [nod_adjele]=get_nod_adjele  
%
%  Variable Description:
%    nod_adjele - matrix containing adjacent elements of each node
%--------------------------------------------------------------------------
% Coded by Sandeep Kshirsagar                                              %
%--------------------------------------------------------------------------
nod_adjele = cell(NNODES,1);
for inn=1:NNODES   % loop for all nodes   
    inc=0;
    for ien=1:NELM  % loop for all elements
       if (find(inn==NCA(ien,1:3)))>=1   
           inc = inc+1;
           nod_adjele{inn}(inc) = ien; 
       end
    end
end