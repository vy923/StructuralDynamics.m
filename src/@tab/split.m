function T = split(T)
%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       T = tab.split(T)
%       Adds rows to tab.val given by T such that discontinuities are properly represented  
%
%       See also:       --
%       Related:        collapse, autofill, compute, validate
%
%   VERSION
%   v1.1 / 28.10.22 / --    better performance, verbose replacements for mapSet/setdiff 
%   v1.0 / 26.10.22 / V.Y.
%  ------------------------------------------------------------------------------------------------
   
% Compute until no change in nrows of T occurs
lvar = 0;

while lvar ~= size(T,1) 
    lvar = size(T,1);

    % Prev RS~=Inf, X~=Prev X, Prev Y~=0, Y==0
    mask = T(1:end-1,4)~=Inf & T(1:end-1,1)~=T(2:end,1) & T(1:end-1,2)~=0 & T(2:end,2)==0;                          
    if any(mask)
        idxN = (1:nnz(mask))' + find(mask);
        idxO = (1:size(T,1))' + [0; cumsum(mask)];
        T = expandLoc(T, mask, idxO, idxN, [1 0 Inf NaN]);
    end
    
    % Prev RS~=Inf, X~=Prev X, Prev Y~=0, LS==-Inf        
    mask = [NaN; T(1:end-1,4)]~=Inf & [NaN; T(1:end-1,1)]~=T(:,1) & [NaN; T(1:end-1,2)]~=0 & T(:,3)==-Inf;         
    if any(mask)
        idxN = (0:nnz(mask)-1)' + find(mask);
        idxO = (1:size(T,1))' + cumsum(mask); 
        T = expandLoc(T, mask, idxO, idxN, [1 0 NaN NaN]);
    end

    % Next LS~=-Inf, X~=Next X, Next Y~=0, Y==0
    mask = T(2:end,4)~=-Inf & T(2:end,1)~=T(1:end-1,1) & T(2:end,2)~=0 & T(1:end-1,2)==0;                           
    if any(mask)
        idxN = [0:nnz(mask)-1]' + find(mask) + 1;
        idxO = 1:(size(T,1) + nnz(mask));
        idxO(idxN) = [];
        T = expandLoc(T, [false; mask], idxO, idxN, [1 0 NaN Inf]);
    end

    % Next LS~=-Inf, X~=Next X, Next Y~=0, RS==-Inf
    mask = [T(2:end,4); NaN]~=-Inf & [T(2:end,1); NaN]~=T(:,1) & [T(2:end,2); NaN]~=0 & T(:,4)==-Inf;               
    if any(mask)
        idxN = [0:nnz(mask)-1]' + find(mask) + 1;
        idxO = 1:(size(T,1) + nnz(mask));
        idxO(idxN) = [];
        T = expandLoc(T, mask, idxO, idxN, [1 0 NaN NaN]);
    end
end 

%  ------------------------------------------------------------------------------------------------

function T = expandLoc(T,mask,idxO,idxN,row)

    Q = nan(size(T));
    Q(idxO,:) = T;
    Q(idxN,:) = T(mask,1).*row;
    T = Q;

%  ------------------------------------------------------------------------------------------------













