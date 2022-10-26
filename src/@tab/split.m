function T = split(T)
    
    lvar = 0;
    while lvar ~= size(T,1)  
        lvar = size(T,1);
    
        % Prev RS ~= Inf, F ~= Prev F, Prev PSD ~= 0, PSD == 0
        mask = T(1:end-1,4)~=Inf & T(1:end-1,1)~=T(2:end,1) & T(1:end-1,2)~=0 & T(2:end,2)==0;                          
        if any(mask)
            idxN = (1:nnz(mask))' + find(mask);
            idxO = (1:size(T,1))' + [0; cumsum(mask)];
            idx  = [idxO; idxN];
            T    = mapSet(idx, 1:max(idx), [T; T(mask,1).*[1 0 Inf NaN]]);
        end
        
        % Prev RS ~= Inf, F ~= Prev F, Prev PSD ~= 0, LS == -Inf        
        mask = [NaN; T(1:end-1,4)]~=Inf & [NaN; T(1:end-1,1)]~=T(:,1) & [NaN; T(1:end-1,2)]~=0 & T(:,3)==-Inf;         
        if any(mask)
            idxN = (0:nnz(mask)-1)' + find(mask);
            idxO = (1:size(T,1))' + cumsum(mask); 
            idx  = [idxO; idxN];
            T    = mapSet(idx, 1:max(idx), [T; T(mask,1).*[1 0 NaN NaN]]);
        end
    
        % Next LS ~= -Inf, F ~= Next F, Next PSD ~= 0, PSD == 0
        mask = T(2:end,4)~=-Inf & T(2:end,1)~=T(1:end-1,1) & T(2:end,2)~=0 & T(1:end-1,2)==0;                           
        if any(mask)
            idxN = [0:nnz(mask)-1]' + find(mask) + 1;
            idxO = setdiff(1:size(T,1)+nnz(mask), idxN)';
            idx  = [idxO; idxN];
            T    = mapSet(idx, 1:max(idx), [T; T(find(mask)+1,1) .* [1 0 NaN Inf]]);
        end
    
        % Next LS ~= -Inf, F ~= Next F, Next PSD ~= 0, RS ==- Inf
        mask = [T(2:end,4); NaN]~=-Inf & [T(2:end,1); NaN]~=T(:,1) & [T(2:end,2); NaN]~=0 & T(:,4)==-Inf;               
        if any(mask)
            idxN = [0:nnz(mask)-1]' + find(mask) + 1;
            idxO = setdiff(1:size(T,1)+nnz(mask), idxN)';
            idx  = [idxO; idxN];
            T    = mapSet(idx, 1:max(idx), [T; T(mask,1) .* [1 0 NaN NaN]]);
        end
    end 