function T = collapse(T,tol,idx)

    if nargin < 3
        idx = [2 3 4];
    end

    mask = tab.symeq(T(2:end-1,idx), T(1:end-2,idx), tol) & tab.symeq(T(2:end-1,idx), T(3:end,idx), tol); 
    T(find(any(mask,2)) + 1, :) = [];