function ids = getblocks(T,tol)

arguments
    T
    tol = 0
end

idx = find(tab.symeq(T(2:end,1), T(1:end-1,1), tol));
ids = [[1; idx+1] [idx; size(T,1)]]; 