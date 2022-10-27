function T = collapse(T,opts)

arguments
    T
    opts.tol = 0
    opts.idx = [2 3 4]
end

mask = tab.symeq(T(2:end-1,opts.idx), T(1:end-2,opts.idx), opts.tol) & ...
       tab.symeq(T(2:end-1,opts.idx), T(3:end,opts.idx), opts.tol);

if any(mask,'all')
    T(find(any(mask,2)) + 1, :) = [];
end