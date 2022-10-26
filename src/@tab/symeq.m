function b = symeq(v,w,tol) 

    if nargin < 3
        tol = 1e4;
    end
    
    b = abs(v-w) < tol * max( eps(), eps(min(abs(v),abs(w))) );