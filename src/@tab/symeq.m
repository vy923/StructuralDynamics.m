function mask = symeq(v,w,tol) 

if nargin < 3
    tol = 0;
end

if tol==0
    mask = v==w;
else
    mask = abs(v-w) <= tol*max(eps(), eps(min(abs(v),abs(w))));
end