%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       mask = tab.symeq(v,w,tol)
%
%       Check for approximate equality between arrays v, w
%       Default eps tolerance is 1e4
%
%   OUTPUTS
%       boolean array |v-w| <= tol*eps(...)
%
%   VERSION
%       v1.0 / 26.10.22 / V.Yotov
%  ------------------------------------------------------------------------------------------------

function mask = symeq(v,w,tol) 

arguments
    v
    w
    tol = 1e4
end

% Direct comparison if tolerance is zero
if tol==0
    mask = v==w;
else
    mask = abs(v-w) <= tol*max(eps(), eps(min(abs(v),abs(w))));
end