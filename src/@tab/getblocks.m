%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       ids = tab.getblocks(T,tol)
%       Computes tab.block property from tab.val given by T, default tol = 0
%
%       See also: symeq
%
%   OUTPUTS
%       [k x 2] array with with [start end] row indices of continuous blocks of T
%
%   VERSION
%       v1.0 / 26.10.22 / V.Yotov
%  ------------------------------------------------------------------------------------------------

function ids = getblocks(T,tol)

arguments
    T
    tol = 0
end

idx = find(tab.symeq(T(2:end,1), T(1:end-1,1), tol));
ids = [[1; idx+1] [idx; size(T,1)]]; 