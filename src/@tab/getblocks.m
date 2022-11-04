function ids = getblocks(T,atol)
%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       ids = tab.getblocks(T,atol)
%       Computes tab.block property from tab.val given by T, default atol = 0
%
%       See also: isapprox
%
%   OUTPUTS
%       [k x 2] array with with [start end] row indices of continuous blocks of T
%
%   VERSION
%   v1.0 / 26.10.22 / V.Y.
%  ------------------------------------------------------------------------------------------------

arguments
    T
    atol = 0
end

idx = find(isapprox(T(2:end,1), T(1:end-1,1), atol));
ids = [[1; idx+1] [idx; size(T,1)]]; 