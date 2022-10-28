%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       T = tab.collapse(T,opts)
%       Removes superfluous rows of T, i.e. having constant Y / LS / RS
%
%       See also:       symeq
%       Related:        split, autofill, compute, validate
%
%   INPUTS
%       T               tab.val array
%       opts{:}
%           tol         default = 0
%           idx         default = [2 3 4], columns of T considered in collapse criterion
%
%   OUTPUTS
%       updated T
%
%   VERSION
%       v1.0 / 26.10.22 / V.Yotov
%  ------------------------------------------------------------------------------------------------

function T = collapse(T,opts)

arguments
    T
    opts.tol = 0
    opts.idx = [2 3 4]
end

% 3-row window equality check for columns opts.idx
mask = tab.symeq(T(2:end-1,opts.idx), T(1:end-2,opts.idx), opts.tol) & ...
       tab.symeq(T(2:end-1,opts.idx), T(3:end,opts.idx), opts.tol);

% Remove superfluous rows
if any(mask,'all')
    T(find(any(mask,2)) + 1, :) = [];
end