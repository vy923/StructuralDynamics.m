function T = collapse(T,opts)
%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       T = tab.collapse(T,opts)
%       Removes superfluous rows of T, i.e. having constant Y / LS / RS
%
%       See also:       isapprox
%       Related:        split, autofill, compute, validate
%
%   INPUTS
%       T               tab.val array
%       opts{:}
%           atol        default = 0
%           idx         default = [2 3 4], columns of T considered in collapse criterion
%
%   OUTPUTS
%       updated T
%
%   VERSION
%       v1.0 / 26.10.22 / V.Yotov
%  ------------------------------------------------------------------------------------------------

arguments
    T
    opts.atol = 0
    opts.idx = [2 3 4]
end

% 3-row window equality check for columns opts.idx
mask = isapprox(T(2:end-1,opts.idx), T(1:end-2,opts.idx), opts.atol) & ...
       isapprox(T(2:end-1,opts.idx), T(3:end,opts.idx), opts.atol);

% Remove superfluous rows
if any(mask,'all')
    T(find(any(mask,2)) + 1, :) = [];
end