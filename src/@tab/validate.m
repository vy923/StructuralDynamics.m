%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       tab.validate(T,flags)
%       Various consistency checks for tab.val given by T
%
%       See also:       --
%       Related:        split, compute, autofill, collapse
%
%   INPUTS
%       T               [n x 4] in the form of tab.val
%
%       flags           [k x 1] flags of checks to perform, default = [1:5]
%           923         T contains NaNs
%           -Inf        zero Y and -Inf slope
%           Inf         nonzero Y and +Inf slope
%           1           nonzero slope in constant Y block
%           2           zero slope in nonconstant Y block
%           3           inconsistent Y left and right slopes
%           4           zero Y and noninf/nonnan slope -> nonzero/nonnan Y
%           5           zero Y1 and Y2 with noninf/nonnan slope
%
%   VERSION
%       v1.0 / 26.10.22 / V.Yotov
%  ------------------------------------------------------------------------------------------------

function validate(T,flags)

if nargin < 2
    flags = 1:5;
end

if any(~isinf(flags))
    LS = T(2:end,3);
    RS = T(1:end-1,4);
    Y1 = T(1:end-1,2);
    Y2 = T(2:end,2);
    maskNNS = ~isnan([LS RS]);                                                      % [nblock x 2] non-NaN slopes
    maskNB = isnan([Y2 Y1]);                                                        % [nblock x 2] NaN block ends
end

for f = flags
switch f
    case 923
        assert( ~any(isnan(T),'all'), ...
                "tab: incomplete/ambiguous input" ) 
    case -Inf 
        assert( ~any(T(:,3)==-Inf & T(:,2)==0), ...
                "tab: zero Y and -Inf slope" )
    case Inf
        assert( ~any(T(:,4)==Inf & T(:,2)~=0 & ~isnan(T(:,2))), ...
                "tab: nonzero Y and +Inf slope" )
    case 1
        assert( ~any(Y2==Y1 & any([LS RS]~=0 & maskNNS, 2) & Y2~=0), ...
                "tab: nonzero slope in constant Y block" )
    case 2
        assert( ~any(Y2~=Y1 & ~any(maskNB,2) & (LS==0 | RS==0)), ...
                "tab: zero slope in nonconstant Y block" )
    case 3 
        assert( ~any(sum([LS RS].*maskNNS, 2)), ...
                "tab: inconsistent Y left and right slopes" )
    case 4 
        assert( ~any([Y1 Y2]==0 & [RS LS]~=Inf & maskNNS(:,[2 1]) & [Y2 Y1]~=0 & ~maskNB, 'all'), ...
                "tab: zero Y and noninf/nonnan slope -> nonzero/nonnan Y" )
    case 5
        assert( ~any(all([Y1 Y2]==0,2) &  any(maskNNS & [LS RS]~=0, 2)), ...
                "tab: zero Y1 and Y2 with noninf/nonnan slope" )
    otherwise
        continue
end
end






