function T = autofill(T,flags,opts)
%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       T = tab.autofill(T,flags,opts)
%
%       See also:       --
%       Related:        compute, split, validate, collapse
%
%   INPUTS
%       T               [n x 4] in the form of tab.val
%       flags           [k x 1] flags what values to fill, default = [0 Inf -1 10]
%           -1          left-right slope pairs
%           0           zero slopes from constant Y
%           inf         +Inf slopes
%           10          constant Y from 0 slopes
%           99          extrapolation slopes
%
%       opts{:}
%           sl          [2 x 1] with extrapolation [LS RS], default = [NaN NaN]
%
%   OUTPUTS
%       T               [n x 4] with filled values according to 'flags'
%
%   VERSION
%   v1.0 / 26.10.22 / V.Y.
%  ------------------------------------------------------------------------------------------------

arguments
    T
    flags = [0 Inf -1 10]
    opts.sl = [NaN NaN]
end

Y1 = T(1:end-1,2);
Y2 = T(2:end,2);

for f = flags
switch f
    case -1 
        % L-R slope pairs
        mask = ~isnan([T(2:end,3) T(1:end-1,4)]);
        if any(mask)
            idxL = find(mask(:,1)) + 1;
            idxR = find(mask(:,2));
            T(idxL-1,4) = -T(idxL,3);
            T(idxR+1,3) = -T(idxR,4);
        end
    case 0     
        % Zero slopes
        idx = find(Y2==Y1);
        if ~isempty(idx)
            T(idx,4) = 0;
            T(idx+1,3) = 0;
        end
    case Inf     
        % +Inf slopes
        T([false; Y1~=0 & Y2==0 & ~isnan(Y1)], 3) = Inf;
        T(Y1==0 & Y2~=0 & ~isnan(Y2), 4) = Inf;
    case 10 
        % Constant Y                                                                         
        maskNB = isnan([Y2 Y1]);   
        mask = sum(maskNB,2)==1 & Y2~=Y1 & T(1:end-1,4)==0;                         % [nblock x 1] one NaN block end / Y1 or Y2 given / zero RS
        mask = mask & maskNB;
        if any(mask)
            idxL = find(mask(:,2));
            idxR = find(mask(:,1));
            T(idxL,2) = T(idxL+1,2);
            T(idxR+1,2) = T(idxR,2);
            Y1 = T(1:end-1,2);                                                      % update preallocated Y1/Y2 on exiting case 10
            Y2 = T(2:end,2);
        end
    case 99
        % Extrapolation slopes
        idx = size(T,1)*[2 4] + [1 0];                                              % extrap slope indices
        ids = size(T,1)*[3 3] + [1 0];                                              % adjacent block slope indices
        ixl = ~isnan(opts.sl);
        ixn = isnan(T(idx));
        T(idx(ixl)) = opts.sl(ixl);                                                 % non-NaN slopes specified in optional arguments
        T(idx(ixn & ~ixl)) = -T(ids(ixn & ~ixl));                                   % defaults to end block slope
        T(idx(T([1,end],2)==0)) = 0;                                                % overwrites with zero if Y == 0
    otherwise
        continue
end
end















