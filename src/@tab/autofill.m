function T = autofill(T,flags,opts)

arguments
    T
    flags = [0 Inf -1 10]
    opts.sl = [NaN NaN]
end

W1 = T(1:end-1,2);
W2 = T(2:end,2);

for f = flags
    switch f
    case -1                                                                         % L-R slope pairs
        mask = ~isnan([T(2:end,3) T(1:end-1,4)]);
        if any(mask)
            idxL = find(mask(:,1)) + 1;
            idxR = find(mask(:,2));
            T(idxL-1,4) = -T(idxL,3);
            T(idxR+1,3) = -T(idxR,4);
        end
    case 0                                                                          % Zero slopes
        idx = find(W2==W1);
        if ~isempty(idx)
            T(idx,4) = 0;
            T(idx+1,3) = 0;
        end
    case Inf                                                                        % +Inf slopes
        T(find(W1~=0 & W2==0 & ~isnan(W1)) + 1, 3) = Inf;
        T(find(W1==0 & W2~=0 & ~isnan(W2)), 4) = Inf;
    case 10                                                                         % Constant W 
        maskNB = isnan([W2 W1]);   
        mask = sum(maskNB,2)==1 & W2~=W1 & T(1:end-1,4)==0;                         % [nband x 1] one NaN band end / W1 or W2 given / zero RS
        mask = mask & maskNB;
        if any(mask)
            idxL = find(mask(:,2));
            idxR = find(mask(:,1));
            T(idxL,2) = T(idxL+1,2);
            T(idxR+1,2) = T(idxR,2);
            W1 = T(1:end-1,2);                                                      % Update on exiting this case
            W2 = T(2:end,2);
        end
    case 99                                                                         % Extrapolation slopes
        idx = size(T,1)*[2 4] + [1 0];                                              % extrap slope indices
        ids = size(T,1)*[3 3] + [1 0];                                              % adjacent band slope indices
        ixl = ~isnan(opts.sl);
        ixn = isnan(T(idx));
        T(idx(ixl)) = opts.sl(ixl);                                                 % non-NaN slopes specified in optional arguments
        T(idx(ixn & ~ixl)) = -T(ids(ixn & ~ixl));                                   % defaults to end band slope
        T(idx(T([1,end],2)==0)) = 0;                                                % overwrites with zero if W == 0
    otherwise
        continue
    end
end















