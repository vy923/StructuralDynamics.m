function validate(T,flags)

    if nargin < 2
        flags = [1:5];
    end
    
    if any(~isinf(flags))
        LS = T(2:end,3);
        RS = T(1:end-1,4);
        W1 = T(1:end-1,2);
        W2 = T(2:end,2);
        maskNNS = ~isnan([LS RS]);                                                      % [nband x 2] non-NaN slopes
        maskNB = isnan([W2 W1]);                                                        % [nband x 2] NaN band ends
    end

    for f = flags
        switch f
            case 923
                assert( ~any(isnan(T),'all'), ...
                        "tabc4: incomplete/ambiguous input" ) 
            case -Inf 
                assert( ~any(T(:,3)==-Inf & T(:,2)==0), ...
                        "tabc4: zero PSD and -Inf slope" )
            case Inf
                assert( ~any(T(:,4)==Inf & T(:,2)~=0 & ~isnan(T(:,2))), ...
                        "tabc4: nonzero PSD and +Inf slope" )
            case 1
                assert( ~any(W2==W1 & any([LS RS]~=0 & maskNNS, 2) & W2~=0), ...
                        "tabc4: nonzero slope in constant W band" )
            case 2
                assert( ~any(W2~=W1 & ~any(maskNB,2) & (LS==0 | RS==0)), ...
                        "tabc4: zero slope in nonconstant W band" )
            case 3 
                assert( ~any(sum([LS RS].*maskNNS, 2)), ...
                        "tabc4: inconsistent W left and right slopes" )
            case 4 
                assert( ~any([W1 W2]==0 & [RS LS]~=Inf & maskNNS(:,[2 1]) & [W2 W1]~=0 & ~maskNB, 'all'), ...
                        "tabc4: zero W and noninf/nonnan slope -> nonzero/nonnan W" )
            case 5
                assert( ~any(all([W1 W2]==0,2) &  any(maskNNS & [LS RS]~=0, 2)), ...
                        "tabc4: zero W1 and W2 with noninf/nonnan slope" )
            otherwise
                continue
        end
    end