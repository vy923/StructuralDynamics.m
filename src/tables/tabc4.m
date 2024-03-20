%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       T = tabc4(T,opts)
%
%       See also:       mapSet, mustBeMemberSCI
%       Related:        interpc4, scalec4, integratec4, opsc4
%
%   INPUTS
%       T               [nband + 1 x c] table with c = 2/3/4 and NaN empty fields
%           c == 2      assumes [f W]
%           c == 3      assumes [f W RS]
%           c == 4      must be [f W LS RS]
%
%       opts{:}
%           type        table type, currently only "PSD"
%               "PSD"   default, assumes log-log and dB/oct slope units
%
%           split       autosplit table on discontinuities
%               true    default, handles +/-inf, discontinuities 
%               false   assumes double points are given explicitly
%
%           compress    
%               true    default
%               false   allows redundant points, e.g. midpoint of constant band
%
%           epstol      eps tolerance scaling for comparisons, default = 1E+4  
%           ls          left slope of first band
%           rs          right slope of final band
%
%   OUTPUTS
%       T               completed and validated [nband + 1 x 4] table
%
%   UPDATES
%       - combine table i.e. with compress + sort
%       - merge function with min/max options for adding profiles, etc.
%       - add autofill in comp. loop for (isinf(RS) | isinf(LS) | W1.*W2==0)
%       - loop slope autofill and table split??
%
%   VERSION
%       v2.1 / 21.03.24 / --        allow zero PSD blocks
%       v2.0 / 23.10.22 / --        handling of zero PSD, +/-inf slopes, discontinuous tables
%       v1.1 / 16.10.22 / --        epstol, extrapolation slopes for end bands
%       v1.0 / 14.10.22 / V.Yotov
%  ------------------------------------------------------------------------------------------------

function T = tabc4(T,opts)

arguments
    T
    opts.type {mustBeMemberSCI(opts.type,["PSD"])} = "PSD"
    opts.split (1,1) {mustBeMember(opts.split,[0 1])} = true
    opts.compress (1,1) {mustBeMember(opts.compress,[0 1])} = true
    opts.epstol (1,1) {mustBePositive} = 1E+4
    opts.ls (1,1) {mustBeReal} = NaN
    opts.rs (1,1) {mustBeReal} = NaN
end

% Initialisation, sorting
switch size(T,2)
    case 2; T = [T T*NaN];
    case 3; T = [T(:,1:2) T(:,3)+[NaN 0]]; 
    case 4 
    otherwise error("tabc4: table must have 2 to 4 columns")
end

if ~issorted(T(:,1))
    h = msgbox(sprintf(['First column of input was not sorted. \nCorrect results ' ...
        'are not guaranteed in case of discontinuous tables.']), 'tabc4', 'warn');
    T = sortrows(T);
end

% Discontinuous-specific checks
    validate(T,[-Inf Inf]);

% Split discontinuous tables
    % timeit(@() splitDiscontinuous(T))
    if opts.split
        T = splitDiscontinuous(T);
    end

% Consistency check
    validate(T, [1:5]);

% Direct completions 
    %timeit(@() autofill(T))
    T = autofill(T, [0 -Inf Inf -1 10]);                                             % Zero slopes -> Inf slopes -> LS-RS pairs -> const. W

% Computation
    F1 = T(1:end-1,1);
    F2 = T(2:end,1);
    logF = log10(F2./F1)/log10(2);
    
    maskBand = true(size(T,1)-1, 1);

    while any(maskBand)
        LS = T(2:end,3); 
        RS = T(1:end-1,4); 
        W1 = T(1:end-1,2);
        W2 = T(2:end,2);

        maskBand(isinf(RS)|isinf(LS)|xor(W1==0,W2==0)) = false;

        maskR = ~isnan(W1 + RS) & maskBand;
        maskL = ~isnan(LS + W2) & maskBand;
        maskW = ~isnan(W1 + W2) & maskBand & (W1.*W2~=0);
        assert(any(maskR|maskL|maskW));                                             % prevent infinite loops

        idxR = find(maskR);
        idxL = find(maskL);
        idxW = find(maskW);

        tmpW2 = W1(idxR) .* 10.^(RS(idxR).*logF(idxR)/10);                          % <--- allow for other computational rules in next versions
        tmpW1 = W2(idxL) .* 10.^(RS(idxL).*-logF(idxL)/10);
        tmpRS = 10*log10(W2(idxW)./W1(idxW)) ./ logF(idxW);

        assert( all([
                symeq(W2(maskR & ~isnan(W2.*maskR)), tmpW2(~isnan(W2(maskR))), opts.epstol)
                symeq(W1(maskL & ~isnan(W1.*maskL)), tmpW1(~isnan(W1(maskL))), opts.epstol)
                symeq(RS(maskW & ~isnan(RS.*maskW)), tmpRS(~isnan(RS(maskW))), opts.epstol)
                ]), ...
            "tabc4: overspecified input");

        maskRS = isnan(T(idxR+1,2));
        maskLS = isnan(T(idxL,2));
        maskWR = isnan(T(idxW,4)); 
        maskWL = isnan(T(idxW+1,3)); 

        T(idxR(maskRS)+1,2) = tmpW2(maskRS);
        T(idxL(maskLS),2) = tmpW1(maskLS);
        T(idxW(maskWR),4) = tmpRS(maskWR);
        T(idxW(maskWL)+1,3) = -tmpRS(maskWL);

        T = autofill(T, [Inf -1]);                                                  % "aggressive" autofill for +/-Inf slopes
        maskBand(maskR|maskL|maskW) = false;
    end

% Fill extrapolation slopes
    idx = size(T,1)*[2 4] + [1 0];                                                  % extrapolation slope indices
    ids = size(T,1)*[3 3] + [1 0];                                                  % adjacent band slope indices
    ops = [opts.ls opts.rs];
    ixl = ~isnan(ops);
    ixn = isnan(T(idx));

    T(idx(ixl)) = ops(ixl);                                                         % non-NaN slopes specified in optional arguments
    T(idx(ixn & ~ixl)) = -T(ids(ixn & ~ixl));                                       % defaults to end band slope
    T(idx(T([1,end],2)==0)) = 0;                                                    % overwrites with zero if W == 0

% Check completeness and clean-up
    validate(T, [923, 1:5, -Inf, Inf]);

    %timeit(@() compress(T,opts.epstol))
    if opts.compress
        T = compress(T,opts.epstol);
    end



%  ------------------------------------------------------------------------------------------------

function b = symeq(v,w,tol) 
    if nargin < 3
        tol = 1e4;
    end
    b = abs(v-w) < tol * max( eps(), eps(min(abs(v),abs(w))) );

%  ------------------------------------------------------------------------------------------------

function T = compress(T,tol,idx)
    if nargin < 3
        idx = [2 3 4];
    end
    mask = symeq(T(2:end-1,idx), T(1:end-2,idx), tol) & symeq(T(2:end-1,idx), T(3:end,idx), tol); 
    T(find(any(mask,2)) + 1, :) = [];

%  ------------------------------------------------------------------------------------------------

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
 
%  ------------------------------------------------------------------------------------------------

function T = autofill(T,flags) 

    if nargin < 2
        flags = [0 Inf -1 10];
    end
    
    W1 = T(1:end-1,2);
    W2 = T(2:end,2);
    
    for f = flags
    switch f
        case -1                                                                         % L-R slope pairs
            mask = ~isnan([T(2:end,3) T(1:end-1,4)]);
            idxL = find(mask(:,1)) + 1;
            idxR = find(mask(:,2));
            T(idxL-1,4) = -T(idxL,3);
            T(idxR+1,3) = -T(idxR,4);
        case 0                                                                          % Zero slopes
            band = find(W2==W1);
            T(band,4) = 0;
            T(band+1,3) = 0;
        case Inf                                                                        % +Inf slopes
            T(find(W1~=0 & W2==0 & ~isnan(W1)) + 1, 3) = Inf;
            T(find(W1==0 & W2~=0 & ~isnan(W2)), 4) = Inf;
        case 10                                                                         % Constant W 
            maskNB = isnan([W2 W1]);   
            mask = sum(maskNB,2)==1 & W2~=W1 & T(1:end-1,4)==0;                         % [nband x 1] one NaN band end / W1 or W2 given / zero RS
            mask = mask & maskNB;
            idxL = find(mask(:,2));
            idxR = find(mask(:,1));
            T(idxL,2) = T(idxL+1,2);
            T(idxR+1,2) = T(idxR,2);
            W1 = T(1:end-1,2);                                                          % Update on exiting this case
            W2 = T(2:end,2);
        case 923                                                                        % Computable NaNs from +/-Inf slopes
                                                                                        % < --------------------------------------
        otherwise
            continue
    end
    end

%  ------------------------------------------------------------------------------------------------

function T = splitDiscontinuous(T)
    
    lvar = 0;
    while lvar ~= size(T,1)  
        lvar = size(T,1);
    
        % Prev RS ~= Inf, F ~= Prev F, Prev PSD ~= 0, PSD == 0
        mask = T(1:end-1,4)~=Inf & T(1:end-1,1)~=T(2:end,1) & T(1:end-1,2)~=0 & T(2:end,2)==0;                          
        if any(mask)
            idxN = (1:nnz(mask))' + find(mask);
            idxO = (1:size(T,1))' + [0; cumsum(mask)];
            idx  = [idxO; idxN];
            T    = mapSet(idx, 1:max(idx), [T; T(mask,1).*[1 0 Inf NaN]]);
        end
        
        % Prev RS ~= Inf, F ~= Prev F, Prev PSD ~= 0, LS == -Inf        
        mask = [NaN; T(1:end-1,4)]~=Inf & [NaN; T(1:end-1,1)]~=T(:,1) & [NaN; T(1:end-1,2)]~=0 & T(:,3)==-Inf;         
        if any(mask)
            idxN = (0:nnz(mask)-1)' + find(mask);
            idxO = (1:size(T,1))' + cumsum(mask); 
            idx  = [idxO; idxN];
            T    = mapSet(idx, 1:max(idx), [T; T(mask,1).*[1 0 NaN NaN]]);
        end
    
        % Next LS ~= -Inf, F ~= Next F, Next PSD ~= 0, PSD == 0
        mask = T(2:end,4)~=-Inf & T(2:end,1)~=T(1:end-1,1) & T(2:end,2)~=0 & T(1:end-1,2)==0;                           
        if any(mask)
            idxN = [0:nnz(mask)-1]' + find(mask) + 1;
            idxO = setdiff(1:size(T,1)+nnz(mask), idxN)';
            idx  = [idxO; idxN];
            T    = mapSet(idx, 1:max(idx), [T; T(find(mask)+1,1) .* [1 0 NaN Inf]]);
        end
    
        % Next LS ~= -Inf, F ~= Next F, Next PSD ~= 0, RS ==- Inf
        mask = [T(2:end,4); NaN]~=-Inf & [T(2:end,1); NaN]~=T(:,1) & [T(2:end,2); NaN]~=0 & T(:,4)==-Inf;               
        if any(mask)
            idxN = [0:nnz(mask)-1]' + find(mask) + 1;
            idxO = setdiff(1:size(T,1)+nnz(mask), idxN)';
            idx  = [idxO; idxN];
            T    = mapSet(idx, 1:max(idx), [T; T(mask,1) .* [1 0 NaN NaN]]);
        end
    end 

%  ------------------------------------------------------------------------------------------------

%{
% Example 1, chained
T = [     20	    nan	        nan          +6	
          50        nan         nan          0
         800	    nan	        nan          nan
        2000        0.026       +6	         nan     ];
tabc4(T)

% Example 2, discontinuous
T = [     32          15        -Inf         NaN
          40          94         NaN        -Inf
          50          69         NaN         NaN
          55          77         NaN         NaN
          63           0         0           NaN
          80          86        -Inf         NaN
          80           0         NaN         NaN
         100           0         NaN         NaN
         125           0         0           Inf
         125          10         NaN         NaN
         160         108         NaN         NaN
         400           0         NaN         NaN
        1000           4         NaN         NaN
        1600           1        -Inf         NaN
        2000           1         NaN         NaN
        2500           4        -Inf        -Inf   ];
tabc4(T)

% Example 3, chained, non-sorted, compression, zeros
T = [     20	    nan	        nan          +6	
          50        nan         nan          0
         250        .05         nan          7
         310        nan         nan          nan   
         320        0           nan          nan
         350        nan         nan          3
         380        nan         nan          3
         400        0.100       nan          nan
         500        nan         -1           nan
         550        nan         -1           nan
         590        nan         -1           nan
         800	    nan	        nan          nan
        2000        0.026       +6	         nan     ];
tabc4(T,compress=false)
tabc4(tabc4(ans))

%  ------------------------------------------------------------------------------------------------

% [OBSOLETE CODE] Validation flags 1:5
    if any(W2==W1 & any([LS RS]~=0 & maskNNS, 2) & W2~=0)                           % const nonzero W band with nonzero RS or LS
        error("tabc4: nonzero slope in constant W band")
    elseif any(W2~=W1 & ~any(maskNB,2) & (LS==0 | RS==0))                           % not const W / numeric W1 and W2 / zero slope
        error("tabc4: zero slope in nonconstant W band")
    elseif any(sum([LS RS].*maskNNS, 2))                                            % inconsistent slopes at band ends
        error("tabc4: inconsistent W left and righ slopes")
    elseif any([W1 W2]==0 & [RS LS]~=Inf & maskNNS(:,[2 1]) & [W2 W1]~=0 & ~maskNB, 'all')  % was ~isinf([RS LS])
        error("tabc4: zero W and noninf/nonnan slope -> nonzero/nonnan W")          
    elseif any(all([W1 W2]==0,2) &  any(maskNNS & [LS RS]~=0, 2))                   % zero W -> zero W with noninf/nonnan slope
        error("tabc4: zero W1 and W2 with non-inf/non-NaN slope")
    end  

% [OBSOLETE CODE] Direct reassignment of computed values, even if non-nan
    T(idxR+1,2) = tmpW2;
    T(idxL,2)   = tmpW1;
    T(idxW,4)   = tmpRS;
    T(idxW+1,3) = -tmpRS;
%}












