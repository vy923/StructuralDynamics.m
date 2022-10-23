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
%           optimise    
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
%       - optimise=0/1 to remove bands with const. PSD/LS/RS
%       - combine table i.e. with optimise + sort
%       - merge function with min/max options for adding profiles, etc.
%       - computational rule affecting ends check?
%       - update names to interp4c, scale4c,...
%
%   VERSION
%       v2.0 / 23.10.22 / --        handling of zero PSD, +/-inf slopes,  discontinuous tables
%       v1.1 / 16.10.22 / --        epstol, extrapolation slopes for end bands
%       v1.0 / 14.10.22 / V.Yotov
%  ------------------------------------------------------------------------------------------------

function T = tabc4(T,opts)

arguments
    T
    opts.type {mustBeMemberSCI(opts.type,["PSD"])} = "PSD"
    opts.split (1,1) {mustBeMember(opts.split,[0 1])} = true
    opts.optimise (1,1) {mustBeMember(opts.optimise,[0 1])} = true
    opts.epstol (1,1) {mustBePositive} = 1E+4
    opts.ls (1,1) {mustBeReal} = NaN
    opts.rs (1,1) {mustBeReal} = NaN
end

% Initialisation, sorting
switch size(T,2)
    case 2; T = [T T*NaN];
    case 3; T = [T(:,1:2) T(:,3)+[NaN 0]]; 
    case 4 
    otherwise error("tabc4: PSD must have 2 to 4 columns")
end

if ~issorted(T(:,1))
    h = msgbox(sprintf(['First column of input was not sorted. \nCorrect results ' ...
        'are not guaranteed in case of discontinuous tables.']), 'tabc4','warn');
    T = sortrows(T);
end

% Discontinuity-specific checks
    if any(T(:,3)==-Inf & T(:,2)==0)                                                % zero PSD, -inf slope
        error("tabc4: zero PSD and -Inf slope")
    elseif any(T(:,4)==Inf & T(:,2)~=0)                                             % nonzero PSD, +inf slope
        error("tabc4: nonzero PSD and +Inf slope")
    end

% Split discontinuous tables
    if opts.split
        T = splitDiscontinuous(T);
    end

% Consistency check
    LS = T(2:end,3);
    RS = T(1:end-1,4);
    W1 = T(1:end-1,2);
    W2 = T(2:end,2);

    maskNNS = ~isnan([LS RS]);                                                      % [nband x 2] non-NaN slopes
    maskNB = isnan([W2 W1]);                                                        % [nband x 2] NaN band ends

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

% Fill zero slopes
    band = find(W2==W1);
    T(band,4) = 0;
    T(band+1,3) = 0;

% Fill +/-inf slopes
    W1 = T(1:end-1,2);
    W2 = T(2:end,2);
    T(find(W1~=0 & W2==0 & ~isnan(W1)) + 1, 3) = Inf;
    T(find(W1==0 & W2~=0 & ~isnan(W2)), 4) = -Inf;

% Fill L-R slopes 
    mask = ~isnan([T(2:end,3) T(1:end-1,4)]);
    idxL = find(mask(:,1)) + 1;
    idxR = find(mask(:,2));
    T(idxL-1,4) = -T(idxL,3);
    T(idxR+1,3) = -T(idxR,4);

% Fill constant Ws
    mask = sum(maskNB,2)==1 & W2~=W1 & T(1:end-1,4)==0;                             % [nband x 1] one NaN band end / W1 or W2 given / zero RS
    mask = mask & maskNB;
    idxL = find(mask(:,2));
    idxR = find(mask(:,1));
    T(idxL,2) = T(idxL+1,2);
    T(idxR+1,2) = T(idxR,2);

% Computation
    clearvars mask* idx*
    symeq = @(v,w) abs(v-w) < opts.epstol * eps(min(abs(v),abs(w)));                % check for equality within tol*eps(...)

    F1 = T(1:end-1,1);
    F2 = T(2:end,1);
    logF = log10(F2./F1)/log10(2);
    nband = size(T,1)-1;

    maskBand = true(nband,1);

    while any(maskBand)
        LS = T(2:end,3); 
        RS = T(1:end-1,4); 
        W1 = T(1:end-1,2);
        W2 = T(2:end,2);

        maskBand(isinf(RS) | isinf(LS) | W1.*W2==0) = false;

        maskR = ~isnan(W1 + RS) & maskBand;
        maskL = ~isnan(LS + W2) & maskBand;
        maskW = ~isnan(W1 + W2) & maskBand;
        assert(any(maskR|maskL|maskW));                                             % prevent infinite loops, should be redundant

        idxR = find(maskR);
        idxL = find(maskL);
        idxW = find(maskW);

        tmpW2 = W1(idxR) .* 10.^(RS(idxR).*logF(idxR)/10);                          % <--- allow for other computational rules!!
        tmpW1 = W2(idxL) .* 10.^(RS(idxL).*-logF(idxL)/10);
        tmpRS = 10*log10(W2(idxW)./W1(idxW)) ./ logF(idxW);

        assert( all([
                symeq( W2(maskR & ~isnan(W2.*maskR)), tmpW2(~isnan(W2(maskR))) )
                symeq( W1(maskL & ~isnan(W1.*maskL)), tmpW1(~isnan(W1(maskL))) )
                symeq( RS(maskW & ~isnan(RS.*maskW)), tmpRS(~isnan(RS(maskW))) )   ...
                ]), ...
            "tabc4: overspecified input");

        T(idxR+1,2) = tmpW2;
        T(idxL,2)   = tmpW1;
        T(idxW,4)   = tmpRS;
        T(idxW+1,3) = -tmpRS;

        maskBand(maskR|maskL|maskW) = false;
    end

% Fill extrapolation slopes
    idx = (nband+1)*[2 4] + [1 0];                                                  % extrapolation slope indices
    ids = (nband+1)*[3 3] + [1 0];                                                  % adjacent band slope indices
    ops = [opts.ls opts.rs];
    ixl = ~isnan(ops);
    ixn = isnan(T(idx));

    T(idx(ixl)) = ops(ixl);                                                         % non-NaN slopes specified in optional arguments
    T(idx(ixn & ~ixl)) = -T(ids(ixn & ~ixl));                                       % defaults to end band slope
    T(idx(T([1,end],2)==0)) = 0;                                                    % overwrites with zero if W == 0

% Test table
    assert(~any(isnan(T),'all'), "tabc4: incomplete input")
    


%  ------------------------------------------------------------------------------------------------

function T = splitDiscontinuous(T)

lvar = 0;
while lvar ~= size(T,1)  
    lvar = size(T,1);
           
    % Condition: (Prev RS ~= Inf) and (F ~= Prev F) and (Prev PSD ~= 0) and (PSD == 0)
        mask = T(1:end-1,4)~=Inf & T(1:end-1,1)~=T(2:end,1) & T(1:end-1,2)~=0 & T(2:end,2)==0;
    if any(mask)
        idxN = (1:nnz(mask))' + find(mask);
        idxO = (1:size(T,1))' + [0; cumsum(mask)];
        idx  = [idxO; idxN];
        T    = mapSet(idx, 1:max(idx), [T; T(mask,1).*[1 0 Inf NaN]]);
    end
            
    % Condition: (Prev RS ~= Inf) and (F ~= Prev F) and (Prev PSD ~= 0) and (LS == -Inf)
        mask = [NaN; T(1:end-1,4)]~=Inf & [NaN; T(1:end-1,1)]~=T(:,1) & [NaN; T(1:end-1,2)]~=0 & T(:,3)==-Inf; 
    if any(mask)
        idxN = (0:nnz(mask)-1)' + find(mask);
        idxO = (1:size(T,1))' + cumsum(mask); 
        idx  = [idxO; idxN];
        T    = mapSet(idx, 1:max(idx), [T; T(mask,1).*[1 0 NaN NaN]]);
    end

    % Condition: (Next LS ~= -Inf) and (F ~= Next F) and (Next PSD ~= 0) and (PSD == 0)
        mask = T(2:end,4)~=-Inf & T(2:end,1)~=T(1:end-1,1) & T(2:end,2)~=0 & T(1:end-1,2)==0;
    if any(mask)
        idxN = [0:nnz(mask)-1]' + find(mask) + 1;
        idxO = setdiff(1:size(T,1)+nnz(mask), idxN)';
        idx  = [idxO; idxN];
        T    = mapSet(idx, 1:max(idx), [T; T(find(mask)+1,1) .* [1 0 NaN Inf]]);
    end
        
    % Condition: (Next LS ~= -Inf) and (F ~= Next F) and (Next PSD ~= 0) and (RS == -Inf)
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
%}












