%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       T = tablePSD(T,opts)
%
%       See also:       --
%       Related:        interpPSD, scalePSD, integratePSD
%
%   INPUTS
%       T               [nband + 1 x c] table with c = 2/3/4 and NaN empty fields
%           c == 2      assumes [f W]
%           c == 3      assumes [f W RS]
%           c == 4      must be [f W LS RS]
%       opts{:}
%           epstol      eps tolerance scaling for comparisons, default = 1E+4   
%           ls          left slope of first band
%           rs          right slope of final band
%
%   OUTPUTS
%       T               completed and validated [nband + 1 x 4] table
%
%   UPDATES
%       - comparison handling for +/-Inf values
%
%   VERSION
%       v1.1 / 16.10.22 / V.Yotov   epstol, extrapolation slopes for end bands
%       v1.0 / 14.10.22 / V.Yotov
%  ------------------------------------------------------------------------------------------------

function T = tablePSD(T,opts)

arguments
    T
    opts.epstol (1,1) {mustBePositive} = 1E+4
    opts.ls (1,1) {mustBeReal} = NaN
    opts.rs (1,1) {mustBeReal} = NaN
end

switch size(T,2)
    case 2; T = [T T*NaN];
    case 3; T = [T(:,1:2) T(:,3)+[NaN 0]]; 
    case 4 
    otherwise error("interpPSD :: PSD must have 2 to 4 columns")
end

% Constants 
    nband = size(T,1)-1;
    F1 = T(1:end-1,1);
    F2 = T(2:end,1);
    LF = log10(F2./F1)/log10(2);

% Consistency check
    LS = T(2:end,3);
    RS = T(1:end-1,4);

    maskNNS = ~isnan([LS RS]);                                                  % [nband x 2] non-NaN slopes
    maskCW = T(2:end,2)==T(1:end-1,2);                                          % [nband x 1] constant PSD bands
    maskNB = isnan([T(2:end,2) T(1:end-1,2)]);                                  % [nband x 2] NaN band ends              
    
    if any(maskCW & any([LS RS]~=0 & maskNNS, 2))                               % const W band with nonzero RS or LS 
        error("interpPSD :: nonzero slope in constant PSD band")
    elseif any(~maskCW & ~any(maskNB,2) & (LS==0 | RS==0))                      % not const W / numeric W1 and W2 / zero slope
        error("interpPSD :: zero slope in nonconstant PSD band")
    elseif any(sum([LS RS].*maskNNS, 2))                                        % inconsistent slopes at band ends
        error("interpPSD :: inconsistent PSD left and righ slopes")
    end

% Zero slopes
    band = find(maskCW);
    T(band,4) = 0;
    T(band+1,3) = 0;

% L-R slopes 
    mask = ~isnan([T(2:end,3) T(1:end-1,4)]);
    idxL = find(mask(:,1)) + 1;
    idxR = find(mask(:,2));
    T(idxL-1,4) = -T(idxL,3);
    T(idxR+1,3) = -T(idxR,4);

% Constant PSDs
    mask = sum(maskNB,2)==1 & ~maskCW & T(1:end-1,4)==0;                        % [nband x 1] one NaN band end / W1 or W2 given / zero RS
    mask = mask & maskNB;
    idxL = find(mask(:,2));
    idxR = find(mask(:,1));
    T(idxL,2) = T(idxL+1,2);
    T(idxR+1,2) = T(idxR,2);

% Computation
    maskBand = ones(nband,1);
    symeq = @(v,w) abs(v-w) < opts.epstol * eps(min(abs(v),abs(w)));            % check for equality within tol*eps(...)

    while any(maskBand)
        LS = T(2:end,3); 
        RS = T(1:end-1,4); 
        W1 = T(1:end-1,2);
        W2 = T(2:end,2);

        maskR = ~isnan(W1 + RS) & maskBand;
        maskL = ~isnan(LS + W2) & maskBand;
        maskW = ~isnan(W1 + W2) & maskBand;
        assert(any(maskR|maskL|maskW));                                         % prevent infinite loops, should be redundant

        idxR = find(maskR);
        idxL = find(maskL);
        idxW = find(maskW);

        tmpW2 = W1(idxR) .* 10.^(RS(idxR).*LF(idxR)/10);
        tmpW1 = W2(idxL) .* 10.^(RS(idxL).*-LF(idxL)/10);
        tmpRS = 10*log10(W2(idxW)./W1(idxW)) ./ LF(idxW);

        assert( all([
                symeq( W2(maskR & ~isnan(W2.*maskR)), tmpW2(~isnan(W2(maskR))) )
                symeq( W1(maskL & ~isnan(W1.*maskL)), tmpW1(~isnan(W1(maskL))) )
                symeq( RS(maskW & ~isnan(RS.*maskW)), tmpRS(~isnan(RS(maskW))) )   ...
                ]), ...
            "interpPSD :: overspecified PSD input");

        T(idxR+1,2) = tmpW2;
        T(idxL,2)   = tmpW1;
        T(idxW,4)   = tmpRS;
        T(idxW+1,3) = -tmpRS;

        maskBand(maskR|maskL|maskW) = 0;
    end

% Extrapolation slopes
    idx = (nband+1)*[2 4] + [1 0];                                              % extrapolation slope indices
    ids = (nband+1)*[3 3] + [1 0];                                              % adjacent band slope indices
    ops = [opts.ls opts.rs];
    ixl = ~isnan(ops);
    ixx = isnan(T(idx));
    T(idx(ixl)) = ops(ixl);                                                     % non-NaN slopes specified in optional arguments
    T(idx(ixx & ~ixl)) = -T(ids(ixx & ~ixl));                                   % defaults to end band slope                                         

% Test if table is full
    assert(~any(isnan(T),'all'),  "interpPSD :: incomplete PSD input")

    

%  ------------------------------------------------------------------------------------------------
%{
PSD =  [    20	        nan	        nan          +6	
            50          nan         nan          0
            800	        nan	        nan          nan
            2000        0.026       +6	         nan        ];
tablePSD(PSD)
%}












