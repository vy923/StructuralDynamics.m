%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       T = tab(T,opts)
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
%           collapse    
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
%       + add autofill in comp. loop for (isinf(RS) | isinf(LS) | W1.*W2==0)
%       X loop slope autofill and table split?
%
%   VERSION
%       v2.1 / 26.10.22 / --        redefined as a separate class in @tab
%       v2.0 / 23.10.22 / --        handling of zero W, +/-inf slopes, discontinuous tables
%       v1.1 / 16.10.22 / --        epstol, extrapolation slopes for end bands
%       v1.0 / 14.10.22 / V.Yotov
%  ------------------------------------------------------------------------------------------------

classdef tab

properties
    val
end

methods 
    
    function tb = tab(T,opts)    
        arguments
            T
            opts.type {mustBeMemberSCI(opts.type,["PSD"])} = "PSD"
            opts.split (1,1) {mustBeMember(opts.split,[0 1])} = true
            opts.collapse (1,1) {mustBeMember(opts.collapse,[0 1])} = true
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
        
        % Check / split / check / autofill
            tab.validate(T,[-Inf Inf]);

            if opts.split
                T = tab.split(T);
            end
            
            tab.validate(T, [1:5]);
            T = tab.autofill(T, [0 -Inf Inf -1 10]);                                        % Zero slopes -> Inf slopes -> LS-RS pairs -> const. W
        
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
        
                maskBand(isinf(RS)|isinf(LS)|W1.*W2==0) = false;
        
                maskR = ~isnan(W1 + RS) & maskBand;
                maskL = ~isnan(LS + W2) & maskBand;
                maskW = ~isnan(W1 + W2) & maskBand;
                assert(any(maskR|maskL|maskW));                                             % prevent infinite loops
        
                idxR = find(maskR);
                idxL = find(maskL);
                idxW = find(maskW);
        
                tmpW2 = W1(idxR) .* 10.^(RS(idxR).*logF(idxR)/10);                          % <--- allow for other computational rules in next versions
                tmpW1 = W2(idxL) .* 10.^(RS(idxL).*-logF(idxL)/10);
                tmpRS = 10*log10(W2(idxW)./W1(idxW)) ./ logF(idxW);
        
                assert( all([
                        tab.symeq(W2(maskR & ~isnan(W2.*maskR)), tmpW2(~isnan(W2(maskR))), opts.epstol)
                        tab.symeq(W1(maskL & ~isnan(W1.*maskL)), tmpW1(~isnan(W1(maskL))), opts.epstol)
                        tab.symeq(RS(maskW & ~isnan(RS.*maskW)), tmpRS(~isnan(RS(maskW))), opts.epstol)
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
        
                T = tab.autofill(T, [Inf -1]);                                              % "aggressive" autofill for +/-Inf slopes
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
            tab.validate(T, [923, 1:5, -Inf, Inf]);
            if opts.collapse
                T = tab.collapse(T,opts.epstol);
            end
    
        % Create object
            tb.val = T;
    end
end

methods(Static)
    T = split(T)
    T = autofill(T,flags)   
    T = collapse(T,tol,idx)
    b = symeq(v,w,tol) 
    validate(T,flags) 
end

end % CLASSDEF



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
%}
%  ------------------------------------------------------------------------------------------------










