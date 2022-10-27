%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       T = tab(T,opts)
%
%       See also:       mapSet, mustBeMemberSCI
%       Related:        tab.interp, tab.scale, tab.integrate, tab.ops
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
%       - update help file
%       X combine table i.e. with compress + sort
%       - merge function with min/max options for adding profiles, etc.
%       + separate function for computing block indices
%       + move extrapolation slopes to autofill
%       + add autofill in comp. loop for (isinf(RS) | isinf(LS) | W1.*W2==0)
%       X loop slope autofill and table split?
%
%   VERSION
%       v2.0 / 26.10.22 / --        classdef in @tab, updated static methods for the constructor
%       v1.2 / 23.10.22 / --        handling of zero W, +/-inf slopes, discontinuous tables
%       v1.1 / 16.10.22 / --        epstol, extrapolation slopes for end bands
%       v1.0 / 14.10.22 / V.Yotov
%  ------------------------------------------------------------------------------------------------

classdef tab

properties
    val
    block
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
            otherwise error("tab: table must have 2 to 4 columns")
        end
        
        if ~issorted(T(:,1))
            h = msgbox(sprintf(['First column of input was not sorted. \nCorrect results ' ...
                'are not guaranteed in case of discontinuous tables.']), 'tab', 'warn');
            T = sortrows(T);
        end
    
    % Check / split / check / autofill / compute fields
        tab.validate(T, [-Inf Inf]);

        if opts.split
            T = tab.split(T);
        end
        
        tab.validate(T, [1:5]);
        T = tab.autofill(T, [0 -Inf Inf -1 10]);                                        % Zero slopes -> Inf slopes -> LS-RS pairs -> const. W
        T = tab.compute(T, opts);
    
    % Fill extrapolation slopes
        T = tab.autofill(T, 99, sl=[opts.ls opts.rs]);
    
    % Check completeness and clean-up
        tab.validate(T, [923, 1:5, -Inf, Inf]);
        if opts.collapse
            T = tab.collapse(T, tol=opts.epstol);
        end

    % Create object
        tb.val = T;
        tb.block = tab.getblocks(T, opts.epstol);
    end
end

methods(Static)
    T = split(T)
    T = autofill(T,flags,opts)   
    T = compute(T,tabopts)
    T = collapse(T,opts)

    mask = symeq(v,w,tol)
    ids = getblocks(T,tol)

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
tab(T); disp(ans.val)

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
tab(T); disp(ans.val)

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
tab(T,collapse=false); disp(ans.val)
tab(ans.val);  disp(ans.val)
%}
%  ------------------------------------------------------------------------------------------------










