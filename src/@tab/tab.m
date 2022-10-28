%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       tb = tab(T,opts)
%
%       Constructor for the tab class
%       Objects hold a piecewise linear representation of a function in the form [X Y LS RS], where 
%       LS and RS are the left and right slopes at points [X Y]
%
%       See also:       mapSet, mustBeMemberSCI
%       Static:         split, autofill, compute, collapse, getblocks, validate, symeq
%       Overloads:      interp1
%
%   INPUTS
%       T               [n+1 x ncols] array with 2, 3 or 4 columns and NaN empty fields
%           2 columns   assumes [X Y]
%           3 columns   assumes [X Y RS]
%           4 columns   must be [X Y LS RS]
%
%       opts{:}
%           type        data type, currently only 'PSD'
%               'PSD'   default, assumes log-log and dB/oct slope units
%
%           split       autosplit table on discontinuities
%               true    default, handles +/-inf, discontinuities 
%               false   assumes double points are given explicitly
%
%           collapse    default = true
%               true    removes redundant points
%               false   allows redundant points, e.g. midpoint of constant band
%
%           epstol      eps tolerance scaling for comparisons, default = 1E+4  
%           ls          left slope of first block
%           rs          right slope of final block
%
%   OUTPUTS
%       tb              initialised object
%
%   UPDATES
%       - merge function with min/max options for adding profiles, etc.
%       - inerp1 -> addpoints -> integral
%
%   VERSION
%       v2.0 / 26.10.22 / --        redefined as 'tab' class, updated methods for the constructor
%       v1.2 / 23.10.22 / --        handling of discontinuities i.e. zero Y, +/-Inf slopes 
%       v1.1 / 16.10.22 / --        epstol, extrapolation slopes for end blocks
%       v1.0 / 14.10.22 / V.Yotov
%  ------------------------------------------------------------------------------------------------

classdef tab

properties
    val
    block
end

methods    
    % CONSTRUCTOR
    function tb = tab(T,opts) 
        arguments
            T
            opts.type {mustBeMemberSCI(opts.type,"PSD")} = "PSD"
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
            otherwise; error("tab: table must have 2 to 4 columns")
        end
        
        if ~issorted(T(:,1))
            warning(['tab: first column of input was not sorted. Correct results ' ...
                'are not guaranteed in case of discontinuous tables.'], 'tab', 'warn');
            T = sortrows(T);
        end
        

        % Check / split / check / autofill / compute fields
        tab.validate(T, [-Inf Inf]);
        if opts.split
            T = tab.split(T);
        end
        tab.validate(T, [1:5]);
        T = tab.autofill(T, [0 -Inf Inf -1 10]);                                        % Zero slopes -> Inf slopes -> LS-RS pairs -> const. Y
        T = tab.compute(T, opts);
    
        % fill extrapolation slopes
        T = tab.autofill(T, 99, sl=[opts.ls opts.rs]);
    
        % check completeness and clean-up
        tab.validate(T, [923, 1:5, -Inf, Inf]);
        if opts.collapse
            T = tab.collapse(T, tol=opts.epstol);
        end

        % create object
        tb.val = T;
        tb.block = tab.getblocks(T, opts.epstol);
    end

    % OVERLOADS
    Y = interp1(T,f,opts)
end

methods(Static)
    % PROPERTIES
    T = split(T)
    T = autofill(T,flags,opts)   
    T = compute(T,opts)
    T = collapse(T,opts)
    ids = getblocks(T,tol)
    validate(T,flags)   
    
    % UTILITY
    mask = symeq(v,w,tol)
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
T = [     32          15        -Inf         nan
          40          94         nan        -Inf
          50          69         nan         nan
          55          77         nan         nan
          63           0         0           nan
          80          86        -Inf         nan
          80           0         nan         nan
         100           0         nan         nan
         125           0         0           Inf
         125          10         nan         nan
         160         108         nan         nan
         400           0         nan         nan
        1000           4         nan         nan
        1600           1        -Inf         nan
        2000           1         nan         nan
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
         400        0.1         nan          nan
         500        nan         -1           nan
         550        nan         -1           nan
         590        nan         -1           nan
         800	    nan	        nan          nan
        2000        0.026       +6	         nan     ];
tab(T,collapse=false); disp(ans.val)
tab(ans.val);  disp(ans.val)
%}
%  ------------------------------------------------------------------------------------------------










