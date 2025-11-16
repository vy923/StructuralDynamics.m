%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       [Y,T] = interpc4(T,f,extrap)
%       Interpolates a table given by T
%
%       See also:       tab
%       Related:        scalec4, integratec4
%
%   VERSION
%       v1.3 / 10.05.24 / --    bugfixes for extrap = true
%       v1.2 / 21.03.24 / --    size(W) == size(f) / tabc4 -> tab / efficiency improvements
%       v1.1 / 16.10.22 / --    extrapolation uses actual end band slopes
%       v1.0 / 14.10.22 / V.Y.
%  ------------------------------------------------------------------------------------------------

function [Y,T] = interpc4(T,f,extrap)

% Require outside original range
    if nargin<3 || isempty(extrap)
        extrap = true;
    end

% Validate/complete T
    if ~isa(T,'tab')
        T = tab(T);
    end
    V = T.val;
    Y = zeros(size(f));

% Add extrap bands if needed
    minf = min(f,[],'all');
    maxf = max(f,[],'all');

    if extrap
        if minf < T.val(1,1)
            V = [minf, nan(1,size(V,2)-1); V];        
        end
        if maxf > V(end,1)
            V(end+1,:) = [maxf, nan(1,size(V,2)-1)];  
        end
        if size(V,1) ~= size(T.val,1)                                                              % Update T if necessary
            T = tab(V,collapse=false);
            V = T.val;
        end
        mask = true(size(f));
    else
        mask = V(1,1)<=f & V(end,1)<=maxf;
    end
    
% Recompute with new bands and interpolate
    Y(mask) = interp1(log(V(:,1)),log(V(:,2)),log(f(mask)),'linear');
    Y(mask) = exp(Y(mask));

% -----
% split table in parts, get query poit intersections with subtable bands,
% special cases being the ends only

%  ------------------------------------------------------------------------------------------------
%{
% Example 1
T =  [  20          0.026
        50          0.06
        150         0.06
        300         0.02
        700         0.02 
        800         0.06
        925         0.06
        2000        0.02    ];

T = [     20	    nan	        nan          +6	
          50        nan         nan          0
         250        .05         nan          nan
         320        nan         nan          3
         300        nan         nan          -4
         350        nan         nan          3
         380        nan         nan          3
         400        0.100       nan          nan
         500        nan         -1           nan
         550        nan         -1           nan
         590        nan         -1           nan
         800	    nan	        nan          nan
        2000        0.026       +6	         nan     ];

f = [1:1:6200]';
T = tabc4(T);
W = interpc4(T,f);

[ax,fig]=xfig(n=1,xy=1,b=1);
    plot(T(:,1),T(:,2),'-k',linewidth=3.2,color=col('k'))
    plot(f,W,linewidth=1.2,color=col('coolgrey'),linestyle='--')

    ax.XLimitMethod = 'padded';
    % ax.YLimitMethod = 'padded';
    legend("Table","Interpolated")

% Example 2
T = [     20	        nan	        0            +6	
          50            nan         nan          0
          800	        nan	        nan          nan
          2000        0.026         +6	         nan        ];

f = [4:.1:6200]';
T = tabc4(T,ls=-15,rs=-3);
[W,Te] = interpc4(T,f);

[ax,fig]=xfig(n=1,xy=1,b=1);
    plot(T(:,1),T(:,2),'-k',linewidth=3.2,color=col('k'))
    plot(f,W,linewidth=1.2,color=col('coolgrey'),linestyle='--')

    ax.XLimitMethod = 'padded';
    % ax.YLimitMethod = 'padded';
    legend("Table","Interpolated")
%}
