%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       [W,T] = interpc4(T,f)
%       Interpolates a c4 table given by T
%
%       See also:       tabc4
%       Related:        scalec4, integratec4
%
%   UPDATES
%       - interpolate accross a table split by discontinuity
%       - handle +/-Inf slopes consistently with mid-band discontinuity
%
%   VERSION
%       v2.0 / 22.10.22 / --        interpolation of discontinuous tables
%       v1.1 / 16.10.22 / --        extrapolation uses actual end band slopes
%       v1.0 / 14.10.22 / V.Yotov
%  ------------------------------------------------------------------------------------------------

function [W,T] = interpc4(T,f)

% Validate/complete T
    T = tabc4(T);
    
% Add extrap bands if needed
    if min(f) < T(1,1)
        T = [min(f), nan(1,size(T,2)-1); T];        
    end
    if max(f) > T(end,1)
        T(end+1,:) = [max(f), nan(1,size(T,2)-1)];  
    end
    
% Recompute with new bands and interpolate
    T = tablePSD(T);                                                                          
    W = interp1(log(T(:,1)), log(T(:,2)), log(f), 'linear');
    W = exp(W);

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

f = [1:1:6200]';
W = interpPSD(T,f);

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
T = tablePSD(T,ls=-15,rs=-3);
[W,Te] = interpPSD(T,f);

[ax,fig]=xfig(n=1,xy=1,b=1);
    plot(T(:,1),T(:,2),'-k',linewidth=3.2,color=col('k'))
    plot(f,W,linewidth=1.2,color=col('coolgrey'),linestyle='--')

    ax.XLimitMethod = 'padded';
    % ax.YLimitMethod = 'padded';
    legend("Table","Interpolated")
%}
