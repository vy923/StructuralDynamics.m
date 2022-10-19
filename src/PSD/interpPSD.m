%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       [W,T] = interpPSD(T,f)
%       Interpolates a PSD table given by T over frequency vector f 
%
%       See also:       tablePSD
%       Related:        scalePSD, integratePSD
%
%   VERSION
%       v1.1 / 16.10.22 / V Yotov   extrapolation uses actual end band slopes
%       v1.0 / 14.10.22 / V.Yotov
%  ------------------------------------------------------------------------------------------------

function [W,T] = interpPSD(T,f)

% Validate/complete T
    T = tablePSD(T);
    
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
