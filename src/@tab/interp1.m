function Y = interp1(T,xq)

Y = nan(size(xq));
b = T.block;

% Sort query array
xq = rand(3,3,2);
[xs,idxQ2S] = sort(xq(:));
[~,idxS2Q] = sort(idxQ2S);

% test mapping
reshape(xs(idxS2Q),size(xq)) - xq

% Get faster intersections





for i = 1:size(b,1)
    if i==1
        mask = T.val(b(i,2)) <= xq;
    elseif i==size(b,1)
        mask = T.val(b(i,1)) >= xq;
    else
        mask = T.val(b(i,1)) <= xq & xq <= T.val(b(i,2)); % OK for mid bands...
    end
    Qi = T.val(b(i,1):b(i,2), :);
    Yi = interp1(log(Qi(:,1)), log(Qi(:,2)), log(xq(mask)), 'linear', 'extrap');
    Y(mask) = exp(Yi);
end


%{
% Add extrap bands if needed
    if min(f) < T(1,1)
        T = [min(f), nan(1,size(T,2)-1); T];        
    end
    if max(f) > T(end,1)
        T(end+1,:) = [max(f), nan(1,size(T,2)-1)];  
    end
    
% Recompute with new bands and interpolate
    T = tabc4(T);                                                                          
    W = interp1(log(T(:,1)), log(T(:,2)), log(f), 'linear');
    W = exp(W);


interp -> addpoints -> integral

%}
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
T = tab(T);
W = interp1(T,f);

[ax,fig]=xfig(n=1,xy=1,b=1);
    plot(T.val(:,1),T.val(:,2),'-k',linewidth=3.2,color=col('k'))
    plot(f,W,linewidth=1.2,color=col('coolgrey'),linestyle='--')

    ax.XLimitMethod = 'padded';
    legend("Table","Interpolated")


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
