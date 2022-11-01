function T = compute(T,opts)
%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       T = tab.compute(T,opts)
%
%       See also:       autofill, isapprox
%       Related:        collapse, split, validate
%
%   INPUTS
%       T               [n x 4] array in the form of tab.val
%       opts{:}         directly passed from tab class constructor
%
%   OUTPUTS
%       T               [n x 4] with filled and consistency-validated values
%
%   VERSION
%       v1.0 / 26.10.22 / V.Yotov
%  ------------------------------------------------------------------------------------------------

X1 = T(1:end-1,1);
X2 = T(2:end,1);
logX = log10(X2./X1)/log10(2);

% Chained calc/validate until no further values can be computed
maskBlock = true(size(T,1)-1, 1);

while any(maskBlock)
    
    LS = T(2:end,3); 
    RS = T(1:end-1,4); 
    Y1 = T(1:end-1,2);
    Y2 = T(2:end,2);

    % Exclude singularities from explicit computation
    maskBlock( isinf(RS) | isinf(LS) | Y1.*Y2==0 ) = false;                            

    % Fields computable from current data, may overlap
    maskR = ~isnan(Y1 + RS) & maskBlock;                                         
    maskL = ~isnan(LS + Y2) & maskBlock;
    maskY = ~isnan(Y1 + Y2) & maskBlock;

    % Prevent infinite loops
    assert(any( maskR | maskL | maskY ));                                             

    idxR = find(maskR);
    idxL = find(maskL);
    idxY = find(maskY);

    % Compute, incl. with overlapping methods
    tmpY2 = Y1(idxR) .* 10.^(RS(idxR).*logX(idxR)/10);                          % <--- allow for other computational rules in next versions
    tmpY1 = Y2(idxL) .* 10.^(RS(idxL).*-logX(idxL)/10);
    tmpRS = 10*log10(Y2(idxY)./Y1(idxY)) ./ logX(idxY);

    % Verify that different comp. methods give consistent results
    assert( all([
            isapprox(Y2(maskR & ~isnan(Y2.*maskR)), tmpY2(~isnan(Y2(maskR))), opts.atol)
            isapprox(Y1(maskL & ~isnan(Y1.*maskL)), tmpY1(~isnan(Y1(maskL))), opts.atol)
            isapprox(RS(maskY & ~isnan(RS.*maskY)), tmpRS(~isnan(RS(maskY))), opts.atol)
            ]), ...
        "tab: overspecified input");

    % Do not overwrite existing non-NaN fields
    maskRS = isnan(T(idxR+1,2));
    maskLS = isnan(T(idxL,2));
    maskYR = isnan(T(idxY,4)); 
    maskYL = isnan(T(idxY+1,3)); 

    T(idxR(maskRS)+1,2) = tmpY2(maskRS);
    T(idxL(maskLS),2) = tmpY1(maskLS);
    T(idxY(maskYR),4) = tmpRS(maskYR);
    T(idxY(maskYL)+1,3) = -tmpRS(maskYL);

    % "aggressive" autofill for +/-Inf slopes
    T = tab.autofill(T, [Inf -1]);  

    % Discard checked entries from next iteration
    maskBlock( maskR | maskL | maskY ) = false;
end





