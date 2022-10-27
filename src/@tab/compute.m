function T = compute(T,opts)                                                    % opts is directly passed from tab()

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
        "tab: overspecified input");

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