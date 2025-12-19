function X = genMat(n,opts)
%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       X = genMat(n,opts)
%
%       See also: genPSD
%
%   INPUTS
%       n                   matrix size
%       opts{:}
%           type
%               "integer"   diagonal, integer
%               "real"      diagonal, real 
%               "chain-K"   real chain, free-free BC
%               < else >    real, PSD
%           range           interval of diag entries
%           eigRange        interval of eigenvalues
%           sparse          output format
%
%   VERSION 
%   v2.0 / xx.xx.xx / --    [-] int/real/complex, [-] symm/skew/hermitian
%   v1.1 / xx.12.25 / --    opts struct / eigRange default [0,1] / nonsparse default
%                           [-] drange, mrange, [-] sparse, [-] genPSD formalisation 
%   v1.0 / 16.12.25 / V.Y.  
%  ------------------------------------------------------------------------------------------------

    arguments 
        n (1,1) double {mustBeInteger,mustBePositive} = 1
        opts.type (1,1) string = "psd"
        opts.range (1,2) double = [0,n]
        opts.eigRange (1,2) double = [0,1]
        opts.sparse (1,1) = false
    end

    switch opts.type
        case "integer"                                                                              % diagonal / rand integer in dRange
            X = randi(opts.range,n,1);                                                                  
        case "real"                                                                                 % diagonal / rand real in dRange
            X = opts.range(1) + rand(n,1)*diff(opts.range);                                                 
        case "chain-K"                                                                              % chain / free-free BC
            tmp = randi(opts.range,n-1,1);
            dsub = -[tmp(1:end); 0];
            dsup = -[0; tmp(1:end)];
            X = [dsub, -dsub-dsup, dsup];
            X = spdiags(X,-1:1,n,n);
        otherwise
            X = genPSD(n, opts.eigRange);                                                           % dense
    end

    if ismember(opts.type, ["integer" "real"])
        if opts.sparse
            X = spdiags(X,0,n,n);                                                                   % sparse diag from vector X
        else
            X = diag(X);
        end
    end
 
%  ------------------------------------------------------------------------------------------------
%{
% [EXAMPLE 1]
    n = 8;
    xi = genMat(n,type="integer")
    xr = genMat(n,type="real",sparse=true)
    xc = genMat(n,type="chain-K")
    xp = genMat(n)
%}
%  ------------------------------------------------------------------------------------------------