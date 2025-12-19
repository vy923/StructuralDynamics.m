function X = sbm(x,sys,dof,opts)
%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       X = SBM(x,sys,dof,opts)
%
%       See also:   submat, autoExtract, s2vars
%
%   INPUTS
%       x           matrix, e.g. 'TAG' or "T.A.Gb"
%       sys         struct with system matrices
%       dof         struct with DOF sets
%       opts{:}
%           w       [rad/s] discrete frequencies omega 
%           fmax    [Hz] 
%
%   OUTPUTS
%       X           submatrix/FRF from sys
%
%   VERSION
%   v1.2 / xx.xx.xx / --    [-] examples
%   v1.1 / 16.12.25 / --    autoExtract recursively computes G/P/X as needed
%   v1.0 / 25.10.25 / V.Y.  from dummySC v1.0
%  ------------------------------------------------------------------------------------------------

    arguments
        x {mustBeText}
        sys
        dof
        opts.w {mustBeNumeric} = []
        opts.fmax {mustBeScalarOrEmpty} = []
    end
    s2vars(opts)

    x = strsplit(string(x),".");
    if numel(x)~=3
        x = char(x);
    end

    % autocomplete
    autoExtract(x,sys,dof,updateSys=true);

    % output
    if any(x(1)==["G" "P"])                                                                             % Default DOF sets of P, G are A x S
        X = submat(sys.(x(1)),dof.A,dof.(x(2)),dof.S,dof.(x(3)));

    elseif any(x(1)==["M" "C" "K"])                                                                     % Default DOF sets of M, C, K are G x G
        X = submat(sys.(x(1)),dof.G,dof.(x(2)),dof.G,dof.(x(3)));

    elseif x(1)=="X"                                                                                    % Default DOF set of X is A x m
        X = submat(sys.(x(1)),dof.A,dof.(x(2)),dof.m,dof.(x(3)));

    elseif any(x(1)==["H" "T"]) && ~isempty(w)
        if ~isempty(fmax)
            assert(isvector(w), 'sbm: output freqs w must be a vector when fmax is specified')
            q = sys.k<=(2*pi*fmax)^2;
            sys.m = sys.m(q);
            sys.c = sys.c(q);
            sys.k = sys.k(q);
            idx = {dof.m,find(q)};
        else
            idx = {[],[]};
        end

        w = reshape(w,[1 1 size(w)]);
        X1q = submat(sys.X,dof.A,dof.(x(2)),idx{:});
      
        switch x(1)
            case 'H'
                X2q = submat(sys.X,dof.A,dof.(x(3)),idx{:});
                X = pagemtimes(X1q./pagetranspose(-w.^2.*sys.m + sqrt(-1)*w.*sys.c + sys.k), X2q');
            case 'T'
                XAq = submat(sys.X,[],[],idx{:});
                PA2 = submat(sys.P,[],[],dof.S,dof.(x(3)));
                G12 = submat(sys.G,dof.A,dof.(x(2)),dof.S,dof.(x(3)));
                X = pagemtimes(X1q./pagetranspose(sys.m - sqrt(-1)./w.*sys.c - sys.k./w.^2), XAq'*PA2) + G12;
        end

    else
        error(['Unsupported requested matrix ' x])
    end

