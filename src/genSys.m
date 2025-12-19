function sys = genSys(dof,opts)
%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       X = genSys(dof,opts)
%
%       See also:           genMat, genPSD, autoExtract, modalDamping, s2vars
%
%   INPUTS
%       dof
%       opts{:}
%           sysType         
%               "chain"     free-free BC 
%               "dense-K"   diag M, full K
%               "dense"     full M, full K
%           massType        
%               "integer"   int diag
%               "real"      real diag
%               "full"      dense
%           massLim         range of diag(M)
%           zetaLim         range of modal damping ratios
%           springLim       range of diag(K) entries
%           nf              nominal NFs [Hz] 
%
%   VERSION
%   v1.1 / xx.12.25 / --    [-] default values, [-] distributions, [-] examples
%   v1.0 / 16.12.25 / V.Y.  
%  ------------------------------------------------------------------------------------------------

    arguments
        dof
        opts.sysType {mustBeMemberSCI(opts.sysType,["chain","dense-K","dense"])} = "chain"
        opts.massType {mustBeMemberSCI(opts.massType,["integer","real","full"])} = "integer"
        opts.massLim (1,2) double = [1 20]                                                          % mass range per DOF
        opts.zetaLim (1,2) double = [1e-3 1e-1]                                                     % modal damping range
        opts.springLim (1,2) double = [1e2 1e4]                                                     % diag K range
        opts.nf double {mustBeVector} = 1:numel(dof.G)                                              % default NFs 1,2,...n Hz / ignored for "chain"
    end
    s2vars(opts)

    % checks
    assert( ~(sysType=="dense-K" && massType=="full"), ...
            "genSys: incompatible inputs" )

    % generate M
    s = numel(dof.S);
    a = numel(dof.A);
    g = numel(dof.G);
    M = genMat(g,type=massType,range=massLim);

    % generate K
    switch sysType
        case "chain"
            K = genMat(g,type="chain-K",range=springLim);
        case "dense-K"
            lam = (cvec(nf)*2*pi).^2;                                                               % eigenvalues
            Lm = chol(M,'lower');
            [Qh,~] = qr(randn(g));                                                                  % random eigenbasis for H
            H = Qh.*lam'*Qh';                                                                       % generalised -> ordinary eigenproblem matrix
            K = Lm'*H*Lm;
        otherwise
            K = genPSD(g,[1 100]);
    end

    % base system
    sys.M = M;
    sys.K = K;
    sys.C = double.empty;
    sys.zeta = min(zetaLim) + diff(zetaLim)*rand(a,1);                                              % random A set modal damping

    % G/P/modal completion
    autoExtract('T',sys,dof,autofill=true);

    % A set matrices
    MAA = submat(M,dof.G,dof.A);
    KAA = submat(K,dof.G,dof.A);                                        
    CAA = modalDamping(MAA,KAA,sys.zeta);                                                           % modal --> physical

    % RBM-consistent expansion
    Y = [eye(a) zeros(a,s); -sys.G' eye(s)];
    Y = submat(Y, [dof.A; dof.S], dof.G);                                                           % reorder to actual G set
    sys.C = Y*submat(CAA,dof.A,dof.G)*Y';                                                           % verify / full(K)-Y*submat(KAA,dof.A,dof.G)*Y' / should be approx 0 


