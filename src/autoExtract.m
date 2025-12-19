function autoExtract(x,sys,dof,opts)
%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       autoExtract(x,sys,dof,opts)
%
%       Computes entris that can be unambiguously extractd sbm but are not explicitly present 
%       in sys; May overwrite sys, dof in caller function
%
%       See also:       sbm, submat, eigsol, s2vars
%
%   INPUTS
%       x               matrix, e.g. 'TAG' or "T.A.Gb"
%       sys             struct with system matrices
%       dof             struct with DOF sets
%       opts{:}
%           updateSys   update sys/dof in caller, default=false
%           autofill    compute P/G/modal, default=false
%
%   VERSION
%   v1.1 / xx.xx.xx / --    [-] examples / [-] damping autofill from c/zeta/C
%   v1.0 / 16.12.25 / V.Y.  recursively computes G/P/X as needed
%  ------------------------------------------------------------------------------------------------

    arguments
        x {mustBeText}
        sys
        dof
        opts.updateSys = false
        opts.autofill = false
    end
    s2vars(opts)

    % flag
    updateSys = autofill || updateSys;                                                                  % ensures correct assignin later

    % autofill G/P/modal
    if autofill && x(1)~='T'
        autoExtract('T',sys,dof,updateSys=true)                                                         % autofill=false prevents infinite recursion
    end
    
    % recursive eval
    switch x(1)
        case 'G'
            if ~isfield('G','sys')
                KAA = submat(sys.K,dof.G,dof.A);
                KAS = submat(sys.K,dof.G,dof.A,dof.G,dof.S);
                sys.G = -pinv(full(KAA))*KAS;
            end

        case 'P'
            if ~isfield('P','sys')
                if ~isfield('G','sys')
                    autoExtract('G',sys,dof,updateSys=true);
                end
                MAA = submat(sys.M,dof.G,dof.A);
                MAS = submat(sys.M,dof.G,dof.A,dof.G,dof.S);
                sys.P = -(MAA*sys.G + MAS);  
            end

        case 'X'
            if ~isfield('X','sys')
                [sys.X,sys.k] = eigsol(sbm('KAA',sys,dof),sbm('MAA',sys,dof));
                sys.k = diag(sys.k);

                mask = ~isinf(sys.k);
                sys.X = sys.X(:,mask);
                sys.k = sys.k(mask);
                sys.zeta = sys.zeta(mask);

                sys.c = 2*sys.zeta.*sqrt(sys.k);
                sys.m = ones(size(sys.k)); 
                dof.m = (1:numel(sys.m))';
            end

        case 'H'
            if ~isfield('X','sys')
                autoExtract('X',sys,dof,updateSys=true);
            end

        case 'T'
            if ~isfield('X','sys')
                autoExtract('X',sys,dof,updateSys=true);                                                % fills modal properties
            end
            if ~isfield('P','sys')
                autoExtract('P',sys,dof,updateSys=true);                                                % auto-checks for G 
            end

        otherwise % do nothing
    end % switch
    
    % update in caller
    if updateSys
        assignin('caller','dof',dof)
        assignin('caller','sys',sys)
    end