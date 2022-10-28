%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       [B,idxO2N,idxN2O] = mapSet(SO,SN,A,opts)
%
%       See also:       mustBeOfSize, mustBeEqualDims
%       Related:        SELASM
%
%   INPUTS
%       SO              set to be mapped
%       SN              target set
%       A               array with at least one dim corresponding to SO
%       opts{:}             
%           dims        dims of A to be mapped to SN, default = 1
%
%   OUTPUTS
%       B               remapped array of same number of dims as A
%       idxO2N          forward  map idxs SO->SN: SO(idxO2N) == SN(sort(idxN2O))
%       idxN2O          backward map idxs SN->SO: SN(idxN2O) == SO(sort(idxO2N))
%
%   UPDATES
%       - do operations on pre-initialised B, replace/add/multiply with A
%
%   VERSION
%       v1.0 / 22.10.22 / V.Yotov
%  ------------------------------------------------------------------------------------------------

function [B,idxO2N,idxN2O] = mapSet(SO,SN,A,B,opts)

arguments
    SO {mustBeVector}
    SN {mustBeVector}
    A = []
    B = []
    opts.dims {mustBeInteger,mustBePositive} = 1
end
if ~isempty(A)
    mustBeOfSize(A,length(SO),opts.dims)
    flag = true;
end

% Intersection indices
    [~,~,idxN2O] = intersect(SO,SN,'stable');
    [~,~,idxO2N] = intersect(SN,SO,'stable');
    
% Remap A if nonempty
if flag
    sz = size(A);
    sz(opts.dims) = length(SN);
    
    cla = repmat({':'},size(sz));
    clb = cla;
    cla(opts.dims) = {idxO2N};
    clb(opts.dims) = {sort(idxN2O)};

    B = zeros(sz); 
    B(clb{:}) = A(cla{:});
end


%  ------------------------------------------------------------------------------------------------
%{
% Example 1
    A = zeros(4); 
    A(1:numel(A)) = 10*[1:numel(A)];
    SO = [9 3 7 4]';
    SN = [8:10 4:-1:1 6]';
    [B,SOtoSN,SNtoSO] = mapSet(SO,SN,A,dims=[1 2]);
        disp(A)
        disp(B)
        disp(SO')
        disp(SN')

%  ------------------------------------------------------------------------------------------------

% [OBSOLETE CODE] Mapping masks 
    if opts.type=='mask'
        msg = "mapSet: required 'mask' output but input sets are not sorted";
        assert(issorted(SO) & issorted(SN), msg);

        maskSOtoSN = false(size(SO));
        maskSOtoSN(SOtoSN) = SO(SOtoSN)~=0;
        maskSNtoSO = false(size(SN));
	    maskSNtoSO(SNtoSO) = SO(SOtoSN)~=0;
    end

%  ------------------------------------------------------------------------------------------------
%}




















