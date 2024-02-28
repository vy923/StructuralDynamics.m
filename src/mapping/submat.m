function A = submat(A,SO,SN)
%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       A = submat(A,SO,SN)
%
%       See also:       mapSet, mustBeVector
%
%   INPUTS
%       A               SO1 x SO2 x ... x SOn 
%       SO/SN           repeating old/new DOF set pairs or empty
%
%   OUTPUTS
%       A               SN1 x SN2 x ... x SNn
%
%   VERSION
%   v1.1 / 24.02.24 / --    single SO, SN for column vector A / skip empty DOF sets / examples
%   v1.0 / 23.02.24 / V.Y.
%  ------------------------------------------------------------------------------------------------

arguments
    A
end
arguments(Repeating)
    SO {mustBeVector(SO,'allow-all-empties')}
    SN {mustBeVector(SN,'allow-all-empties')}
end

% Argument sizes
    if isvector(A) && numel(SO)==1
        SO{2} = [];
        SN{2} = [];
    end  
    na = ndims(A);
    ns = numel(SO);
    flag = ns==[1 na];

% Break
    assert(any(flag),"submat: SO, SN pairs must be 1 or ndims(A)")

% Skip empty SO sets on remap
    for d = find(~cellfun(@isempty,SO))                                                             % d = 1 / 1:na
        A = mapSet(SO{d},SN{d},A,dims=d:max(flag.*[na d]));                                         % dims = 1:na / d
    end


%  ------------------------------------------------------------------------------------------------
%{
% Example 1: a wrapper for mapSet
    A = zeros(4); 
    A(1:numel(A)) = 10*[1:numel(A)];
    SO = [9 3 7 4]';
    SN = [8:10 4:-1:1 6]';

    [B,O2N,N2O] = mapSet(SO,SN,A,dims=[1 2]);

    disp(A)
    disp(B)
    disp(SO')
    disp(SN')

    disp(submat(A,SO,SN))
    disp(submat(A,SO,SN,SO,SN))

% Example 2: skipping dof sets
    A = reshape(1:24,[3 2 4])
    submat(A,[],[], 3:4,3:-1:1, 1:4,2:3)
   
% Example 3: vector inputs
    A = cvec(10:-1:1);
    SO = 2:11;
    SN = 15:-2:-5;

    %   Column vec
    [   submat(A,SO,SN), ... 
        submat(A,SO,SN,[],[])       ]  

    %   Row vec
    [   submat(A',[],[],SO,SN)  
        submat(A',1,1,SO,SN)   
        submat(A,SO,SN)'            ] 

    %   Dims errors
        submat(A,[],[],SO,SN)                                                                       
        submat(A,SO,SN,[],[],[],[])
        submat(A.',SO,SN,[],[])
%}
%  ------------------------------------------------------------------------------------------------







