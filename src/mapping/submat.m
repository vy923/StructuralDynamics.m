function A = submat(A,SO,SN)
%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       A = submat(A,SO,SN)
%
%       See also:       mapSet, mustBeVector
%
%   INPUTS
%       A               SO1 x SO2 x ... x SOn 
%       SO/SN           repeating old/new DOF set pairs
%
%   OUTPUTS
%       A               SN1 x SN2 x ... x SNn
%
%   VERSION
%   v1.0 / 23.02.24 / V.Y.
%  ------------------------------------------------------------------------------------------------

arguments
    A
end
arguments(Repeating)
    SO {mustBeVector}
    SN {mustBeVector}
end

% Validate input sizes
    na = ndims(A);
    ns = numel(SO);
    flag = ns==[1 na];
    assert(any(flag),"submat: SO, SN pairs must be 1 or ndims(A)")

% Remap array
    for d = 1:ns                                                                                    % d = 1 / 1:na
        A = mapSet(SO{d},SN{d},A,dims=d:max(flag.*[na d]));                                         % dims = 1:na / d
    end

%  ------------------------------------------------------------------------------------------------
%{
% Example 1: a wrapper for mapSet

    A = zeros(4); 
    A(1:numel(A)) = 10*[1:numel(A)];
    SO = [9 3 7 4]';
    SN = [8:10 4:-1:1 6]';

    [B,SOtoSN,SNtoSO] = mapSet(SO,SN,A,dims=[1 2]);
        disp(A)
        disp(B)
        disp(SO')
        disp(SN')

    disp(submat(A,SO,SN))
    disp(submat(A,SO,SN,SO,SN))
%}
%  ------------------------------------------------------------------------------------------------







