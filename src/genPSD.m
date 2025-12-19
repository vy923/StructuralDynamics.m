function X = genPSD(n,eigRange)
%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       X = genPSD(n,eigRange)
%
%   INPUTS
%       n           matrix size
%       eigRange    default = [0,1], uniform i.i.d.
%
%   VERSION
%   v1.0 / 16.12.25 / V.Y.  
%  ------------------------------------------------------------------------------------------------

    arguments 
        n (1,1) double {mustBeInteger} = 1
        eigRange (1,2) double = [0,1]
    end

    lam = eigRange(1) + rand(n,1)*diff(eigRange);                                                   % eigs uniform i.i.d
    [Q,~] = qr(randn(n));                                                                           % SO(n) with Haar measure
    X = Q.*lam.'*Q';

%  ------------------------------------------------------------------------------------------------
%{
% [EXAMPLE 1]
    n = 12;
    x = genPSD(n); eig(x)
    x = genPSD(n,[10 1e3]); eig(x)
    x = genPSD(n,[-10 5]); eig(x)
%}
%  ------------------------------------------------------------------------------------------------