function [x,lam] = eigsol(K,varargin)
%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       [x,lam] = eigsol(K,varargin)
%       Wrapper for built-in eig() that forces eigenvector scaling  
%
%   VERSION
%       v1.2 / 30.10.22 --          clean-up
%       v1.1 / 16.12.21 / V.Yotov
%  ------------------------------------------------------------------------------------------------

switch nargin
    case 1
        [x,lam] = eig(K);                               % Switch between ordinary/full eigenv. problem
    case 2 
        [x,lam] = eig(K,varargin{1});                   % Same as eig(K,M)
    otherwise
        [x,lam] = eigs(K,varargin{:});                  % Call with (K,M,nEigs,'opts') or (K,nEigs,'opts')
end

[x,lam] = eigscale(K,x,lam);
