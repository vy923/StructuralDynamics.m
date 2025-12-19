function [x,lam] = eigsol(K,varargin)
%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       [x,lam] = eigsol(K,varargin)
%       Wrapper for built-in eig() that forces eigenvector scaling  
%
%   VERSION
%   v1.3 / 16.12.25 / --    auto handling of sparse inputs
%   v1.2 / 30.10.22 / --    clean-up
%   v1.1 / 16.12.21 / V.Y.
%  ------------------------------------------------------------------------------------------------

flag = issparse(K) || (nargin > 1 && issparse(varargin{1}));

if flag && nargin <= 2 
    K = full(K);                                        % k < n mandatory for eigs()  
    if nargin == 2
        varargin{1} = full(varargin{1});
    end
end

switch nargin
    case 1
        [x,lam] = eig(K);                               % Switch between ordinary/full eigenv. problem
    case 2 
        [x,lam] = eig(K,varargin{1});                   % Same as eig(K,M)
    otherwise
        [x,lam] = eigs(K,varargin{:});                  % Call with (K,M,nEigs,'opts') or (K,nEigs,'opts')
end

[x,lam] = eigscale(K,x,lam);
