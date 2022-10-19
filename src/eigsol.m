%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       [x,lam] = eigsol(K,varargin)
%       Wrapper for built-in eig() that forces eigenvector scaling   
%
%   VERSION
%       v1.1 / 16.12.21 / V.Yotov
%  ------------------------------------------------------------------------------------------------

function [x,lam] = eigsol(K,varargin)

    n = size(K,1);                                                          % No verification for square matrices

    if nargin == 1          [x,lam] = eig(K);                               % Switch between ordinary/full eigenv. problem
    elseif nargin == 2      [x,lam] = eig(K,varargin{1});                   % Same as eig(K,M)
    else                    [x,lam] = eigs(K,varargin{:});                  % Call with (K,M,nEigs,'opts') or (K,nEigs,'opts')
    end

    [x,lam] = eigscale(K,x,lam);
