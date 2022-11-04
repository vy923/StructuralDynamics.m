function [x,lam,xscale] = eigscale(K,x,lam)
%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       [x,lam,xscale] = eigscale(K,x,lam)
%       Reordering and scaling of eigensolutions
%
%       See also: eigsol
%
%   VERSION
%   v1.1 / 16.12.21 / V.Y.
%  ------------------------------------------------------------------------------------------------

q = size(x,2);

[tmpLam,ind] = sort(diag(lam));                                         % Sort the eigenvalues and find sorting index
lam(1:q+1:end) = tmpLam;                                                % Reorder eigenvalue matrix: faster than calling diag()
x = x(:,ind);                                                           % Reorder eigenvector matrix

xscale = (sum((K'*x).*x)'./tmpLam).^-0.5;                               % Eigenvector scaling factor / 2x faster than diag(x'*K*x) / works for truncated eigenbasis too
x = x.*xscale';                                                         % Scale to diagonalise M to I and K to lambda    

