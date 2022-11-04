function [C,ccmod] = modalDamping(M,K,zeta)
%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       [C,ccmod] = modalDamping(M,K,zeta)
%
%       See also:   cvec, eigsol, eigscale
%
%   OUTPUTS
%       C           structural damping matrix
%       ccmod       critical modal damping vector
%
%   INPUTS
%       K           stiffness matrix
%       M           mass matrix
%       zeta        scalar or vector of modal damping values
%
%   VERSION
%   v1.0 / 23.06.22 / V.Y.
%  ------------------------------------------------------------------------------------------------

% Make zeta n x 1 vec
switch length(zeta)                                     
    case 1     
        zeta = zeta*ones(length(M),1);
    otherwise  
        zeta = cvec(zeta,false);
end

% Get modal crit damping from 2*sqrt(k*m). 
% Note Q'*M*Q == I from eigsol
[Q,~] = eigsol(K,M); 
ccmod = 2*sqrt(diag(Q'*K*Q));

% Project to physical space
invQ = pinv(Q);
C = invQ'.*(zeta.*ccmod).'*invQ;








