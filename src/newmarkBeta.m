function [x,xd,xdd] = newmarkBeta(M,C,K,f,dt,beta,gamma,x0)
%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       [x,xd,xdd] = newmarkBeta(M,C,K,f,dt,beta,gamma,x0)
%
%       See also:       --
%       Related:        SRS, forcedAcce1D
%
%   OUTPUTS
%       x/xd/xdd        disp/vel/acce [dof x timesteps] 
%
%   INPUTS
%       dt              time step
%       beta            default = 0.25
%       gamma           default = 0.5
%       x0              initial displacement, e.g. 0
%
%   VERSION
%   v1.1 / 23.06.22 / V.Y.
%  ------------------------------------------------------------------------------------------------

% Constants
a0 = 1/(beta*dt^2);
a1 = gamma/(beta*dt);
a2 = 1/(beta*dt);
a3 = 1/(2*beta) - 1;
a4 = gamma/beta - 1;
a5 = (dt/2)*(gamma/beta-2);
a6 = dt*(1-gamma);
a7 = gamma*dt;

Khat = K + a0*M + a1*C;

% Preallocation
N   = length(M);
nt  = size(f,2);
x   = zeros(N,nt);
xd  = zeros(N,nt);
xdd = zeros(N,nt);

% Initialise solution
x(:,1) = x0;

% Integration
if isvector(K) 
    % diagonal M, C, K given as vecotrs
    for i = 1:nt-1
        fhat        = f(:,i+1) + ...
                      M.*(a0*x(:,i) + a2*xd(:,i) + a3*xdd(:,i)) + ...
                      C.*(a1*x(:,i) + a4*xd(:,i) + a5*xdd(:,i));
        x(:,i+1)    = Khat.\fhat;
        xdd(:,i+1)  = a0*(x(:,i+1)-x(:,i)) - a2*(xd(:,i)) - a3*(xdd(:,i));
        xd(:,i+1)   = xd(:,i) + a6*xdd(:,i) + a7*xdd(:,i+1);
    end
else 
    % general M, C, K matrices
    for i = 1:nt-1
        fhat        = f(:,i+1) + ...
                      M*(a0*x(:,i) + a2*xd(:,i)+ a3*xdd(:,i)) + ...
                      C*(a1*x(:,i) + a4*xd(:,i)+ a5*xdd(:,i));
        x(:,i+1)    = Khat\fhat;
        xdd(:,i+1)  = a0*(x(:,i+1)-x(:,i)) - a2*(xd(:,i)) - a3*(xdd(:,i));
        xd(:,i+1)   = xd(:,i) + a6*xdd(:,i) + a7*xdd(:,i+1);
    end
end






