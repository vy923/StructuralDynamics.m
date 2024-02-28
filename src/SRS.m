function S = SRS(xddb,dt,opts)
%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       S = SRS(xddb,dt,opts)
%
%       See also:       mustBeEqualDims, mustBeStartString, v2string, newmarkBeta     
%       External:       v2struct
%       Related:        squareWaveSRS
%
%   OUTPUTS
%       S
%           .f          natural frequencies vector of SDOF systems
%           .side.meth  output SRS, e.g. S.abs.vel
%           .x/xd/xdd   rel.disp/rel.vel/abs.acc [ndof x timesteps x loadcases] 
%
%   INPUTS
%       xddb            [m/s] base shake acceleration
%       dt              [s] time increments in xdd
%
%       opts{:}         name-value pairs
%           t           [s] time vector for xdd, gets resampled to 1/dt Hz
%           methInterp  time vec interpopation method, default = 'linear'
%           methSRS     any combination of {'disp','vel','acc'}, default = 'acc'
%           signSRS     any combination of {'+','-','abs'}, default = 'abs'
%           f0          [Hz] start frequency, default = 100
%           f1          [Hz] end frequency, default = 10000
%           zeta        modal damping of SDOF systems, default = 0.05
%
%           noct        sample points per octave, default = 25
%           ndof        total sample points, overrides nPerOct
%
%           xout        vector with output derivative orders, default = []
%               []      output variable X is empty
%               [0 1 2] compute x/xd/xdd e.g. [0 2] requests only x/xdd
%
%   EXT. CALLS
%       v2struct        unpack opts to workspace  
%
%   VERSION
%   v1.4 / 28.02.24 / --    added square wave + n-period sine examples
%   v1.3 / 27.06.22 / --    updated argument validation
%   v1.2 / 24.06.22 / V.Y.
%  ------------------------------------------------------------------------------------------------

arguments
    xddb (:,:) {mustBeReal}
    dt (1,1) {mustBePositive}
    opts.f0 (1,1) {mustBePositive} = 100
    opts.f1 (1,1) {mustBePositive} = 10000
    opts.zeta (1,1) {mustBePositive} = 0.05
    opts.noct (1,1) {mustBeInteger,mustBePositive} = 25
    opts.ndof (1,1) {mustBeInteger,mustBePositive}
    opts.t (:,1) {mustBeReal,mustBeEqualDims(opts.t,xddb,1)}
    opts.methInterp {mustBeMember(opts.methInterp,["linear","spline","makima"])} = 'linear'
    opts.methSRS {mustBeStartString(opts.methSRS,["acceleration","velocity","displacement"])} = 'acc'
    opts.signSRS {mustBeStartString(opts.signSRS,["positive","negative","absolute"])} = 'abs'
    opts.xout {mustBeMember(opts.xout,0:2)} = []
end

v2struct(opts);                                                                                 % Alternative: cellfun(@(n) assignin('caller',n,getfield(opts,n)),fieldnames(opts))
v2string(methInterp,methSRS,signSRS);                                                           % convert class of arguments to string

% Sampling points  
    if ~exist('ndof','var')
        ndof = ceil(noct*log(f1/f0)/log(2)) + 1;                                                % [sampling rate] x [number of octaves] + 1
    end
    S.f = exp(linspace(log(f0),log(f1),ndof))';                                                 % [ndof x 1] natural frequencies

% Create SDOF systems
    w = 2*pi*S.f;
    M = ones(ndof,1);
    K = w.^2;
    C = 2*zeta*w;                                                                               % 2*sqrt(k*m), general matrix version in modalDamping.m

% Resample base shake if time vector is given 
    if exist('t','var')
        xddb = interp1(t,xddb,[t(1):dt:t(end)]',methInterp);                                    
    end

% Loop over loadcases
    for lcase = size(xddb,2):-1:1                                                               % backwards to initialise correct output array sizes directly

        f = -xddb(:,lcase)'.*M;                                                                 % [ndof x timesteps] base shake equivalent force
        [x,xd,xdd] = newmarkBeta(M,C,K,f,dt,0.25,0.5,zeros(ndof,1));                            % integration
        xdd = xdd + xddb(:,lcase)';                                                             % absolute acceleration
     
    % Collect required transient outputs
        fieldNamesXout = ["x","xd","xdd"];

        for i = 1:numel(xout)                                                                   % Repeated string construction >> faster than using 3-dim x* arrays
            iFieldName = fieldNamesXout(xout(i)+1);
            S.(iFieldName)(:,:,lcase) = eval(iFieldName);                                       % Initialised directly to final size
        end
    
    % Output structure prep
        fnamesMeth = {'displacement','velocity','acceleration'};                                % Order must be the same as for fieldNamesXout
        fnamesSign = {'positive','negative','absolute'};
        evalExprSign = ["max(abs(X),[],2)"];
        evalExprMeth = ["X.*(X>=0)","X.*(X<=0)","X"];
    
    % Calculate and assign required SRS
        for i = 1:numel(fnamesSign)
            for k = 1:numel(fnamesMeth)
                if isSubstring(signSRS,fnamesSign(i)) && isSubstring(methSRS,fnamesMeth(k))     % arg validation already implies match='start'
                    evalExpr = strrep(evalExprSign,"X",evalExprMeth(i));
                    evalExpr = strrep(evalExpr,"X",fieldNamesXout(k));
                    S.(fnamesSign{i}(1:3)).(fnamesMeth{k}(1:3))(:,lcase) = eval(evalExpr);
                    %disp(evalExpr)
                end
            end
        end 

        clearvars x* -except xout xddb

    end % for lcase 


%  ------------------------------------------------------------------------------------------------
%{
% [EXAMPLE 1], Square wave SRS

dt = 1e-6;
sf = 200;                                                                                       % square wave frequency [Hz]
n  = 2;                                                                                         % sq. wave periods
k  = 1;                                                                                         % zero signal at end of excitation k times longer than input sq. wave
f  = [1e1 1e4];                                                                                 % frequency range
zeta = [.01:.02:.07];                                                                           % modal damping

padw = @(x,k) [cvec(x(1:end-1)); zeros(length(x)*k,1)];                                         % pad with zeros k times the length of input
swav = @(f,dt,n) square(2*pi*(0:sf*dt:n));                                                      % n-period square wave of freq. sf and sampling rate dt
wv = padw(swav(sf,dt,n),k);

t = tiledlayout(3,1);
[ax,fig] = xfig(nexttile(t,1,[2 1]),xy=1,b=1,c=1);

    fig.Position = [100 100 698 392].*[1 1 1.2 1.8];
    for z = zeta
        S = SRS(wv,dt,noct=200,zeta=z,f0=f(1),f1=f(2));
        plot(S.f,S.abs.acc)
    end

    legend(arrayfun(@(i)string(['$\zeta=' sprintf('%g',i) '$']),zeta))
    ax.XLabel.String = 'Frequency, Hz';
    ax.YLabel.String = 'SRS';

[ax,fig] = xfig(nexttile(t,3),b=1,c=1);

    plot(0:dt:size(wv,1)*dt-dt,wv);
    ax.XLimitMethod = 'tight';
    ax.YLimitMethod = 'padded';
    ax.YLabel.String = 'Acceleration';
    ax.XLabel.String = 'Time, s';


% [EXAMPLE 2], n-period sine SRS

dt = 1e-6;
sf = 200;                                                                                       % frequency [Hz]
n = 0.5;                                                                                        % number of periods
k = 1;                                                                                          % zero signal padding
w = 2*pi*sf;
x = 0:dt*w:2*pi*n;

S = SRS([sin(cvec(x)); zeros(k*length(x),1)],dt,noct=200,f0=10);

ax = xfig(1,b=1,x=1,g=2,c=1);
    plot(S.f,S.abs.acc);
%}
%  ------------------------------------------------------------------------------------------------



























