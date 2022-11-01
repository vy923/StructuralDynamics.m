function [y,t,f,R] = sineSweep(f0,f1,opts)
%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       [y,t,f,R] = sineSweep(f0,f1,opts)
%
%       See also:       --
%       External:       sinintx, cosintx, fresnelsx, fresnelcx, v2struct
%       Related:        rotatingForce, forcedAcce1D
%
%   OUTPUTS
%	    y               3 x n time signal matrix [a,v,u]
%	    t               time in seconds
%       f               instantaneous frequency as a function of time
%       R               structure with ramp data
%           R.tr        time
%           R.nR        ramp data points number
%           R.gs        smoothstep function from 0 to +inf
%
%   INPUTS
%	    f0              start [Hz]
%       f1              end [Hz]
%
%       opts{:}         name-value pairs
%           output      vector with output integral orders, default = -1
%               [-1] 	no drift correction
%               [0 1 2] compute a, v, u, e.g. [0 2] requests only a and u
%
%           tmax        sweep time in seconds, overrides sweepRate 
%                       if rampOffset = false, tmax = ramp + clean sweep 
%                       if rampOffset = true, tmax = clean sweep only
%
%           sweepRate   octaves per minute, default = 2
%
%           sweepType   default = 'log'
%               'log'   logarithmic sweep
%           	'lin'   linear sweep
%
%           sampleRate  points per period, default = 20
%
%           rampOffset  default = 0
%               true    ramp up is complete exactly at frequency f0
%               false   ramp up starts at f0, clean sweep starts at higher fr.
%
%           rampCycles  num. cycles for g(t) to go from 0 to 1, rounded down to 
%                       floor(sampleRate*rampCycles), default = 0
%
%           rampOrder   clamp function g(t) default = -1
%               n < 0	half sine
%               n >= 0  2n+1-st order smoothstep polynomial
%
%           limCond     drift correction condition, default = 0      
%               inf     velocity v(inf)==0
%               ~=inf   displacement u(tlim)==0, where tlim is the last
%                       acceleration zero crossing
%
%           phaseShift  [rad], default = 0
%
%   EXT. MODDED
%       (sin/cos)intx 	trigonometric integrals (Erik Koene / modded VY)
%       fresnel(sx/cx) 	normalised Fresnel integrals (John D'Errico / modded VY)
%
%   KNOWN BUGS
%       quadgk error with tr == 0 and out ~= -1
%
%   VERSION
%       v2.1 / 24.06.22 / --        updated opts unpacking to workspace
%       v2.0 / 01.03.22 / V.Yotov
%  ------------------------------------------------------------------------------------------------

arguments
    f0 (1,1) double
    f1 (1,1) double
    opts.output {mustBeMember(opts.output,-1:2)} = -1
    opts.tmax {mustBeReal} = []
    opts.sweepRate (1,1) {mustBeReal} = 2
    opts.sampleRate (1,1) {mustBeInteger} = 20
    opts.sweepType {mustBeMember(opts.sweepType,{'lin','log'})} = 'log'
    opts.rampOffset {mustBeMember(opts.rampOffset,[0 1])} = false
    opts.rampCycles (1,1) {mustBeReal} = 0
    opts.rampOrder {mustBeInteger} = -1
    opts.limCond (1,1) double = 0
    opts.phaseShift (1,1) double = 0
end
v2struct(opts);                                                                             % Legacy unpack: cellfun(@(n) assignin('caller',n,getfield(opts,n)),fieldnames(opts))

% Parameters
    r = sampleRate;                                                                         % to simplify notation
    oct = log(f1/f0)/log(2);                                                                % number of octaves in sweep
    phs = phaseShift;                                                                       % phase shift in radians 
    
    if isempty(tmax)
        tmax = oct/sweepRate*60;                                                            % [s]
    else
        sweepRate = oct/tmax*60;                                                            % [oct/min]
    end

    lin = strcmp(sweepType,'lin');                                                          % bool switch 
	nRoffs = rampOffset*floor(rampCycles*r);                                                % Additional sampling points if rampOffset==true

% Time points, frequencies, ramp function
    if lin       
        a    = 0.5*(f1-f0)/tmax;
        c 	 = pi*f0^2/2/a; 
        phi  = @(t) 2*pi*(a*t.^2+f0*t);                                                     % phase [rad]
        ninv = @(x) (-f0*sqrt(r)+sqrt(4*a*x+r*f0^2))/(2*a*sqrt(r));                         % inverse of cumulative number of samples function n(t)    
        t	 = ninv([-nRoffs:r/(2*pi)*phi(tmax)])';                                         % time point vector
        f	 = (f1-f0)/tmax*t+f0;                                                           % instantaneous frequency 

        f0   = f(1);                                                                        % Updates f0/c/phi when rampOffset==true
        c    = pi*f0^2/2/a; 
        phi  = @(t) 2*pi*(a*t.^2+f0*t) + phs;
    else
        b	 = (f1/f0)^(1/tmax);                                                            % exponent base
        a    = 2*pi*f0/log(b);
        phi  = @(t) a*(b.^t-1);
        ninv = @(x) log(1+2*pi/a*x/r)/log(b);
        t	 = ninv([-nRoffs:r/(2*pi)*phi(tmax)])';
        f  	 = f0*b.^t;
        
        a    = 2*pi*f(1)/log(b);                                                            % Updates a/phi when rampOffset==true
        phi  = @(t) a*(b.^t-1) + phs;                                                       % add required phase shift
    end 
	
% Check if requested ramp offset is possible
    if ~isreal(ninv(-nRoffs))
        error('ERROR: Reduce number of ramp cycles')
    end

% Clean-up time vector
    t = t-t(1);                                                                             % start time = 0, as it is negative if rampOffset==true
	t(abs(t)<eps(1e2)|t<0) = 0;                                                             % clear numerical zeros
    nR = min(floor(rampCycles*r+1),numel(t));                                               % ramp points, incl. 0 and tr
    tr = t(nR);                                                                             % ramp time

% Generate ramp function
    g  = smoothstep(rampOrder,tr);
    Rf = optAnon(@(x)g(x).*sin(phi(x)));

% Ramp data, only for output  
    R.gs = @(x) (x<tr).*g(x) + (x>=tr);                                                     % Smoothstep extended to +inf  
    R.tr = tr;
    R.nR = nR;

% ------ no corrections for v, u
if ismember(-1,output)

    y = [ Rf(t(1:nR-1));  sin(phi(t(nR:end))) ];

% ------ corrections for v or u
else 
    
% Antiderivative without integration constant
    if lin  
        intA = @(t)( cos(c)*fresnelsx((2*a*t+f0)/sqrt(a)) - ...
                     sin(c)*fresnelcx((2*a*t+f0)/sqrt(a)) ) / (2*sqrt(a));
    else
        intA = @(t)( cos(a)*sinintx(a*b.^t) - sin(a)*cosintx(a*b.^t) ) / log(b);
    end                      
    optq = {'reltol',1e-8,'MaxIntervalCount',1e6};                                          % [quadgk] parameters
    opti = {'reltol',eps(1e6)};                                                             % [integral] parameters

% Boundary conditions
    h = msgbox(sprintf('Check for warnings, error bound << 1.0 is acceptable'), ...
        'BC parameter calc','warn');
    
    w = @(x,k)-k.*sin(pi*x/tr).^2;                                                          % offset function
    
    if limCond == inf
        cr  = intA(inf)-intA(tr);                                                           % [tr,inf], sin(phi(x))
        cl  = quadgk(Rf,0,tr,optq{:});                                                      % [0,tr], g(x)*sin(phi(x))
        k   = -(cr+cl)/quadgk(@(x)g(x).*w(x,1),0,tr,optq{:});                               % Linear in k => cl + cr + k*int(g(x)*w(t,1)) = 0
    else
        tlim = t(end-rem(numel(t)-1,r/2));                                                  % last accel. zero crossing
        cr  = quadgk(@(x)(tlim-x).*Rf(x),0,tr,optq{:});                                     % double/Cauchy [tr,inf], sin(phi(x)), 
        cl  = quadgk(@(x)(tlim-x).*sin(phi(x)),tr,tlim,optq{:});                            % double/Cauchy [0,tr], g(x)*sin(phi(x))
        k   = -(cr+cl)/quadgk(@(x)(tlim-x).*g(x).*w(x,1),0,tr,optq{:}); 
    end   

% Ramped signal updated with w(x,k) term
    RF = optAnon(@(x)g(x).*(sin(phi(x))+optAnon(@(x)w(x,k))));

% Definite integrals
    c1 = quadgk(RF,0,tr,optq{:}) - intA(tr);                                                % As intA1(tr) = v(tr) = intA1L(tr) + c1
    intA1 = @(t)intA(t) + c1;
    
    c2 = quadgk(@(x)(tr-x).*RF(x),0,tr,optq{:});                                            % contribution from int(sin(phi(t))-RF(t))
    if lin
        c2 = c2 - (cos(2*pi*(a*tr.^2+f0*tr))/(4*pi*a)+(tr+f0/2/a).*intA1(tr));  
        intA2 = @(t)(cos(2*pi*(a*t.^2+f0*t))/(4*pi*a)+(t+f0/2/a).*intA1(t)) + c2;
    else
        intA2 = @(t) intN1(intA1,t,5*r,opti{:}) + c2;                                       % Correct if t(1)==t(nR) as intN1 computes F(t(:))-F(t(1))
    end

% Compute outputs
    y = zeros(length(t),3);
	y(:,1) = [RF(t(1:nR-1)); sin(phi(t(nR:end))) ];
    if ismember(1,output)
	    y(:,2) = [intN1(RF,t(1:nR-1),5*r,opti{:}); intA1(t(nR:end)) ];  
    end
    if ismember(2,output)
	    y(:,3) = [intN2(RF,t(1:nR-1),opti{:}); intA2(t(nR:end)) ];      
    end
    
end % if ismember(-1,out)

% Delete temp message box
    try delete(h); end

% Display solution settings
	fprintf('\nSweep freq: %.2f-%.2f Hz \n', [f(1) f1])
    fprintf('Sweep rate: %g oct/min\n', sweepRate)
    fprintf('Sweep time: %.2f s\n', t(end))
    fprintf('Samples:    %.0f \n', numel(t))
    fprintf('Ramp up:    %.2f-%.2f Hz \n', [f(1) f(nR)])
    fprintf('rampCycles: %.2f \n', (nR-1)/r)
    fprintf('rampTime:   %.2f s \n', tr)
    fprintf('rampFunc:   S%.0f \n', rampOrder) 	% --- UPDATE
    fprintf('limitCond:  %.0f \n', limCond)     % --- UPDATE
    
end % sineSweep



%% Local functions

% Unfold nested numeric functions
function g = optAnon(f,varargin)
    g = matlabFunction(sym(f),varargin{:});
end

% Returns a half sine or n-th order smoothstep ramp function mapped to 
% [0,xmax] for numeric and [0,t] for symbolic handles g and s
function [g,s] = smoothstep(n,xmax)
    syms s(x,t)
    if xmax == 0                                                                            % no ramping
        g = @(x)1.;
        return
    elseif n < 0                                                                            % sinusoid
        s(x,t) = sym(@(x,t)sin(.5*pi*x./t).^2);
    else                                                                                    % polynomial
        k = arrayfun(@(k)(-1)^k*nchoosek(n+k,k)*nchoosek(2*n+1,n-k),[0:n]);
        s(x,t) = sym(@(x,t)sum(k.*(x/t).^(n+1+[0:n])));
    end
	g = matlabFunction(s(x,xmax));                                                          % numeric function g(x)
end

% Vectorised integration with array limits T0, T1 
function y = intArr(fun,T0,T1,varargin)
    f = @(x)fun((T1-T0).*x+T0);                                                             % maps all inervals [T0,T1] to [0,1]
    y = (T1-T0).*integral(f,0,1,'ArrayValued',1,varargin{:});
end

% Numerical integration by division into r-point long blocks. Returns
% vector of definite integral values y = F(t(1:n))-F(t(1))
function y = intN1(fun,t,r,varargin)
    n = length(t);
	y = zeros(n,1);
    m = '[%.0f/%.0f] blocks';
	h = waitbar(0,sprintf(m,0,n/r),'name','Integrating...');   
	for i = 1:ceil((n-1)/r)
        k = (i-1)*r+1;                                                                      % i-th block start index
        s = min(r,n-k);  
        y(k+1:k+s) = y(k) + intArr(fun,t(k),t(k+1:k+s),varargin{:});                        % accumulate integral for i-th block
        if rem(i,ceil((n-1)/r/25))==0
            waitbar(i/ceil(n/r),h,sprintf(m,i,n/r))
        end
	end
	delete(h)
end

% Numerical double integration on [0,T1] using Cauchy's formula 
function y = intN2(fun,T1,varargin)
    h = msgbox(sprintf(['%.0f evaluations [ Cauchy ]' repmat(' ',1,30) '\n'], ...
        numel(T1)),'Integrating...');
    y = intArr(@(x)(T1-x).*fun(x),0,T1,varargin{:});
    delete(h)
end



%% Examples
%{
% ------------------------------------------------------
f0 = 45.0;
f1 = 46.0;
vars = {
'output'       [0 1 2]
'sampleRate'   20
'tmax'         20.0
'rampOffset'   true
'rampCycles'   20
'rampOrder'    -1
'sweepType'    'log'
'limCond'      0
}'; 

[y,t,f] = sineSweep(f0,f1,vars{:});

if 1 % --- for verification
a=y(:,1);v=y(:,2);u=y(:,3);
figure(2)
    clf
    ax = axes();
    set(ax,'LooseInset',ax.TightInset)
    hold on
	grid on
    %plot(t,a)
    plot(t,v)
    %plot(t,v-cumtrapz(t,a),'--')
    %plot(t,a)
    %plot(t,cumtrapz(t,cumtrapz(t,a)),'--')
legend
end

% ------------------------------------------------------
% Sweeps from f0 to f1 Hz, 5 sec total time, no ramping
% ------------------------------------------------------

f0 = 0.2;
f1 = 1.6;
r = 15; 
T = 5;

[y,t] = sineSweep(f0,f1,'sweeptype','lin','samplerate',r,'tmax',T);
    clf
    xlabel('Time [s]')
    ylabel('Amplitude')
    hold on
        m = numel(t)-1;
        p1 = plot(t(end-m:end),y(end-m:end),'r');
        scatter(t(end-m:end),y(end-m:end),'rs')

[y,t] = sineSweep(f0,f1,'sweeptype','log','samplerate',r,'tmax',T);
        m = numel(t)-1;
        p2 = plot(t(end-m:end),y(end-m:end),'k');
        scatter(t(end-m:end),y(end-m:end),'ks')
    legend([p1 p2],'lin','log','Location','southwest');

% -------------------------------------------------------------------
% Plot frequency vs time for 10-2000Hz sweeps at 4 oct/min using
% the default sampling of 20 points per period, no ramping
% -------------------------------------------------------------------

f0 = 10;
f1 = 2e3;
sr = 4;

[~,t,f] = sineSweep(f0,f1,'sweeptype','lin','sweepRate',4);
    clf
    xlabel('Time [s]')
    ylabel('Frequency [Hz]')
    % set(gca,'yscale','log')
    hold on
        p1 = plot(t,f,'r');

[~,t,f] = sineSweep(f0,f1,'sweeptype','log','sweepRate',sr);
        p2 = plot(t,f,'k');
    legend([p1 p2],'lin','log','location','southeast');
%}


%% Diagnostic code
%{
% -------------------------------------------------------------------
% Test vectorised integration by using random array limits and some 
% nontrivial test function, compare against arrayfun/loops
% -------------------------------------------------------------------

T0  = [-2 3].*[-1:.5:2]';                                                   % lower integration limits array
T1  = 0.5-rand(size(T0));                                                   % upper integration limits array
fun = @(x)3*x.^2+tanh(x)./x;                                                % integrand
fss = @(x)fun((T1-T0).*x+T0);                                               % shifted and stretched
Fvec = (T1-T0).*integral(fss,0,1,'ArrayValued',1,'reltol',eps(1e6));        % vectorised computation
Fref = arrayfun(@(x0,x1)integral(fun,x0,x1),T0,T1);                         % reference soln., neater than loop
    disp([Fref, Fvec, Fvec./Fref])

% -------------------------------------------------------------------
% Recursive vectorised computation of double integrals over same var., 
% slightly faster than arrayfun/loops of quad2d
% -------------------------------------------------------------------

v1 = [-3:1:4]';
v0 = 10*(.5-rand(size(v1)));
fun = @(x)3*x.^2+tanh(x).*x;
intArr = @(f,T0,T1)(T1-T0).*integral(@(x)f((T1-T0).*x+T0),0,1,'ArrayValued',1,'reltol',eps(1e6));
intRec = @(fun,T0,T1)intArr(@(x)intArr(fun,0,x),T0,T1); 
Fvec = intRec(fun,v0,v1);
Fref = arrayfun(@(v0,v1)quad2d(@(t,x)fun(x),v0,v1,0,@(x)x),v0,v1);
    disp([Fvec Fref])
%}


%% Extra plots
%{
% -------------------------------------------------------------------
% Offset function w(x) application options for ramp g(x)
% -------------------------------------------------------------------
n = -1;	% n<0 for sinusoid
tr = 1;	% max time
k = -.2; % max offset constant

w = @(x)k*sin(pi*x/tr).^2;

if n<0	g = @(x)sin(0.5*pi*x/tr).^2;
else	k = arrayfun(@(k)(-1)^k*nchoosek(n+k,k)*nchoosek(2*n+1,n-k),[0:n]);
        g = @(x)sum(k.*(x./tr).^(n+1+[0:n]));
end

clf 
hold on
col={'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE','#A2142F'};
co=[7 3 4];

fplot(@(x)g(x),[0 tr],'color',col{co(1)},'linestyle','-.')
fplot(@(x)-g(x),[0 tr],'color',col{co(1)},'linestyle','-.')
fplot(@(x)g(x)+w(x),[0 tr],'color',col{co(2)},'linestyle','-')
fplot(@(x)-g(x)+w(x),[0 tr],'color',col{co(2)},'linestyle','-')
fplot(@(x)g(x).*(w(x)+1),[0 tr],'color',col{co(3)},'linestyle','-')
fplot(@(x)g(x).*(w(x)-1),[0 tr],'color',col{co(3)},'linestyle','-')

p1 = fplot(@(x)0*x,[0 tr],'color',col{co(1)},'linestyle','-.','linewidth',0.5);
p2 = fplot(@(x)w(x),[0 tr],'color',col{co(2)},'linestyle','-','linewidth',0.5);
p3 = fplot(@(x)g(x).*w(x),[0 tr],'color',col{co(3)},'linestyle','-','linewidth',0.5);

legend([p1 p2 p3],'\pmg(x)','\pmg(x)+w(x)','g(x)(\pm1+w(x))','location','northwest')
legend([p1 p3],'\pmg(x)','g(x)(\pm1+w(x))','location','northwest')
title(sprintf(['Normalised S%.0f ' 'ramp function with offset'],n))
%}



