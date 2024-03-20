%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       W2SRS computes and visualises SRS S.D. for 1-DOF randomly excited systems
%
%       See also: XFIG, INTERPC4, TABC4
%
%   UPDATES
%       [-] separate useful functions as standalone ones in Mechanics
%       [-] tabc4 -> @tab when subroutines are ready
%       [-] plotting of xdd/xd/x TFs simultaneously with appropriate linestyle/legend
%       [-] window box with RMS of table and output PSD
%       [-] surface plot of SRS-generating TFs(?)
%
%   VERSION
%   v1.0 / 10.11.22 / V.Y.
%  ------------------------------------------------------------------------------------------------

clearvars 

% System, input, SRS 
    fn = 11.25;
    z = 0.01;
    Q = 10;

    range = [1 1e4];                                                    % TF / plotting range
    rs    = [20 2e3];                                                   % SRS integration limits
    order = 0;                                                          % integrals of xdd
    ratio = 10.^linspace(-2, 2, 5);                                     % spectra of W->SDOF(fn,zeta)->SDOF(k*fn,Q)

    T = [ 20    nan    nan  +6	    
          50    nan    nan  0
          800   nan    nan  nan
          2000  0.04   +6   nan  ];

% Plot options
    op.x = {'xy',1,'g',1,'b',0};                                        % xfig
    op.f = {range}; %{range,'meshdensity',500};                         % fplot
    op.b = true;                                                        % integration limit box
    
% [ FUNCTIONS ]
    W   = @(f)interpc4(T,f);                                            % input spectrum
    H   = genTF(fn,z);                                                  % TF H(f,n) for integral order n
    HW  = genTF2W(fn,z,W);                                              % output spectrum HW(f,fn,n)
    HX  = genTF2W([],0.5/Q,HW);
    SRR = @(f)sqrt(integral(@(x)HX(x,f),rs(1),rs(2)));                  % S.D. of SRS

% [ PLOTS ]
    warning off
    t = tiledlayout(2,2);
    
    % TFs xdd -> xdd/xd/x
    ax = xfig(nexttile(t,1),op.x{:});
        p.H = arrayfun(@(n)fplot(@(x)abs(H(x,n)), op.f{:}, color=col('bc')), order);

        drawBox(ax,rs,range,op.b)
        xlabel('$f$ [Hz]');
        ylabel('TF');
        legend(p.H, '$H_{\zeta,f_n}$');

    % Response PSD
    ax = xfig(nexttile(t,3),op.x{:});
        p.HW = arrayfun(@(n)fplot(@(x)HW(x,fn,n), op.f{:}, color=col('k')), order);
        p.W  = fplot(W, range, color=col('applegreen'));

        lgdv = [p.W p.HW(1)];
        if ~isempty(ratio)
            p.HX = arrayfun(@(k)fplot(@(x)HX(x,k*fn,0), op.f{:}, color=col('gainsboro')), ratio);
            ax.Children = [p.HX p.W p.HW];  % reorder plotted objects
            lgdv = [lgdv, p.HX(1)];
        end

        drawBox(ax,rs,range,op.b)
        xlabel('$f$ [Hz]');
        ylabel('PSD');
        legend(lgdv, ["$W$" "${|H_{\zeta,\,f_n}|}^2 W$" "${|H_{Q,\,k f_n}|}^2 {|H_{\zeta,\,f_n}|}^2 W$"]);

        % sqrt(integral(W,rs(1),rs(2)))
        % sqrt(integral(@(x)HW(x,fn),rs(1),rs(2)))
        % SRR(fn)
        
    % SRS of response PSD
    ax = xfig(nexttile(t,2,[1 2]),op.x{:});
        p.SRR = fplot(SRR, range, color=col('rc'));

        drawBox(ax,rs,range,op.b)
        xlabel('$f$ [Hz]');
        ylabel('SRS [S.D.]');
        legend(p.SRR,'$\displaystyle\int_{f_0}^{f_1}{{|H_{Q,\,f}|}^2 {|H_{\zeta,\,f_n}|}^2 W}$')
        
    % Set style and export
    % exportgraphics(gcf,'fig.pdf','contenttype','vector')
    warning on
    


%  ------------------------------------------------------------------------------------------------

% Draw boxes on integration patch
function drawBox(ax,rs,range,req)
    if any(rs~=range) && req
        yy = ax.YLim;
        line( rs([1 1]), yy, linestyle=':', color=col('l') );
        line( rs([2 2]), yy, linestyle=':', color=col('l') );
        patch('xdata', rs([1 2 2 1]) ,'ydata', yy([1 1 2 2] ), ...
                facealpha=.06, edgealpha=0, facecolor=col('coolgrey'));
        ax.Children = ax.Children([end,1:end-1]);
    end
end
    
% Generate handles from TF
function g = genTF(fn,zeta)
    function y = G(f,n)
        if nargin<2, n=0; end
        y = TF(f,fn,zeta,n);
    end
    g = @G;
end

% Generate handles for TF2W
function g = genTF2W(fn,zeta,W)
    function y = G(f,fnx,n)
        if nargin<3 || isempty(n),   n=0;    end
        if nargin<2 || isempty(fnx), fnx=fn; end
        y = TF2W(f,fnx,zeta,W,n);
    end
    g = @G;
end

% 1-DOF transfer function xdd -> xdd/xd/x
function y = TF(f,fn,zeta,n)
    if nargin<4 || isempty(n), n=0; end
    i = sqrt(-1);
    y = -(2*zeta*(i*f/fn) + 1) ./ ((i*f/fn).^2 + 2*zeta*(i*f/fn) + 1) ./ (2*pi*i*f).^n;
end

% Random response 
function y = TF2W(f,fn,zeta,W,n)
    if nargin<5 || isempty(n), n=0; end
    y = abs(TF(f,fn,zeta,n)).^2 .* W(f);
end


%  ------------------------------------------------------------------------------------------------
%{
% System
    f = 20; 
    z = 0.01;
    m = 1;
    w = 2*pi*f;
    k = m*w^2;
    c = 2*sqrt(k*m)*z;

% Laplace transforms
    FS = @(s) -(c*s + k) ./ (m*s.^2 + c*s + k);
    FN = @(s) -(2*z*w*s + w^2) ./ (s.^2 + 2*z*w*s + w^2);
    fs = @(f) FS(sqrt(-1)*2*pi*f);
    fn = @(f) FN(sqrt(-1)*2*pi*f);

% [Check] xdd/xd/x equal at freq. 0.5/pi
    [y,t] = sineSweep(.499/pi,.501/pi,output=0:2,rampCycles=2,tmax=200);
    xfig;plot(t,y);

% [Check] phase
    xfig(x=1,g=1,gm='y');
    arrayfun(@(n)fplot(@(x)rad2deg(angle(H(x,n))),range),0:2);
%}
%  ------------------------------------------------------------------------------------------------



