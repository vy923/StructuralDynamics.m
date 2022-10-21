% Square wave response, SRS

dt = 1e-5;
sf = 150;                                                           % square wave frequency [Hz]
n = 2;                                                              % sq. wave periods
k = 1;                                                              % zero signal at end of excitation k times longer than input sq. wave
padw = @(x,k) [cvec(x); zeros(length(x)*k,1)];                      % pad with zeros k times the length of input
swav = @(f,dt,n) square(0:f*dt:2*pi*n);                             % n-period square wave of freq. f and sampling rate dt

[ax,fig] = xfig(xy=1,g=1,b=1,n=1);
fig.Position = [100 100 698 392].*[1, 1, .7, .8];                   % second vector is for easier scaling, first vector is position on screen(x,y), scaling(x,y)

wv = padw(swav(sf,dt,n),k).*[1 1];

S = SRS(wv,dt,noct=200,zeta=.05,f0=1,f1=1e3);
    plot(S.f,S.abs.acc)
S = SRS(wv,dt,noct=200,zeta=.03,f0=1,f1=1e3);
    plot(S.f,S.abs.acc)
S = SRS(wv,dt,noct=200,zeta=.01,f0=1,f1=1e3);
    plot(S.f,S.abs.acc)

    legend('$\zeta=0.01$','$\zeta=0.03$','$\zeta=0.05$')
    ax.XLabel.String = 'Frequency, Hz';
    ax.YLabel.String = 'SRS, g';
    ax.Title.String = '\textbf{SRS of a 1 g 50 Hz square wave}';

    exportgraphics(fig,'figVec.pdf','contenttype','vector')

[ax,fig] = xfig(n=2,b=1);
fig.Position = [100 100 698 392].*[1, 1, 1.0, 1.0]; 

    plot([0:sf/2*dt:2*pi*n]/(4*pi/0.08),wv);

    ax.XLimitMethod = 'tickaligned';
    ax.YLimitMethod = 'padded';
    ax.Title.String = '\textbf{50 Hz square wave input, 2 periods}'
    ax.YLabel.String = 'Acceleration, g';
    ax.XLabel.String = 'Time, s';

    exportgraphics(fig,'figVec1.pdf','contenttype','vector')
    print('-painters','-dsvg','figVec1');                           % .svg option 
    