function x = timec4(W,fr,df)
%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       x = timec4(W,opts)
%
%       See also: tabc4, interpc4, cvec
%
%   INPUTS
%       W               handle for PSD table evaluation function
%       fr              freq range, default = [20 2000]
%       df              sampling frequency, default = diff(fr)/100
%
%   OUTPUTS
%       x.               
%           x           time series
%           t           time point vector
%           f           frequency vector
%           W           PSD vector
%
%   VERSION
%   v1.0 / 25.11.22 / V.Y.
%  ------------------------------------------------------------------------------------------------

if nargin < 2 || isempty(fr)
    fr = [20 2000];
end
if nargin < 3 || isempty(df)
    df = (fr(2)-fr(1))/100;                                             % 100 sampling points
end

N   = max([1 round(fr(2)/df)]);                                         % FD number of points
fs  = linspace(0,fr(2),N)';

Wfs = [zeros(nnz(fs<.8*fr(1)),1); cvec(W(fs(fs>=.8*fr(1))))];           % <----- remove .8
phi = 2*pi*rand(size(fs));                                              % random phase angles
an  = sqrt(Wfs*df).*exp(-sqrt(-1)*cvec(phi));                           % IFFT coefficient vector
an  = [an; zeros(size(an))];

x.t = 0.5/fr(2)*(0:2*N-1)';
x.x = sqrt(2)*real(fft(an));
x.W = 2/df*abs(ifft(x.x)).^2;
x.W = x.W(1:end/2);
x.f = fs;


%  ------------------------------------------------------------------------------------------------
%{
T =  [  20    0.006
        50    0.06
        150   0.06
        300   0.02
        700   0.02 
        800   0.06
        925   0.06
        2000  0.03];

W  = @(f)interpc4(T,f);
f  = [20 2000];
x  = timec4(W,f,0.1);

q = tiledlayout(5,4);

ax = xfig(nexttile(q,1,[3 4]),xy=1); 
    % fplot(W,f,color=col('k'),linestyle='--')
    plot(x.f(x.f>=f(1)),x.W(x.f>=f(1)), color=col('k'),linewidth=2.5)
    %xlabel('$f$')
    ax(1).YLimitMethod = 'padded';

    wind = 4*ceil(2*f(2)/f(1));
    [p,fp] = pwelch(x.x,wind,[],f(1):1:f(2),2*f(2));
    plot(fp,2*p,color=col('c2'))

    legend(['iFFT, $\sigma=' sprintf('%.2f',sqrt(integral(W,f(1),f(2)))) '$'], ...
           ['TD Welch']);

ax = xfig(nexttile(q,4,[2,3])); 
    plot(x.t,x.x,color=col('gainsboro'))
    legend('time series')
    %xlabel('$t$')

ax = xfig(nexttile(q,19,[2,1])); 
    nbins = 100;
    [fi,xi] = hist(x.x,nbins,'Normalization','probability');
    bh = barh(xi,fi/numel(x.x)*nbins/100,facecolor=col('gainsboro'),edgecolor=col('w'))

    pd = fitdist(x.x,'normal')
    pdxv = linspace(min(x.x),max(x.x),1e2);
    plot(pdf(pd,pdxv),pdxv,color=col('atomictangerine'));

    legend(['PDF, $\sigma=' sprintf('%.2f',std(x.x)) '$'],...
           ['Fitted normal']);

    exportgraphics(gcf,'fig.pdf','contenttype','vector')

%  ------------------------------------------------------------------------------------------------ 

% Other examples
    W0 = @(f,fu,a) (a./(a.^2 + f.^2)) ./ atan(fu/a);
    W  = @(f)W0(f,40,4);

    W  = @(f)2+sin(2*pi*log(f));


% SRS
    td = interp1(x.t,x.x,linspace(0,x.t(end),20*numel(x.t)),'makima')';
    S = SRS(td,.5/f(2)/20,f0=1)
    xfig(1,g=0,x=1);
        plot(S.f,S.abs.acc)

%  ------------------------------------------------------------------------------------------------
%}
















