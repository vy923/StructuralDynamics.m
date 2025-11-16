%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       [X,Xb] = integratec4(T)
%
%       Integrates a c4 table given by T and returns total and band RMS values in X, Xb
%       Note that for PSD tables X = sqrt(sum(Xb.^2))
%
%       See also:       tabc4
%       Related:        scalec4, splitc4
%
%   UPDATES
%       - optional argument f which automatically calls scalePSD, etc.
%       - robust handling of Inf slopes
%
%   VERSION
%       v2.0 / 23.10.22 / --        integrates discontinuous tables
%       v1.0 / 14.10.22 / V.Yotov
%  ------------------------------------------------------------------------------------------------

function [X,Xb] = integratec4(T)

% Validate/complete T
    T = tabc4(T);

% Preallocatons
    W1 = T(1:end-1,2);
    RS = T(1:end-1,4);
    F1 = T(1:end-1,1);
    F2 = T(2:end,1);

% Computation
    C = 1 + RS/(10*log10(2));
    A = W1.*F1./C .* ((F2./F1).^C - 1);

    X  = sqrt(sum(A));
    Xb = sqrt(A);


%  ------------------------------------------------------------------------------------------------
%{
% Example
T =  [  20          0.0032
        50          0.02
        150         0.02
        200         0.01
        500         0.01 
        2000        0.0006  ];
[X,Xb] = integratec4(T)

% In progress
T1  = T.*[1, 9.81^2];                                               % in SI units
T2  = [20:1:2000; interpc4(T1,20:1:2000)]';                         % in SI units, interpolated @ 1 Hz (optional step, you can say T2 = T1 to skip)
Q   = @(p)T2.*[ones(size(T2,1),1) (T2(:,1).*(2*pi)).^-p];           % Q(0)/Q(2)/Q(4) is ASD/VSD/DSD table

disp([integratec4(Q(0)) integratec4(Q(2)) integratec4(Q(4))])       % RMS acc/vel/disp in SI units

% W   = @(f,p)interpc4(Q(p),f);                                       % ignore
% RMS = @(p)sqrt(integral(@(x)W(x,p),20,2000));                       % ignore

%}








