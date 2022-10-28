%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       [X,Xb] = integratec4(T,f)
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
%}