%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       T = scalec4(T,opts)
%
%       See also:       tabc4, interpc4, integratec4
%
%   INPUTS
%       T               -------
%       opts{:}
%           rms         scalar or vector of RMS values 
%
%   OUTPUTS
%       T               [nband + 1 x 4 x nlevel] list of scaled PSD tables
%
%   VERSION
%       v1.0 / 16.10.22 / V.Yotov
%  ------------------------------------------------------------------------------------------------

function T = scalec4(T,level,opts)

arguments
    T
    level %{mustBeVector,mustBePositive} = 1
    opts.rms %{mustBeVector,mustBePositive} = []
end


    T = tablePSD(T);                                                            % Validation


% integratePSD function    

% ASD.*10.^(level_dB/10)

% - min / max / sum PSD (from Xb formula...)
% - SPL with zero values?

%{
T = [   32      34
        50      69
        63      77
        40      94
        80      86
        100     114
        125     112
        1000    4
        160     108
        200     111
        250     78
        315     55
        400     43
        600     33
        630     19
        800     11
        1250    2
        1600    1
        2000    1
        2500    0  ]

integratePSD([f',interpPSD(T,f)'])

%}


















