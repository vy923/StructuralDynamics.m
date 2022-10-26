%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       T = mergec4(T1,T2)
%
%       See also:       tabc4
%       Related:        interpc4
%
%   VERSION
%       v1.0 / 24.10.22 / V.Yotov
%  ------------------------------------------------------------------------------------------------

function T = mergec4(T1,T2)




% discountinuity single point match -> split single point
% discontinuity double point match -> merge as per min/max/sum rule
% 