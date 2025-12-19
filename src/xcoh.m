function X = xcoh(A,B,dim)
%  ------------------------------------------------------------------------------------------------
%   DESCRIPTION
%       X = XCOH(A,B,dim)
%
%   INPUTS
%       A, B        arrays of compatible sizes
%       dim         dim to evaluate correlation along
%
%   VERSION
%   v1.0 / 25.10.25 / V.Y.  from dummySC v1.1
%  ------------------------------------------------------------------------------------------------

arguments
    A
    B
    dim {mustBeInteger} = 1
end

X = abs( dot(A,B,dim) ); 
X = X ./ ( vecnorm(A,2,dim) .* vecnorm(B,2,dim) );
X = X .^ 2;
