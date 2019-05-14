function [Fresnel] = Fresnel(lambda, z, r, s)

% This function returns Fresnel quadratic phase wave
%
% Inputs:
% - lambda: wavelength [m]
% - z: distance [m]
% - r: array NxN with radius in pixel
% - s: spatial sampling [m]
%
% Outputs:
% - Fresnel: array NxN with the quadratic wave front

Fresnel = exp( 1i*pi/(lambda*z) .* (r*s).^2);
    
end
 