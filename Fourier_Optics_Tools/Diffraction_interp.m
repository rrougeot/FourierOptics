function [ Wd ] = Diffraction_interp(alpha, beta, x, y, dx, dy, z, s_a, Axe, Amplitude )

% This function creates the diffracted wavefront in plane A 
% Use a linear interpolation from pre-coomputed radial Hankel Transformation
%
% Inputs:
% - alpha, beta: anglular coordinates of the point source [rad]
% - x, y: NxN meshgrid, in pixel
% - dx, dy: lateral misalignment [m]
% - z: distance plane O - plane A [m]
% - s_a: spatial sampling [m]
% - Axe: radial coordinates of pre-computation [m]
% - Amplitude: complex amplitude from pre-computation
%
% Outputs:
% - Wd: array NxN including the diffracted wave-front

r_eff = sqrt((alpha*z + s_a.*x - dx).^2 + (beta*z + s_a.*y - dy).^2);

Wd = interp1(Axe, Amplitude, r_eff, 'linear',1);

end

