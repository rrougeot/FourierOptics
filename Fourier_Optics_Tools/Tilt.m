function [T] = Tilt(lambda, alpha, beta, x, y, s)

% This function returns  the planar tilt wave
%
% Inputs:
% - lambda: wavelength [m]
% - alpha, beta: angles [rad]
% - x, y: array NxN pixel meshgrid
% - s: spatial sampling [m]
%
% Outputs:
% - T: array NxN with the tilt

T = exp( -2*1i*pi/lambda * alpha * s .* x) ...
    .* exp( -2*1i*pi/lambda * beta * s .* y);
    
end