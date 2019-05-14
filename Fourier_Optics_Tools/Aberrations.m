function [Wab] = Aberrations(sigma, s, Rp, x, y)

% This function builds the aberration from Zernike polynomials
%
% Inputs:
% - sigma: vector 1x8 with the RMS level of the aberrations [lambda unit]
%          tilt_x, tilt_y, defocus, ast_x, ast_y, coma_x, coma_y, spherical
% - s: spatial sampling in plane A [m]
% - lambda: wavelength [m]
% - Rp: radius of the pupil [m]
% - x, y: meshgrid of plane A, in pixel
%
% Outputs:
% - Wab: array NxN including the wave-front perturbation

r = sqrt(x.^2 + y.^2);
r = r*s/Rp;
t = atan2(y,x);

Wab  = zeros(size(r));

if(sigma(1) ~= 0)
    
        % Tilt in x-axis
        Wab = Wab + sigma(1)* 2*r.*cos(t);
        
elseif(sigma(2) ~= 0)
    
        % Tilt in y-axis
        Wab = Wab + sigma(2)* 2*r.*sin(t);

elseif (sigma(3) ~= 0)
    
        % Defocus
        Wab = Wab + sigma(3)* sqrt(3)*(2*r.^2 - 1);
        
elseif(sigma(4) ~= 0)
    
        % Astigmatism in x-axis
        Wab = Wab + sigma(4)* sqrt(6)*r.^2 .*cos(2*t);

elseif(sigma(5) ~= 0)
    
        % Astigmatism in y-axis
        Wab = Wab + sigma(5)* sqrt(6)*r.^2 .* sin(2*t);        
        
elseif(sigma(6) ~= 0)
    
        % Coma in x-axis
        Wab = Wab + sigma(6)* sqrt(8)*(3*r.^3 - 2*r) .* cos(t);

elseif(sigma(7) ~= 0)
    
        % Coma in y-axis
        Wab = Wab + sigma(7)* sqrt(8)*(3*r.^3 - 2*r) .* cos(t);        
        
elseif(sigma(8) ~= 0)
    
        % spherical
        Wab = Wab + sigma(8)* sqrt(5)*(6*r.^4 - 6*r.^2 +1);     
end

Wab = exp( 2*1i*pi.* Wab);

end
        
        