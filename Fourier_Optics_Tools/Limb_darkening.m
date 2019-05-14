function [ B ] = Limb_darkening(a,b)

% Limb darkening function from Van Hamme 1993 @550nm
%
% Inputs:
% - a, b: angular coordinates of the points on the solar disk, in Rsun unit
%         a, b must be of the same size

% Outputs:
% - B : the limb darkening attenuation, same size as a, b

    r = sqrt(a.^2+b.^2);
   
    B = zeros(size(r));
    B(r<1) = 1 - 0.762*(1-sqrt(1-r(r<1).^2)) - 0.232*(1-r(r<1).^2).*log(sqrt(1-r(r<1).^2));
    
end

