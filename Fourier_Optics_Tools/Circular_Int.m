function [I_int]=Circular_Int(vx,vy,I,Ntheta_sun)

% This function performs the circular integration
%
% Inputs:
% - vx: vector of pixel in x
% - vy: vector of pixel in y
% - I: 2D array of the 2D image
% - Ntheta_sun: angular step size for integration
%
% Outputs:
% - I_int: radial vector with integration

N = length(vx);
t = (1:Ntheta_sun)*2*pi/Ntheta_sun;
r = vx((N/2+1):(N-2))';

I_int=zeros(size(r));

for i=1:Ntheta_sun
    [Xi,Yi]=pol2cart(t(i),r);
    I_int = I_int + Interpgrid(vx,vy,I,Xi,Yi); %/Ntheta_sun;
end

end