function [x, y, dS, Ntot] = Adaptive_sampling(type, Nr_sun, Ntheta_sun, g, CDF)

% This function creates the adaptive sampling from the radial CDF
% Three types of sampling: 
%   0: Regular even radial sampling along the segment [0,1]
%   1: the radial sampling points, along the segment [0,1]
%   2: the two dimension sampling points of the unitary disk
%
% Inputs:
% - type: 0, 1 or 2, to choose the type
% - Nr_sun: number of points along the segment ]0,1], excluding 0
% - Ntehta_sun: number of rotation for the circular integration, for type 0 or 1
% -       number of points in the most external ring, for type 2
% - g: >=1, defines the slope to decrease the number of points per ring
% - CDF: the CDF to be used
%
% Outputs:
% - x: x-coordinates in [0,1]. If type 1), the radial points
% - y: x-coordinates in [0,1]. If type 1), 0
% - dS: the elementary surfaces r dr for type 1), r dr dtheta for type 2)
% - Ntot: total number of points

if (type == 0)
    
    x = (0:1:(Nr_sun-1))/Nr_sun;
    y = zeros(1, Nr_sun);
    dS = x/Nr_sun * (2*pi/Ntheta_sun);
    Ntot = Nr_sun;
    
else

%% Adaptive radial sampling from the CDF
   
    x_cdf = 0:1:(length(CDF)-1);
    x_cdf = x_cdf/length(CDF);

    dy_int = 1/Nr_sun;
    y_int = (0:1:Nr_sun) * dy_int;
    x_int = interp1(CDF, x_cdf, y_int, 'linear', 0);

    r = zeros(1,Nr_sun+1);
    dr = zeros(1,Nr_sun+1);

    for i=2:(Nr_sun+1)    
        r(i) = (x_int(i-1) +  x_int(i))/2;
        dr(i) = x_int(i)-x_int(i-1);   
    end

    if(type == 1)
        
        x = r;
        y = zeros(1, Nr_sun+1);
        dS = r.*dr * (2*pi/Ntheta_sun);
        Ntot = Nr_sun+1;
        
    
    elseif(type == 2)
    
%% Adaptive 2D sampling

    x = 0;
    y = 0;
    dS = 0;
    Ntot = 1;
    
    for i=(Nr_sun+1):(-1):2
    
        N = ceil(Ntheta_sun * (r(i)/r(Nr_sun+1))^g);
        if(N <3)
            N = 3;
        end
        l = 0:1:N;
        a = pi/N.*l;

        x = [x, r(i)*cos(a)];
        y = [y, r(i)*sin(a)];
        
        dS = [dS, r(i) * dr(i) * (pi/N) * ones(1,length(a))];
        Ntot = Ntot + (N+1);
    end

    dS(abs(y)<10e-10) = dS(abs(y)<10e-10)/2;
    
    end
    
end

    