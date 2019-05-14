function [Wr, h] = Scattering(Harvey_param, seed, lambda, n, N, s_a)

% This function builds the roughness wave perturbation
%
% Inputs:
% - Harvey_param: 1x3 vectors with parameters of Harvey scatter model
%     [A [m3] amplitude, B [m] correlation width, C slope]
% - seed: the seed to generate random numbers
% - lambda: wavelength [m]
% - N: size of the arrays
% - s_a: spatial sampling [m]
%
% Outputs:
% - Wr: array NxN including the wave-front perturbation
% - h: array NxN including the surface roughness

    A = Harvey_param(1);
    B = Harvey_param(2);
    C = Harvey_param(3);

    rng(seed)
    h_uc = randn(N,N); % normal noise iid 
    TF_h_uc = fftshift(fft2((h_uc)))/N;
    
    s_f = 1/(N*s_a); % frequency sampling 
    v = ((-1/2):(1/N):(1/2 - 1/N))*N;
    [fx,fy] =  meshgrid(v,v);
    fr = s_f * sqrt(fx.^2 + fy.^2);
    
    Psd = Harvey_PSD(A,B,C,fr); % from Harvey model
  
    G = sqrt(Psd) / s_a; % factor 1/s_a results from discrete representation 
    
    h = fftshift(fft2( TF_h_uc .* fftshift(G) ))/N;
    h = real(h); % remove imaginary residuals
  
    Wr = exp( 2*1i*pi*(n-1)/lambda .* h); 
    
end
 