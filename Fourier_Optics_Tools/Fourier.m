function [F] = Fourier(A, N)

% This function returns Fourier Transform of A
% Matlab 2D FFT routines
% The intensity is conserved (square modulus)
%
% Inputs:
% - A: array NxN
% - N: size of the array A
%
% Outputs:
% - F: Fourier trasnform of A

F = fftshift(fft2(fftshift(A)))/N;
    
end