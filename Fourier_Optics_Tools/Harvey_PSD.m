function [ PSD ] = Harvey_PSD(A,B,C,t)

% Compute the PSD from Harvey et al. model
%
% Inputs:
% - A, B, C: parameters of Harvey model
% - t: matrix or vector with the coordinates (frequency) [m-1]
%
% Outputs:
% - PSD: array or matrix of same size of t

 K = gamma((C+1)/2)/gamma(C/2)  /2/sqrt(pi);
 
 PSD = K * A * B ./(1+(B*t).^2).^((C+1)/2);
 
end

