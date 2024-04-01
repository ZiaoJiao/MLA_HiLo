function[u2]=propTF(u1,L,lambda,z);
% free space propagation - transfer function approach (Fresnel approximation)
% same x and y side length and uniform sampling
% u1 - source plane field
% L - source and observation plane side length
% lambda = wavelength
% z - propagation distance
% u2 - observation plane field

[M,N]=size(u1);                             % get input field array size
dx=L/M;                                     % sample interval
k=2*pi/lambda;                              % wavenumber

fx=-1/(2*dx):1/L:1/(2*dx)-1/L;
fy=fx;                                      % frequency coordinates
[FX,FY]=meshgrid(fx,fy);            

H=exp(j*k*z-j*pi*lambda*z*(FX.^2+FY.^2));   % transfer function
H=fftshift(H);                              % shift transfer function
U1=fft2(fftshift(u1));                      % shift and fourier transfer
                                            % of source field
U2=H.*U1;                                   % spectrum of observation field
u2=ifftshift(ifft2(U2));                    % inverse fourier transfer, and
                                            % center the observation field
end







