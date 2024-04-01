function[u2]=propIR(u1,L,lambda,z);
% free space propagation - impulse response approach (Fresnel approximation)
% sampe x and y side length and uniform sampling
% u1 - source plane field
% L - source and observation plane side length
% lambda = wavelength
% z - propagation distance
% u2 - observation plane field

[M,N]=size(u1);                             % get input field array size
dx=L/M;                                     % sample interval
k=2*pi/lambda;                              % wavenumber

x=-L/2:dx:L/2-dx;
y=x;                                        % spatial coordinates
[X,Y]=meshgrid(x,y);            

h=1/(j*lambda*z)*exp(j*k*z+j*k/(2*z)*(X.^2+Y.^2));   % impulse response
H=fft2(fftshift(h))*dx^2;                            % transfer function
U1=fft2(fftshift(u1));                      % shift and fourier transfer
                                            % of source field
U2=H.*U1;                                   % spectrum of observation field
u2=ifftshift(ifft2(U2));                    % inverse fourier transfer, and
                                            % center the observation field
end