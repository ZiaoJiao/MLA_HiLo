function[uout]=focuscyl(uin,L,lambda,zf,pitchx)
% converging or diverging phase-front
% uniform sampling assume
% uin - input field
% L - side length
% lambda - wavelength
% zf - focal distance (+ converge, - diverge)
% uout - output field
% pitchx and pitchy - lens pitch

[M,N]=size(uin);
dx=L/M;
k=2*pi/lambda;

x=-L/2:dx:L/2-dx;
y=x;
[X,Y]=meshgrid(x,y);

uout=uin.*exp(-j*k/(2*zf)*((X-pitchx).^2));
end