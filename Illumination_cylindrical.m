clear all
close all
clc

%% Sample (MLA imaging)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L=4000; % side length unit is um
M=1000; % # of samples
dx=L/M; % sample intervals
x=-L/2:dx:L/2-dx; y=x; % coord
[X,Y]=meshgrid(x,y);

%% Parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w=25; % x half width (lens radius)
lambda=0.520; % wavelength
lens_pitch=120; % distance bettween neighbour microlens
f=2500; % focal length 
z=2500; % propagation distance
k=2*pi/lambda; % wavenumber
lz=lambda*z;

%% Fresnel propagation
Iin=zeros(size(x));
Ipattern=zeros(size(x));
Ipattern_ana=zeros(size(x));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m=-20:1:20
    uin_loop=rect((X-m*lens_pitch)/w); % Incident field
    Iin_loop=(abs(uin_loop)).^2; % Incident intensity
    uout_mla_loop=focuscyl(uin_loop,L,lambda,f,m*lens_pitch); % Field after microlens array 
    upattern_loop=propTF(uout_mla_loop,L,lambda,z); % Field of illumination pattern at MLA focal plane
    Ipattern_loop=(abs(upattern_loop)).^2; % Intensity of illumination pattern at MLA focal plane
    Iin=Iin+Iin_loop;
    Ipattern=Ipattern+Ipattern_loop;
    %%% Analytical
    Ipattern_loop_ana=(w^2/lz)^2.*(jinc(w/lz*sqrt((X-m*lens_pitch).^2))).^2;
    Ipattern_ana=Ipattern_ana+Ipattern_loop_ana;
end

figure(1);
imagesc(x,y,Iin);   
axis xy; axis square;
colormap('Parula')
xlabel('x (um)'); ylabel('y (um)'); title('Incident Intensity')

figure(2);
imagesc(x,y,Ipattern/max(max(Ipattern)));
colormap Parula
axis xy; axis square;
xlabel('x (um)'); ylabel('y (um)'); title('Illumination Intensity created by MLA')
caxis([0,1])

figure(3);
imagesc(x,y,Ipattern_ana/max(max(Ipattern_ana)));
colormap Parula
axis xy; axis square;
xlabel('x (um)'); ylabel('y (um)'); title('Illumination Intensity created by MLA-Analytical result')
caxis([0,1])



%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NA_obj=0.5;                 % NA of objective lens
f0=NA_obj/lambda;         % coherent cutoff freq

fu=-1/(2*dx):1/L:1/(2*dx)-(1/L);
fv=fu;
[Fu,Fv]=meshgrid(fu,fv);


% %% Incoherent Illumination
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H=circ(sqrt(Fu.^2+Fv.^2)/f0); % coherent transfer function
% OTF=ifft2(abs(fft2(fftshift(H))).^2);
% OTF=abs(OTF/OTF(1,1));
% figure(3)                   % Check OTF
% imagesc(fu,fv,nthroot(fftshift(OTF),2));
% shading interp
% axis xy; axis square;
% ylabel('fu (cycle/um)'); xlabel('fv (cycle/um)');title('OTF')
% 
% 
% G_pattern=fft2(fftshift(Ipattern));        % Fourier transform of upattern
% G_pattern_sample=G_pattern.*OTF;           % Illumination pattern spatial map on sample plane
% I_pattern_sample=ifftshift(ifft2(G_pattern_sample));
% % remove residual imaginary parts, values<0
% I_pattern_sample=real(I_pattern_sample); mask=I_pattern_sample>0; I_pattern_sample=mask.*I_pattern_sample;
% figure(4)                   %illumination pattern
% imagesc(x,y,I_pattern_sample);
% axis xy; axis square;
% ylabel('x (um)'); xlabel('y (um)'); title('Illumination Intensity at focus plane')

%% Volumetric illumination distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dz=5;                        % slice interval
vh=1000;                         % volume height
vz=0:dz:vh-dz;            % height coordicate

I_pattern_volume=zeros(M,M,vh/dz);
I_pattern_volume_ana=zeros(M,M,vh/dz);
index=1

for wd=0:dz:vh-dz
    W=seidel_5(0,0,-(lambda/NA_obj)*Fu,-(lambda/NA_obj)*Fv,wd,0,0,0,0,0);
    H=circ(sqrt(Fu.^2+Fv.^2)/f0).*exp(-j*k*W);
    %I_psf_volume(:,:,index)=abs(ifftshift(ifft2(fftshift(H)))).^2
    OTF=ifft2(abs(fft2(fftshift(H))).^2);
    OTF=abs(OTF/OTF(1,1));
    G_pattern=fft2(fftshift(Ipattern));   % Fourier transform of upattern
    G_pattern_slice=G_pattern.*OTF;      % Illumination pattern spatial map slice
    I_pattern_slice=ifftshift(ifft2(G_pattern_slice));
    I_pattern_volume(:,:,index)=I_pattern_slice;
    % Analytical
    G_pattern_ana=fft2(fftshift(Ipattern_ana));
    G_pattern_slice_ana=G_pattern_ana.*OTF;
    I_pattern_slice_ana=ifftshift(ifft2(G_pattern_slice_ana));
    I_pattern_volume_ana(:,:,index)=I_pattern_slice_ana;
    index=index+1;
    
end

% remove residual imaginary parts, values<0
I_pattern_volume=real(I_pattern_volume); 
mask=I_pattern_volume>0; 
I_pattern_volume=mask.*I_pattern_volume;

I_pattern_volume_ana=real(I_pattern_volume_ana); 
mask=I_pattern_volume_ana>0; 
I_pattern_volume_ana=mask.*I_pattern_volume_ana;

[X,Y,VZ]=meshgrid((-L/2:dx:L/2-dx)/10,(-L/2:dx:L/2-dx)/10,(0:dz:vh-dz)/10);
xslice = [];   
yslice = [0];
zslice = [0];
figure(4)
slice(X,Y,VZ,I_pattern_volume/max(max(max(I_pattern_volume))),xslice,yslice,zslice)
colormap Parula
axis equal; axis xy;
shading interp
xlabel('x (um)'); ylabel('y (um)'); zlabel('z (um)')
title('3D Illumination Intensity - Numerical result')
caxis([0,1])

figure(5)
slice(X,Y,VZ,I_pattern_volume_ana/max(max(max(I_pattern_volume_ana))),xslice,yslice,zslice)
colormap Parula
axis equal; axis xy;
shading interp
xlabel('x (um)'); ylabel('y (um)'); zlabel('z (um)')
title('3D Illumination Intensity - Analytical result')
caxis([0,1])








