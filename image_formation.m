clear all
close all
clc

load D:\BaiduSyncdisk\Microlens_array_HiLo\simulation\revision\I_pattern_volume_squarewave.mat
load D:\BaiduSyncdisk\Microlens_array_HiLo\simulation\revision\sample.mat

%%SAMPLING%%%%%%%%%%
L=4000; % side length unit is um
M=1000; % # of samples
dx=L/M; % sample intervals
x=-L/2:dx:L/2-dx; y=x; % coord
[X,Y]=meshgrid(x,y);


object = I_pattern_volume .* map2
%object = map2

%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NA_obj=0.5;                 % NA of objective lens
lambda = 0.580
f0=NA_obj/lambda;         % coherent cutoff freq
k=2*pi/lambda; % wavenumber

fu=-1/(2*dx):1/L:1/(2*dx)-(1/L);
fv=fu;
[Fu,Fv]=meshgrid(fu,fv);

image = zeros(1000,1000)

dz=5;                        % slice interval
vh=1000;                         % volume height
vz=0:dz:vh-dz;            % height coordicate


%I_psf_volume=zeros(M,M,vh/dz)
index=1

for wd=0:dz:vh-dz
    W=seidel_5(0,0,-(lambda/NA_obj)*Fu,-(lambda/NA_obj)*Fv,wd,0,0,0,0,0);
    H=circ(sqrt(Fu.^2+Fv.^2)/f0).*exp(-j*k*W);
    %I_psf_volume(:,:,index)=abs(ifftshift(ifft2(fftshift(H)))).^2
    OTF=ifft2(abs(fft2(fftshift(H))).^2);
    OTF=abs(OTF/OTF(1,1));
    G_object=fft2(fftshift(object(:,:,index)));   % Fourier transform of object
    G_object=G_object.*OTF;      
    I_image=ifftshift(ifft2(G_object));
    image = image + I_image;    


    index=index+1;
    
end

% remove residual imaginary parts, values<0
image=real(image); 
mask=image>0; 
image=mask.*image;


figure(1);
imagesc(x/10,y/10,image);   
axis xy; axis square;
colormap('jet')
xlabel('x (um)'); ylabel('y (um)'); title('Incident Intensity')
