clc;
clear;
one=ones(1000,1000);
for i=1:200
    map2(:,:,i)=one;
end
% figure(1)
% imshow(map(:,:,1));

L=40;
K=400;
for i=1:200
    map2(:,:,i)=0;
end

for z=1:200
    for y=1:1000
        for x=1:1000
            if((301<=x)&&(x<=301+K)&&(321<=y)&&(y<=321+L))
                map2(x,y,z)=1;
            end
            
            if((301<=x)&&(x<=301+K)&&(401<=y)&&(y<=401+L))
                map2(x,y,z)=1;
            end
            
            if((301<=x)&&(x<=301+K)&&(481<=y)&&(y<=481+L))
                map2(x,y,z)=1;
            end
            
            if((301<=x)&&(x<=301+K)&&(561<=y)&&(y<=561+L))
                map2(x,y,z)=1;
            end
            
            if((301<x)&&(x<301+K)&&(641<y)&&(y<641+L))
                map2(x,y,z)=1;
            end
        end
    end
end

L=4000; % side length unit is um
M=1000; % # of samples
dx=L/M; % sample intervals
x=-L/2:dx:L/2-dx; y=x; % coord
[X,Y]=meshgrid(x,y);

dz=50;                        % slice interval
vh=10000;                         % volume height
vz=0:dz:vh-dz;            % height coordicate

[X,Y,VZ]=meshgrid(x/10,y/10,vz/100);
xslice = [0];   
yslice = [0];
zslice = [0];
% volume map
figure(1)
slice(X,Y,VZ,map2/max(max(max(map2))),xslice,yslice,zslice)
colormap jet
caxis([0 1])
axis equal; axis xy;
shading interp

%figure(2)
%imshow(map2(:,:,21));

%set(gcf,'unit','millimeters','position',[-0.125 0.125 -0.125 0.125])
% axis(xmin xmax ymin ymax zmin zmax)


