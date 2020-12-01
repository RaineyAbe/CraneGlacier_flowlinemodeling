%Script to grab terminus info from shapefiles
%Rainey Aberle
%February 2020

figure;
glacierName = 'Crane Glacier';

%Display Landsat image
    disp('Loading Landsat image... ');
    cd /Users/raineyaberle/Desktop/Research/Crane_modeling/BedandSurface/OIB
    landsat = dir('*B8.TIF');
    [LS,R] = geotiffread(landsat.name);
    [ny,nx] = size(LS);

    %Polar stereographic coordinates of image boundaries (from image metadata)
    x = linspace(-2.500500e6,-2.2219700e6,nx);
    y = linspace(1.154400e6, 1.434600e6,ny);

    imagesc(x,y,flipud(LS)); hold on;

set(gca,'FontName','Times New Roman','FontSize',18,'YDir','normal'); 
title([glacierName,' Terminus Positions']);
xlabel('Easting (m)'); ylabel('Northing (m)'); grid on; axis equal;

%pre-2016
    disp('Loading Terminus Data...');
    cd /Users/raineyaberle/Desktop/Research/Crane_modeling/terminus/pre-2016
    file = dir('*.shp');
    pre16 = shaperead(file.name); dates_pre16 = extractfield(pre16,'Date');
    mapshow(pre16,'Color','b','LineWidth',1); 

%2016-2018
    cd /Users/raineyaberle/Desktop/Research/Crane_modeling/terminus/2016-2018
    file = dir('*.shp'); 
    post16 = shaperead(file.name); dates_post16 = extractfield(post16,'Date');
    mapshow(post16,'Color','r','LineWidth',1); 

%Zoom in
xlim([-2.418e6 -2.396e6]); ylim([1.258e6 1.278e6]);
%%
%Grab terminus positions along centerline

    %Define the centerline
    cd /Users/raineyaberle/Desktop/Research/Crane_modeling
    load('LarsenB_centerline.mat');
    cl_x = centerline.x; cl_y = centerline.y;
    
    plot(cl_x,cl_y,'-g','DisplayName','Approximate Flowline','LineWidth',1); 
    legend; 
    hold off;
    
    %Grab points where termini intersect flowline
    term_coords = ([extractfield(pre16,'X'); extractfield(pre16,'Y')]);
        %This returns ALL x and y spliced together
        %How to separate them and grab a point?
        %For each terminus line, interp terminus location at center
        %flowline & grab date
    
    
    
    


