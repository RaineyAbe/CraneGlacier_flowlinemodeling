% Script to load and display terminus coordinates in shapefiles
% RKA 2021

clear all; close all;

% define homepath
homepath = '/Users/raineyaberle/Desktop/Research/CraneGlacier_flowlinemodeling/';
cd([homepath,'data/terminus/regional/']);

% add path to matlab functions
addpath([homepath,'matlabFunctions/']);

% load Landsat image
cd([homepath,'data/Imagery/']);
landsat = dir('*B8.TIF');
[LS.im,LS.R] = readgeoraster(landsat.name); [LS.ny,LS.nx] = size(LS.im);
% polar stereographic coordinates of image boundaries
LS.x = linspace(min(LS.R.XWorldLimits),max(LS.R.XWorldLimits),LS.nx); 
LS.y = linspace(min(LS.R.YWorldLimits),max(LS.R.YWorldLimits),LS.ny);

% load all terminus coordinate shapefiles
cd([homepath,'data/terminus/regional/']);
names = [{'Crane'},{'Drygalski'},{'Flask'},{'HekGreen'},{'Jorum'}]; % glacier names
% set up figures 
figure(1); clf; hold on; set(gcf,'position',[384 80 999 717]); 
figure(2); clf; hold on; set(gcf,'position',[400 80 999 717]); 
% loop through names
for i=1:length(names)
    cd(char(names(i))); % enter folder 
    files = dir('*.shp'); % load file
    % loop through files to grab centerline
    for j=1:length(files)
        if contains(files(j).name,'CL')
            file = shaperead(files(j).name);
            cl.lon = file.X; % Lon
            cl.lat = file.Y; % Lat
            [cl.X,cl.Y] = wgs2ps(cl.lon,cl.lat,'StandardParallel',-71,'StandardMeridian',0);
            % define x as distance along centerline
            x = zeros(1,length(cl.X));
            for k=2:length(x)
                x(k) = sqrt((cl.X(k)-cl.X(k-1))^2+(cl.Y(k)-cl.Y(k-1))^2)+x(k-1); %(m)
            end
        end
    end
    % initialize f
    clear f; f(1).lon = NaN; f(1).lat = NaN; f(1).date = NaN;
    % loop through files to grab terminus positions
    for j=1:length(files)
        if ~contains(files(j).name,'CL')
            file = shaperead(files(j).name);
            % loop through file
            for k=1:length(file)
                f(length(f)+1).lon = file(k).X; % Lon
                f(length(f)).lat = file(k).Y; % Lat
                % convert to polar stereographic coordinates
                [f(length(f)).X,f(length(f)).Y] = wgs2ps(f(length(f)).lon,f(length(f)).lat,'StandardParallel',-71,'StandardMeridian',0);
                f(length(f)).x = x(dsearchn([cl.X' cl.Y'],[nanmean(f(length(f)).X) nanmean(f(length(f)).Y)])); % x
            end
        end 
    end
    % Plot
    col = parula(length(f)); % color scheme for plotting
    figure(1);
    subplot(2,3,i); hold on; grid on;
        xlabel('Easting (km)'); ylabel('Northing (km)'); title(char(names(i)));
        set(gca,'fontsize',14,'fontname','Arial','linewidth',2);
        for j=1:length(f)
            plot(f(j).X/10^3,f(j).Y/10^3,'color',col(j,:),'linewidth',2);
            drawnow 
        end
    figure(2); hold on;
    subplot(2,3,i); hold on; grid on;
        ylabel('km'); title(char(names(i)));
        set(gca,'fontsize',14,'fontname','Arial','linewidth',2); 
        ct=1;
        for j=2:length(f)
            plot(f(j).x/10^3,'o','color',col(ct,:),'markersize',10,'linewidth',2);
            ct=ct+1;
        end
    cd('../'); % exit folder
end




    


