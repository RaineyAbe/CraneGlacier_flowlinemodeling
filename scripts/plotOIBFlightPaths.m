% Rainey Aberle
% Fall 2020

% Plot flight paths for all OIB granules in working folder

clear all; close all;

homepath = '/Users/raineyaberle/Desktop/Research/';
addpath([homepath,'general_matlabcodes']);
addpath([homepath,'Tate_RadarInterpretation/IceBridgeCopy/functions']);

% Load Crane centerline
cd([homepath,'Crane_modeling']);
cl.x = load('Crane_centerline.mat').x; 
cl.y = load('Crane_centerline.mat').y;

% Load OIB NC files
cd([homepath,'Tate_RadarInterpretation']);
files = dir('IRMC*');
col = parula(length(files)+1);

figure(1); clf
    plot(cl.x,cl.y,'-m','linewidth',2,'displayname','Crane centerline');
    set(gca,'fontsize',14,'linewidth',2);
    xlabel('Easting (m)'); ylabel('Northing (m)');
    hold on; grid on; legend('Location','west');

for i=1:length(files)
    f = load_L1B(files(i).name);
    [f.Easting,f.Northing] = wgs2ps(f.Longitude,f.Latitude,...
        'StandardParallel',-71,'StandardMeridian',0);
    figure(1); plot(f.Easting,f.Northing,'color',col(i,:),'linewidth',2,...
        'displayname',files(i).name(21:23));
    plot(f.Easting(1),f.Northing(1),'.','color',col(i,:),'markersize',20,...
        'handlevisibility','off');
end
