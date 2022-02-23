% Script to create a timeseries of width-averaged geometry, surface speed, and 
% terminus positions along the Crane Glacier centerline for 1999
% (pre-ice shelf collapse) to 2017 (post-ice shelf collapse)
%
% Rainey Aberle (raineyaberle@u.boisestate.edu)
% Spring 2022
%
% Outline:
%   0. Initial setup
%   1. Load centerline coordinates and width
%   2. Width-averaged and centerline surface elevations (GOTOPO30 & OIB)
%   3. Width-averaged surface velocities (ITS_LIVE & ERS)
%       a. Load velocity datasets, average over glacier width segments
%       b. Create a complete pre-collapse velocity profile
%   4. Terminus positions (Landsat-derived; Dryak and Enderlin, 2020)
%   5. Width-averaged bed elevation profile
%
% To-Do:
%   - edit section 5
% -------------------------------------------------------------------------

% 0. Initial setup

clear all; close all;

% define homepath
homepath = '/Users/raineyaberle/Research/MS/CraneGlacier_flowlinemodeling/';

% add path to required functions
addpath([homepath,'functions/hugheylab-nestedSortStruct'],...
    [homepath,'functions/']);

%% 1. Load centerline coordinates, width, and fjord end

% -----Centerline-----
cl.X = load([homepath,'inputs-outputs/Crane_centerline.mat']).x; 
cl.Y = load([homepath,'inputs-outputs/Crane_centerline.mat']).y;
% Define x as distance along centerline
x = zeros(1,length(cl.X));
for i=2:(length(cl.X))
    x(i)=sqrt((cl.X(i)-cl.X(i-1))^2+(cl.Y(i)-cl.Y(i-1))^2)+x(i-1);
end
% Convert to lon/lat, load geoidheight at each pt along centerline
[cl.lon, cl.lat] = ps2wgs(cl.X, cl.Y, 'StandardParallel', -71, 'StandardMeridian', 0);
h_geoid = geoidheight(cl.lat, cl.lon);

% -----Width-----
width = load([homepath, 'inputs-outputs/calculatedWidth.mat']).width;

% -----Fjord end-----
fjord_end = shaperead([homepath,'data/terminus/fjord_end.shp']);

%% 2. Width-averaged surface elevations (GOTOPO30 & OIB)

k=0; % counter for number of files

% -----GOTOPO30 (~1996)-----
[GT.h, GT.R] = readgeoraster([homepath,'data/surface_elevations/gt30w120s60_aps.tif']);
GT.h(GT.h==-9999) = NaN; % set no-data values to NaN
% extract x and y coordinates
[GT.ny, GT.nx] = size(GT.h); % number of x and y points
GT.X = linspace(GT.R.XWorldLimits(1),GT.R.XWorldLimits(2),GT.nx);
GT.Y = linspace(GT.R.YWorldLimits(1),GT.R.YWorldLimits(2),GT.ny);
% save info in structure
k = k+1; 
h(k).h_centerline = interp2(GT.X,GT.Y,flipud(double(GT.h)),cl.X,cl.Y);
h(k).h_width_ave = zeros(1,length(cl.X)); % initialize speed variables
h(k).date = '1996'; % observation date
h(k).units = "m"; % speed units
h(k).source = "GTOPO30"; % data source
% loop through centerline points
for j=1:length(cl.X)
    % interpolate speed along each width segment
    h(k).h_width_ave(j) = mean(interp2(GT.X,GT.Y,flipud(double(GT.h)),width.segsx(j,:),width.segsy(j,:)),'omitnan'); % [m]
end    
h(k).numPts = length(h(k).h_centerline(~isnan(h(k).h_centerline))); % number of valid data points in centerline profile

% -----ASTER (2001-2002)-----
Afn = dir([homepath,'data/surface_elevations/ASTER/AST14DEM*_aps.tif']);
% loop through files
for i=1:length(Afn)
    A(i).fn = Afn(i).name;
    [A(i).h, A(i).R] = readgeoraster([homepath,'data/surface_elevations/ASTER/',A(i).fn]);
    A(i).h = double(A(i).h);
    A(i).h(A(i).h==-9999) = NaN; % set no data values to NaN
    [A(i).ny, A(i).nx] = size(A(i).h); % number of x and y points
    A(i).X = linspace(A(i).R.XWorldLimits(1),A(i).R.XWorldLimits(2),A(i).nx);
    A(i).Y = linspace(A(i).R.YWorldLimits(1),A(i).R.YWorldLimits(2),A(i).ny);
    
    % save info in structure
    k = k+1;
    h(k).h_centerline = interp2(A(i).X, A(i).Y, flipud(double(A(i).h)), cl.X, cl.Y);
    h(k).h_width_ave = zeros(1,length(cl.X)); % initialize speed variables
    h(k).date = A(i).fn(17:26); % observation date
    h(k).units = "m"; % elevation units
    h(k).source = "ASTER"; % data source
    % loop through centerline points
    for j=1:length(width.W)
        % interpolate speed along each width segment
        h(k).h_width_ave(j) = mean(interp2(A(i).X, A(i).Y, flipud(A(i).h), width.segsx(j,:), width.segsy(j,:)),'omitnan');
    end    
    h(k).numPts = length(h(k).h_width_ave(~isnan(h(k).h_width_ave))); % number of data points

end

% -----OIB (2009-12, 2016-18)-----
OIB_files = dir([homepath,'data/surface_elevations/OIB_L2/C*.csv']);
% loop through files
for i=1:length(OIB_files)
    
    % load file
    OIB(i).fn = OIB_files(i).name; % file name
    file = readmatrix([homepath,'data/surface_elevations/OIB_L2/',OIB(i).fn]); % full csv file
    file(file==-9999) = NaN; % replace no data values with NaN
    % extract coordinates
    OIB(i).lat = file(:,1);
    OIB(i).lon = file(:,2);
    % convert to Antarctic polar stereographic 
    [OIB(i).X, OIB(i).Y] = wgs2ps(OIB(i).lon, OIB(i).lat, 'StandardParallel',-71,'StandardMeridian',0);
    
    % extract surface elevations closest to each centerline point
    maxDist = 1500; % maximum distance from centerline
    OIB(i).h = zeros(1,length(cl.X)); % initialize variable
    % loop through centerline points
    for j=1:length(cl.X)
        I = dsearchn([OIB(i).X OIB(i).Y],[cl.X(j) cl.Y(j)]); % index of closest point
        if pdist([OIB(i).X(I) OIB(i).Y(I); cl.X(j) cl.Y(j)],'euclidean')<=maxDist
            OIB(i).h(j) = file(I,5)-file(I,7)-h_geoid(j);
        else
            OIB(i).h(j)= NaN;
        end
    end

    % Save in h structure
    k = k+1;
    h(k).h_centerline = OIB(i).h;
    h(k).h_width_ave = NaN;
    h(k).date = OIB(i).fn(18:25);
    h(k).units = 'm'; 
    h(k).source = 'OIB'; 
    h(k).numPts = length(h(k).h_centerline(~isnan(h(k).h_centerline))); % number of valid data points in centerline profile
    h(k).h_width_ave_adj = NaN;
    
end

% -----remove outliers-----
max_slope = 0.2; % [m/m]
% loop through h profiles
for i=1:length(h)
    % if not from OIB
    if ~contains(h(i).source, 'OIB')
        % initialize variables
        h(i).slope_width_ave = zeros(1,length(x));
        h(i).h_width_ave_adj = zeros(1,length(x));
        % calculate slope
        h(i).slope_width_ave(1) = (h(i).h_width_ave(2) - h(i).h_width_ave(1)) / (x(2) - x(1)); % forward difference
        h(i).slope_width_ave(2:end-1) = (h(i).h_width_ave(3:end) - h(i).h_width_ave(1:end-2)) ./ (x(3:end) - x(1:end-2)); % central difference
        h(i).slope_width_ave(end) = (h(i).h_width_ave(end) - h(i).h_width_ave(end-1)) / (x(end) - x(end-1)); % backward difference
        % loop through centerline points
        for j=1:length(cl.X)
            if abs(h(i).slope_width_ave(j)) < max_slope && h(i).h_width_ave(j) > 0
                h(i).h_width_ave_adj(j) = h(i).h_width_ave(j);
            else 
                h(i).h_width_ave_adj(j) = NaN;
            end
        end
    else
        h(i).slope_width_ave = NaN;
        h(i).h_width_ave_adj = NaN;
    end
end

% -----calculate h spatial mean for pre-collapse profiles----- 
% loop through h profiles
for i=1:length(h)
    % if not during or after 2002
    if ~contains(h(i).date, '201') && ~contains(h(i).date, '2002')
        % add data values to vector
        h_med(i, :) = h(i).h_width_ave_adj;
    end
end
% calculate spatial median
h_med = median(h_med, 1, 'omitnan');
% apply median filter
h_med_medfilt = medfilt1(h_med, 15, 'omitnan');
h_med_medfilt(1) = h_med_medfilt(2);
% add to h variable
k = k+1;
h(k).h_centerline = NaN;
h(k).h_width_ave = h_med_medfilt;
h(k).date = 'Pre-collapse median profile';
h(k).units = 'm';
h(k).source = 'ASTER, GTOPO30';
h(k).numPts = length(h(k).h_width_ave(~isnan(h(k).h_width_ave)));
% calculate slope
h(k).slope_width_ave = zeros(1,length(x));
h(k).slope_width_ave(1) = (h(k).h_width_ave(2) - h(k).h_width_ave(1)) / (x(2) - x(1)); % forward difference
h(k).slope_width_ave(2:end-1) = (h(k).h_width_ave(3:end) - h(k).h_width_ave(1:end-2)) ./ (x(3:end) - x(1:end-2)); % central difference
h(k).slope_width_ave(end) = (h(k).h_width_ave(end) - h(k).h_width_ave(end-1)) / (x(end) - x(end-1)); % backward difference

% -----plot-----
figure(1); clf;
set(gcf,'position',[200 70 1050 650]);
col = parula(length(h)+1); % color scheme for plotting lines
% centerline
ax1 = subplot(2,2,1);
hold on; grid on; legend;
set(gca,'fontsize',12,'linewidth',1);
ylabel('elevation [m]');
title('h_{centerline}');
% width-averaged 
ax2 = subplot(2,2,2);
hold on; grid on; 
set(gca,'fontsize',12,'linewidth',1);
title('h_{width-averaged}');
% slope
ax3 = subplot(2,2,3);
hold on; grid on;
set(gca,'fontsize',12,'linewidth',1);
xlabel('distance along centerline [km]');
ylabel('[m/m]');
title('| dh/dx |');
% adjusted surface
ax4 = subplot(2,2,4); 
hold on; grid on;
set(gca,'fontsize',12,'linewidth',1);
xlabel('distance along centerline [km]');
ylabel('elevation [m]');
title('h_{width-averaged, adjusted}');
for i=1:length(h)-1
    plot(ax1, x/10^3, h(i).h_centerline, 'linewidth', 1, 'displayname', num2str(h(i).date(1:4)), 'color', col(i,:));
    if ~isnan(h(i).h_width_ave)
        plot(ax2, x/10^3, h(i).h_width_ave, 'linewidth', 1, 'color',col(i,:));
    end
    plot(ax3, x/10^3, abs(h(i).slope_width_ave), '.', 'markersize', 1, 'color',col(i,:));
    plot(ax3, x/10^3, max_slope*ones(1, length(x)), '--k', 'linewidth', 2);
    plot(ax4, x/10^3, h(i).h_width_ave_adj, 'linewidth', 1, 'color',col(i,:));
end 
plot(ax4, x/10^3, h_med_medfilt, '-k', 'linewidth', 2);

% -----save h and figure-----
save([homepath,'inputs-outputs/surfaceElevationObs_1996-2018.mat'], 'h');
disp('h saved to file');
saveas(figure(1), [homepath,'figures/surfaceElevation.png'], 'png');
disp('figure saved to file');

%% 3. Width-averaged surface velocities
% -------------------------------------------------------------------------
%   a. Load velocity datasets, average over glacier width segments
% -------------------------------------------------------------------------

k=0; % counter for files

% -----ERS (1994)-----
k=k+1;
% file name
ERS_files{1} = [homepath,'data/surface_velocities/ENVEO_velocities/LarsenFleming_s19940120_e19940319.1.0_20170928/LarsenFleming_s19940120_e19940319.tif'];
% load file
[ERS(1).A, ERS(1).R] = readgeoraster(ERS_files{1});
[ERS(1).ny, ERS(1).nx, ~] = size(ERS(1).A); % dimension sizes
ERS(1).X = linspace(ERS(1).R.XWorldLimits(1),ERS(1).R.XWorldLimits(2),ERS(1).nx); % X [m]
ERS(1).Y = linspace(ERS(k).R.YWorldLimits(1),ERS(1).R.YWorldLimits(2),ERS(1).ny); % Y [m]
ERS(1).ux = ERS(1).A(:,:,1); ERS(1).ux(ERS(1).ux==single(1e20)) = NaN; % Easting velocity [m/d]
ERS(1).uy = ERS(1).A(:,:,2); ERS(1).uy(ERS(1).uy==single(1e20)) = NaN; % Northing velocity [m/d]
ERS(1).u = sqrt(ERS(1).ux.^2 + ERS(1).uy.^2); % velocity magnitude [m/d]
% save info in structure
U(k).date = '1994'; % observation date
U(k).units = "m/s"; % speed units
U(k).source = "ERS"; % data source
% interpolate speed along centerline
U(k).U_centerline = interp2(ERS(1).X,ERS(1).Y,flipud(ERS(1).u),cl.X,cl.Y)./24/60/60; % [m/s];
% initialize width-averaged speed 
U(k).U_width_ave = zeros(1,length(cl.X)); 
% loop through centerline points
for j=1:length(cl.X)
    % interpolate speed along each width segment
    U(k).U_width_ave(j) = mean(interp2(ERS(1).X,ERS(1).Y,flipud(ERS(1).u),width.segsx(j,:),width.segsy(j,:)),'omitnan')./24/60/60; % [m/s]
end    
U(k).numPts = length(U(k).U_width_ave(~isnan(U(1).U_width_ave))); % number of valid data points in centerline profile

% -----ERS (1995)-----
k=k+1;
% file name
ERS_files{2} = [homepath,'data/surface_velocities/ENVEO_velocities/glacapi_iv_LB_ERS_1995_v2.tif'];
% load file
[ERS(2).A, ERS(2).R] = readgeoraster(ERS_files{2});
[ERS(2).ny, ERS(2).nx, ~] = size(ERS(2).A); % dimension sizes
ERS(2).X = linspace(ERS(2).R.XWorldLimits(1),ERS(2).R.XWorldLimits(2),ERS(2).nx); % X [m]
ERS(2).Y = linspace(ERS(2).R.YWorldLimits(1),ERS(2).R.YWorldLimits(2),ERS(2).ny); % Y [m]
ERS(2).ux = ERS(2).A(:,:,1); ERS(2).ux(ERS(2).ux==single(3.4028235e+38)) = NaN; % Easting velocity [m/d]
ERS(2).uy = ERS(2).A(:,:,2); ERS(2).uy(ERS(2).uy==single(3.4028235e+38)) = NaN; % Northing velocity [m/d]
ERS(2).u = sqrt(ERS(2).ux.^2 + ERS(2).uy.^2); % velocity magnitude [m/d]
% save info in structure
U(k).date = '1995'; % observation date
U(k).units = "m/s"; % speed units
U(k).source = "ERS"; % data source
% interpolate speed along centerline
U(k).U_centerline = interp2(ERS(2).X,ERS(2).Y,flipud(ERS(2).u),cl.X,cl.Y)./24/60/60; % [m/s];
% initialize width-averaged speed
U(k).U_width_ave = zeros(1,length(cl.X)); 
% loop through centerline points
for j=1:length(cl.X)
    % interpolate speed along each width segment
    U(k).U_width_ave(j) = mean(interp2(ERS(2).X,ERS(2).Y,flipud(ERS(2).u),width.segsx(j,:),width.segsy(j,:)),'omitnan')./24/60/60; % [m/s]
end    
U(k).numPts = length(U(k).U_width_ave(~isnan(U(2).U_width_ave))); % number of valid data points in width-averaged profile

% -----ITS_LIVE (1999-2017)-----
% file names
ITS_LIVE_files = dir([homepath,'data/surface_velocities/ITS_LIVE/ANT*.nc']);
% Loop through files
for i=1:length(ITS_LIVE_files)
    % load coordinates
    IL(i).X = ncread([homepath,'data/surface_velocities/ITS_LIVE/',ITS_LIVE_files(i).name],'x'); % X [m]
    IL(i).Y = ncread([homepath,'data/surface_velocities/ITS_LIVE/',ITS_LIVE_files(i).name],'y'); % Y [m]
    % determine subset over which to extract velocity data
    Ix = find(IL(i).X >= min(width.segsx(:)) & IL(i).X <= max(width.segsx(:)));
    Iy = find(IL(i).Y >= min(width.segsy(:)) & IL(i).Y <= max(width.segsy(:)));
    startloc = [Ix(1) Iy(1)]; counts = [range(Ix)+1 range(Iy)+1];
    % crop coordinates to subset
    IL(i).X = IL(i).X(Ix); IL(i).Y = IL(i).Y(Iy);
    % read in velocity
    IL(i).U = ncread([homepath,'data/surface_velocities/ITS_LIVE/',ITS_LIVE_files(i).name],'v', startloc, counts)'; % [m/y]
    IL(i).U(IL(i).U==-32767) = NaN; % replace no data values with NaN
    % convert to m/s
    IL(i).U = IL(i).U./3.1536e7;
    % read in velocity error
    IL(i).U_err = ((ncread([homepath,'data/surface_velocities/ITS_LIVE/',ITS_LIVE_files(i).name],'v_err', startloc, counts))')./3.1536e7; % u error [m/s]
    % save info in structure
    k=k+1;
    U(k).date = ITS_LIVE_files(i).name(11:14); % observation date
    U(k).units = "m/s"; % speed units
    U(k).source = "ITS-LIVE"; % speed data source
    % interpolate speed along centerline
    U(k).U_centerline = interp2(IL(i).X, IL(i).Y, IL(i).U, cl.X, cl.Y); % [m/s];
    % initialize width-averaged speed
    U(k).U_width_ave = zeros(1,length(cl.X));
    % loop through centerline points
    for j=1:length(width.W)
        % interpolate speed and average along each width segment
        U(k).U_width_ave(j) = mean(interp2(IL(i).X, IL(i).Y, IL(i).U, width.segsx(j,:), width.segsy(j,:)),'omitnan');
        U(k).U_err(j) = mean(interp2(IL(i).X ,IL(i).Y, IL(i).U_err, width.segsx(j,:), width.segsy(j,:)),'omitnan');
    end    
    U(k).numPts = length(U(k).U_width_ave(~isnan(U(k).U_width_ave))); % number of data points
end

% -----plot-----
col = parula(length(U)+1); % color scheme for plotting lines
figure(2); clf; 
hold on; grid on; legend('location','east outside');
set(gca,'fontsize',12,'linewidth',1);
xlabel('distance along centerline [km]');
ylabel('speed [m/yr]');
for i=1:length(U)
    plot(x/10^3,U(i).U_width_ave*3.1536e7,'color',col(i,:),'displayname',num2str(U(i).date),'linewidth',1);
end

% -------------------------------------------------------------------------
%   b. Create a complete pre-collapse velocity profile
% -------------------------------------------------------------------------

% intersection of centerline with fjord_end
% I_end = 184;

% grab velocity profiles for pre-collapse (t1) and post-collapse (t2)
Ut1 = U(1).U_width_ave; % 1994
Ut1(find(~isnan(U(2).U_width_ave),1,'first'):end) = U(2).U_width_ave(find(~isnan(U(2).U_width_ave),1,'first'):end); % 1995
Ut2 = U(21).U_width_ave;

% find intersection points of Ut1 and Ut2
I1 = find(x > 23.4e3, 1, 'first');
I2 = find(x > 51.9e3, 1, 'first');

% normalize t2 profile from 0 to 1
Ut2_norm = normalize(Ut2,'range');

% rescale to t1 range
Ut1_fill = rescale(Ut2_norm, Ut1(1), Ut1(I2));
% use pre-collapse observed speed before gap
Ut1_fill(1:I1) = Ut1(1:I1);
% % use mean observed speed after gap
% Ut1_fill(I2:end) = mean([Ut1(I2:end); Ut2(I2:end)],1);
% % apply median filter from I1:end
% Ut1_fill(I1:end) = medfilt1(Ut1_fill(I1:end),10);
% % use mean speeds before intersection
% Ut1_fill(1:I1) = mean([Ut1(1:I1); Ut2(1:I1)],1);

% -----plot-----
figure(3); clf; hold on;
set(gca,'fontsize',12,'linewidth',1);
plot(x/10^3, Ut1*3.1536e7, 'linewidth', 2, 'displayname', 't_1');
plot(x/10^3, Ut2*3.1536e7, 'linewidth', 2, 'displayname', 't_2');
plot(x/10^3, Ut1_fill*3.1536e7, '--', 'linewidth', 2, 'displayname', 't_1 filled');
xlabel('distance along centerline [km]');
ylabel('speed [m/yr]');
grid on; legend('location','northwest');

% -----save-----
% save info in structure
k=k+1;
U(k).U_width_ave = Ut1_fill; 
U(k).date = 'pre-collapse';
U(k).units = 'm/s';
U(k).source = 'ERS, ITS_LIVE';
U(k).U_centerline = NaN;
U(k).numPts = NaN;
U(k).U_err = NaN;
% save U variable
save([homepath,'inputs-outputs/speedsWidthAveraged.mat'],'U');
disp('U saved to file');
% save figures
saveas(figure(2), [homepath,'figures/centerlineSpeeds.png'], 'png');
disp('figure 2 saved');
saveas(figure(3), [homepath,'figures/pre-collapse_speed.png'], 'png');
disp('figure 3 saved');

%% 4. Terminus positions (Landsat-derived; Dryak and Enderlin, 2020)

% -----load files-----
term.X = load([homepath,'inputs-outputs/LarsenB_centerline.mat']).centerline.termx;
term.Y = load([homepath,'inputs-outputs/LarsenB_centerline.mat']).centerline.termy;
term.x = x(dsearchn([cl.X cl.Y],[term.X' term.Y']));
term.date = load([homepath,'inputs-outputs/LarsenB_centerline.mat']).centerline.termdate;

% -----plot-----
figure(4); clf; 
grid on; hold on; 
set(gca,'fontsize',12,'linewidth',2);
plot(term.date,term.x/10^3,'.-b','markersize',15);
xlabel('date'); 
ylabel('distance along centerline [km]');

%% 5. Width-averaged bed elevation profile

% -------------------------------------------------------------------------
%   a. Load bed elevation (OIB picks) and bathymetry (Rebesco, 2006)
% -------------------------------------------------------------------------

b = load([homepath,'inputs-outputs/observedBed.mat']).HB.hb0;
bathym = load([homepath,'inputs-outputs/observedBed.mat']).bathym;

% -------------------------------------------------------------------------
%   b. Average thickness over width to adjust bed profile assuming a
%      parabolic bed shape
% -------------------------------------------------------------------------

% load one surface elevation profile
h_2018 = load([homepath,'inputs-outputs/surfaceElevationObs.mat']).h(1).surface;
h_2018(find(isnan(h_2018),1,'first'):end) = 0;

% calculate thickness
H = h_2018-hb;

% calculate parabolic cross-sectional area 
% y = a(x-h)^2+k --> Hn = a(x-(W/2))^2-H, where a = 4H/(W^2)
Hn = NaN*zeros(length(W),round(max(W/10)));
for i=1:length(W)
    xi = 0:10:W(i);
    a = 4*H(i)/(W(i)^2);
    Hn(i,1:length(xi)) = -(a.*(xi-W(i)/2).^2-H(i));
end

% adjust H to average over width
% Area of an ellipse 
H_adj = nanmean(Hn,2)'; 

% adjust hb using thickness and surface
hb_adj = h_2018-H_adj;

% plot results
figure(6); clf; hold on; legend;
set(gcf,'position',[200 200 1000 800]); title('Width-Averaged Bed');
xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)');
set(gca,'linewidth',2,'fontsize',18); grid on;
plot(x/10^3,hb,'--k','linewidth',2,'DisplayName','hb_{cl}');
plot(x/10^3,hb_adj,'-k','linewidth',2,'DisplayName','hb_{adj}');
plot(x/10^3,h_2018,'-b','linewidth',2,'DisplayName','h');

% save results
if save_hb_adj
    save('delineatedBedWidthAveraged.mat','hb_adj','x');
    disp('hb_adj saved');
end
