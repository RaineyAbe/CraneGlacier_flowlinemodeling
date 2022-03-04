% Script to create timeseries of width-averaged geometry, surface speed, and 
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
%   4. Glacier bed elevation (OIB picks)
%   5. Width-averaged thicknesses
%   6. Terminus positions (Landsat-derived; Dryak and Enderlin, 2020)
%
% To-Do:
%   - surface: ASTER
%   - only read portion of ITS_LIVE to speed up process 
%   - smooth pre-collapse velocity profile at transitions
%   - edit sections 4+
% -------------------------------------------------------------------------

% 0. Initial setup

clear all; close all;

% define homepath
homepath = '/Users/raineyaberle/Research/MS/CraneGlacier_flowlinemodeling/';

% add path to required functions
addpath([homepath,'functions/hugheylab-nestedSortStruct'], [homepath,'functions/']);

%% 1. Load centerline coordinates and width

% -----Centerline-----
cl.x = load([homepath,'inputs-outputs/Crane_centerline.mat']).x; 
cl.y = load([homepath,'inputs-outputs/Crane_centerline.mat']).y;
% Define x as distance along centerline
x = zeros(1,length(cl.x));
for i=2:(length(cl.x))
    x(i)=sqrt((cl.x(i)-cl.x(i-1))^2+(cl.y(i)-cl.y(i-1))^2)+x(i-1);
end
% Convert to lon/lat, load geoidheight at each pt along centerline
[cl.lon, cl.lat] = ps2wgs(cl.x, cl.y, 'StandardParallel', -71, 'StandardMeridian', 0);
h_geoid = geoidheight(cl.lat, cl.lon);

% -----Width-----
width = load([homepath, 'inputs-outputs/calculatedWidth.mat']).width;

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
h(k).h_centerline = interp2(GT.X,GT.Y,flipud(double(GT.h)),cl.x,cl.y);
h(k).h_width_ave = zeros(1,length(cl.x)); % initialize speed variables
h(k).date = '1996'; % observation date
h(k).units = "m"; % speed units
h(k).source = "GTOPO30"; % data source
% loop through centerline points
for j=1:length(cl.x)
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
    h(k).h_centerline = interp2(A(i).X, A(i).Y, flipud(double(A(i).h)), cl.x, cl.y);
    h(k).h_width_ave = zeros(1,length(cl.x)); % initialize speed variables
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
    OIB(i).h = zeros(1,length(cl.x)); % initialize variable
    % loop through centerline points
    for j=1:length(cl.x)
        I = dsearchn([OIB(i).X OIB(i).Y],[cl.x(j) cl.y(j)]); % index of closest point
        if pdist([OIB(i).X(I) OIB(i).Y(I); cl.x(j) cl.y(j)],'euclidean')<=maxDist
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
        for j=1:length(cl.x)
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

% -----ERS (1994)-----
[ERS(1).A, ERS(1).R] = readgeoraster([homepath,'data/surface_velocities/ENVEO_velocities/LarsenFleming_s19940120_e19940319.1.0_20170928/LarsenFleming_s19940120_e19940319.tif']);
[ERS(1).ny, ERS(1).nx, ~] = size(ERS(1).A); % dimension sizes
ERS(1).X = linspace(ERS(1).R.XWorldLimits(1),ERS(1).R.XWorldLimits(2),ERS(1).nx); % X [m]
ERS(1).Y = linspace(ERS(1).R.YWorldLimits(1),ERS(1).R.YWorldLimits(2),ERS(1).ny); % Y [m]
ERS(1).ux = ERS(1).A(:,:,1); ERS(1).ux(ERS(1).ux==single(1e20)) = NaN; % Easting velocity [m/d]
ERS(1).uy = ERS(1).A(:,:,2); ERS(1).uy(ERS(1).uy==single(1e20)) = NaN; % Northing velocity [m/d]
ERS(1).u = sqrt(ERS(1).ux.^2 + ERS(1).uy.^2); % velocity magnitude [m/d]
% save info in structure
U(1).date = '1994'; % observation date
U(1).units = "m/s"; % speed units
U(1).source = "ERS"; % data source
% interpolate speed along centerline
U(1).U_centerline = interp2(ERS(1).X,ERS(1).Y,flipud(ERS(1).u),cl.x,cl.y)./24/60/60; % [m/s];
% initialize width-averaged speed 
U(1).U_width_ave = zeros(1,length(cl.x)); 
% loop through centerline points
for j=1:length(cl.x)
    % interpolate speed along each width segment
    U(1).U_width_ave(j) = mean(interp2(ERS(1).X,ERS(1).Y,flipud(ERS(1).u),width.segsx(j,:),width.segsy(j,:)),'omitnan')./24/60/60; % [m/s]
end    
U(1).numPts = length(U(1).U_width_ave(~isnan(U(1).U_width_ave))); % number of valid data points in centerline profile

% -----ERS (1995)-----
[ERS(2).A, ERS(2).R] = readgeoraster([homepath,'data/surface_velocities/ENVEO_velocities/glacapi_iv_LB_ERS_1995_v2.tif']);
[ERS(2).ny, ERS(2).nx, ~] = size(ERS(2).A); % dimension sizes
ERS(2).X = linspace(ERS(2).R.XWorldLimits(1),ERS(2).R.XWorldLimits(2),ERS(2).nx); % X [m]
ERS(2).Y = linspace(ERS(2).R.YWorldLimits(1),ERS(2).R.YWorldLimits(2),ERS(2).ny); % Y [m]
ERS(2).ux = ERS(2).A(:,:,1); ERS(2).ux(ERS(2).ux==single(3.4028235e+38)) = NaN; % Easting velocity [m/d]
ERS(2).uy = ERS(2).A(:,:,2); ERS(2).uy(ERS(2).uy==single(3.4028235e+38)) = NaN; % Northing velocity [m/d]
ERS(2).u = sqrt(ERS(2).ux.^2 + ERS(2).uy.^2); % velocity magnitude [m/d]
% save info in structure
U(2).date = '1995'; % observation date
U(2).units = "m/s"; % speed units
U(2).source = "ERS"; % data source
% interpolate speed along centerline
U(2).U_centerline = interp2(ERS(2).X,ERS(2).Y,flipud(ERS(2).u),cl.x,cl.y)./24/60/60; % [m/s];
% initialize width-averaged speed
U(2).U_width_ave = zeros(1,length(cl.x)); 
% loop through centerline points
for j=1:length(cl.x)
    % interpolate speed along each width segment
    U(2).U_width_ave(j) = mean(interp2(ERS(2).X,ERS(2).Y,flipud(ERS(2).u),width.segsx(j,:),width.segsy(j,:)),'omitnan')./24/60/60; % [m/s]
end    
U(2).numPts = length(U(2).U_width_ave(~isnan(U(2).U_width_ave))); % number of valid data points in width-averaged profile

% -----ITS_LIVE (1999-2018)-----
% grab file names
ILfiles = dir([homepath,'data/surface_velocities/ITS_LIVE/ANT*.nc']);
% Loop through files
for i=1:length(ILfiles)

    X = ncread([homepath,'data/surface_velocities/ITS_LIVE/',ILfiles(i).name],'x'); % X [m]
    Y = ncread([homepath,'data/surface_velocities/ITS_LIVE/',ILfiles(i).name],'y'); % Y [m]
    u = (ncread([homepath,'data/surface_velocities/ITS_LIVE/',ILfiles(i).name],'v')')./3.1536e7; % u [m/s]
    u(u==-32767) = NaN; %Replace no data values with NaN
    u_err = ((ncread([homepath,'data/surface_velocities/ITS_LIVE/',ILfiles(i).name],'v_err'))')./3.1536e7; % u error [m/s]
    % save info in structure
    U(i+2).date = ILfiles(i).name(11:14); % observation date
    U(i+2).units = "m/s"; % speed units
    U(i+2).source = "ITS-LIVE"; % speed data source
    % interpolate speed along centerline
    U(i+2).U_centerline = interp2(X,Y,u,cl.x,cl.y); % [m/s];
    % initialize width-averaged speed
    U(i+2).U_width_ave = zeros(1,length(cl.x));
    % loop through centerline points
    for j=1:length(width.W)
        % interpolate speed along each width segment
        U(i+2).U_width_ave(j) = mean(interp2(X,Y,u,width.segsx(j,:),width.segsy(j,:)),'omitnan');
        U(i+2).U_err(j) = mean(interp2(X,Y,u_err,width.segsx(j,:),width.segsy(j,:)),'omitnan');
    end    
    U(i+2).numPts = length(U(i+2).U_width_ave(~isnan(U(i+2).U_width_ave)));  % number of data points
end

% -----plot-----
col = parula(length(U)+1); % color palette for plotting lines
figure(2); clf; hold on;
set(gca,'fontsize',12,'linewidth',1);
legend('location','west');
xlabel('distance along centerline [km]');
ylabel('speed [m/yr]')
for i=1:length(U)
    plot(x/10^3,U(i).U_width_ave*3.1536e7,'color',col(i,:),'displayname',num2str(U(i).date),'linewidth',1);
end

% -----save----- 
save([homepath,'inputs-outputs/speedsWidthAveraged.mat'],'U');
disp('U saved to file');

% -------------------------------------------------------------------------
%   b. Create a complete pre-collapse velocity profile
% -------------------------------------------------------------------------

% grab velocity profiles for pre-collapse (t1) and post-collapse (t2)
Ut1 = U(1).U_width_ave; % 1994
Ut1(find(~isnan(U(2).U_width_ave),1,'first'):end) = U(2).U_width_ave(find(~isnan(U(2).U_width_ave),1,'first'):end); % 1995
Ut2 = U(end).U_width_ave;

% normalize t2 profile from 0 to 1
Ut2_norm = normalize(Ut2,'range');

% rescale to t1 range
Ut1_fill = rescale(Ut2_norm, Ut1(1), Ut1(end));

% -----plot-----
figure(3); clf; hold on;
set(gca,'fontsize',12,'linewidth',1);
plot(x/10^3, Ut1*3.1536e7, '-b', 'displayname', 't_1', 'linewidth', 2);
plot(x/10^3, Ut1_fill*3.1536e7, '--b', 'displayname', 't_1 filled', 'linewidth', 2);
plot(x/10^3, Ut2*3.1536e7, '-m', 'displayname', 't_2', 'linewidth', 2);
xlabel('distance along centerline [km]');
ylabel('speed [m/yr]');
grid on; legend;

% -----save (append)-----
i = length(U)+1;
U(i).U_width_ave = Ut1_fill; 
U(i).date = 'pre-collapse';
U(i).units = 'm/s';
U(i).source = 'ERS/ITS_LIVE';
U(i).U_centerline = NaN;
U(i).numPts = NaN;
U(i).U_err = NaN;
save([homepath,'inputs-outputs/speedsWidthAveraged.mat'],'U','-append');
disp('U appended and saved to file');

%% 2. Bed elevation (OIB)

if bedElev_save
    
    % Grab Tate's observed bed & surface, convert coordinates,
    % interpolate to centerline
    
    % Load Tate's OIB file
    obs = load('OIBPicks_2018_005.mat').CraneOIBPicks_2018;
    b0 = obs.hb_geoid;
    x_b0 = obs.Easting;
    y_b0 = obs.Northing;
    
    % Interpolate bed elevation for each point along the centerline
    % within a certain distance
    maxDist = 5e3; % m
    b = zeros(1,length(cl.X));
    for i=1:length(cl.X)
        Ib = dsearchn([x_b0 y_b0],[cl.X(i) cl.Y(i)]);
        dist = sqrt((x_b0(Ib)-cl.X(i))^2+(y_b0(Ib)-cl.Y(i))^2);
        if dist<=maxDist
            b(i) = b0(Ib);
        else
            b(i)=NaN;
        end
    end
    % set end values as constant
    b(170:end) = b(170);
        
    % save b and coordinates in structure
    b.b0 = b; B.X = cl.X'; B.Y = cl.Y';
    
    % load bathymetry observations
    bathym = load('bathymetryData.mat').cl_trough;
    bathym = -1.*bathym; % Data reported in depth beneath the surface (convert to negative values for plotting)
    
    % save bed elevations
    save([homepath,'inputs-outputs/observedBed.mat'],'B','bathym');
    disp('bed elevations saved');

else
    
    % load bed elevations
    hb = load([homepath,'inputs-outputs/observedBed.mat']).B.b0;
    bathym = load([homepath,'inputs-outputs/observedBed.mat']).bathym; 
    
end

%% 3. Terminus positions (Landsat-derived; Dryak and Enderlin, 2020)

termX = load([homepath,'inputs-outputs/LarsenB_centerline.mat']).centerline.termx;
termY = load([homepath,'inputs-outputs/LarsenB_centerline.mat']).centerline.termy;
termx = x(dsearchn([cl.X cl.Y],[termX' termY']));
termdate = load([homepath,'inputs-outputs/LarsenB_centerline.mat']).centerline.termdate;
termx = feval(fit(termdate',termx','poly2'),termdate');
% reformat surface dates to match terminus dates
hdates = zeros(1,length(h));
% set up figure
figure(1); clf
    plot(termdate,termx,'.-b','displayname','term','markersize',15);
    hold on; xlabel('date'); ylabel('distance along centerline (m)');
    set(gca,'fontsize',12,'linewidth',2);
    title('interpolated terminus positions'); legend;
% interpolate terminus position to surface date
for i=1:length(h)
    DV  = datevec(h(i).date);  % [N x 6] array
    DV  = DV(:, 1:3);   % [N x 3] array, no time
    DV2 = DV;
    DV2(:, 2:3) = 0;    % [N x 3], day before 01.Jan
    DV3 = cat(2, DV(:, 1), datenum(DV) - datenum(DV2));
    hdates(i) = DV3(1)+DV3(2)/365;
    % interpolate surface dates to terminus dates
    h(i).termx = interp1(termdate,termx,hdates(i));
    h(i).Iterm = dsearchn(x',h(i).termx);
    % plot result
    if i==1
        figure(1); plot(hdates(i),h(i).termx,'.-m','displayname','surf','MarkerSize',15);
    else
        plot(hdates(i),h(i).termx,'.-m','HandleVisibility','off','MarkerSize',15);
    end
end
clear DV*

%% 5. Width-averaged thickness and bed elevation profile
% Assume a parabolic bed and take mean thickness at each width segment

close all;

save_hb_adj = 1; % = 1 to save adjusted hb

cd([homepath,'inputs-outputs/']);

% load 2009 surface
h_2018 = load('surfaceElevationObs.mat').h(1).surface;
h_2018(find(isnan(h_2018),1,'first'):end) = 0;

% load width
W = load('calculatedWidth.mat').width.W;

% use bathymetry observations for centerline bed where they exist
%hb(~isnan(bathym)) = bathym(~isnan(bathym));

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

%% 6. Width-averaged speed
    
save_U_widthavg = 1; % = 1 to save results

% create line segments extending to glacier outline perpendicular to each
% centerline point

    % load glacier outline
    cd([homepath,'inputs-outputs/']);
    ol = load('glacierOutline.mat').ol;
    % use a staggered grid of centerline coordinates for calculating slope
    cl.Xm = (cl.X(2:end)+cl.X(1:end-1))./2;
    cl.Ym = (cl.Y(2:end)+cl.Y(1:end-1))./2;

    % initialize variables
    segs.l = 7000; % length of lines extending from centerline (m)
    dx = 10; % distance between points of line (m)
    segs.x = zeros(length(cl.X),2); % glacier width segments (m)
    segs.y = zeros(length(cl.X),2);
    m = zeros(1,length(cl.X)); m_rev = zeros(1,length(cl.X)); % slopes (unitless)
    segs.xn = zeros(length(cl.X),segs.l*2/dx+1); % segments with higher spatial resolution (m)
    segs.yn = zeros(length(cl.X),segs.l*2/dx+1);    

    figure(1); clf; hold on;
    plot(cl.X/10^3,cl.Y/10^3,'-m','linewidth',2,'displayname','centerline');
    plot(ol.x/10^3,ol.y/10^3,'-c','linewidth',2,'DisplayName','outline');
    set(gca,'fontsize',16,'linewidth',2); grid on; legend('Location','northwest');
    set(gcf,'position',[441 211 712 586]);
    xlabel('Easting (km)'); ylabel('Northing (km)');
    for i=2:length(cl.X)
        % Calculate slope angle at each point
        if i==length(cl.X)
            m(i) = m(i-1); 
            m_rev(i) = -1/m(i);
        else
            m(i) = (cl.Ym(i)-cl.Ym(i-1))/(cl.Xm(i)-cl.Xm(i-1));
            m_rev(i) = -1/m(i);
        end 

        % Calculate distance to points  
        delta_x = sqrt((segs.l)^2/(m_rev(i)^2+1));
        delta_y = sqrt((segs.l)^2/(m_rev(i)^2+1))*m_rev(i);

        % Add dy/dx to the point on the centerline
        segs.x(i,2) = cl.X(i)+delta_x; segs.x(i,1) = cl.X(i)-delta_x;
        segs.y(i,2) = cl.Y(i)+delta_y; segs.y(i,1) = cl.Y(i)-delta_y;        

        % Increase number of points in segment
        segs.xn(i,:) = segs.x(i,1):(segs.x(i,2)-segs.x(i,1))*(dx/(segs.l*2)):segs.x(i,2);
        segs.yn(i,:) = segs.y(i,1):(segs.y(i,2)-segs.y(i,1))*(dx/(segs.l*2)):segs.y(i,2);

        % Plot results on figure
        figure(1);
        if i==2
            plot(segs.xn(i,:)/10^3,segs.yn(i,:)/10^3,'-b','linewidth',1,'displayname','perpendicular lines');            
        else
            plot(segs.xn(i,:)/10^3,segs.yn(i,:)/10^3,'-b','linewidth',1,'handlevisibility','off');
        end

        % Once all other points are complete, go to the first point
        % and use the 2nd point slope to calculate segment 
        % (b/c the slope is 0 at i=1)
        if i==length(cl.X)
            m(1) = m(2); m_rev(1) = -1/m(1);
            % Calculate distance to points  
            delta_x = sqrt((segs.l)^2/(m_rev(1)^2+1));
            delta_y = sqrt((segs.l)^2/(m_rev(1)^2+1))*m_rev(1);        
            % Add dy/dx to the point on the centerline
            segs.x(1,2) = cl.X(1)+delta_x; segs.x(1,1) = cl.X(1)-delta_x;
            segs.y(1,2) = cl.Y(1)+delta_y; segs.y(1,1) = cl.Y(1)-delta_y; 
            % Increase number of points in segment
            segs.xn(1,:) = segs.x(1,1):(segs.x(1,2)-segs.x(1,1))*(dx/(segs.l*2)):segs.x(1,2);
            segs.yn(1,:) = segs.y(1,1):(segs.y(1,2)-segs.y(1,1))*(dx/(segs.l*2)):segs.y(1,2);            
            % Plot results
            plot(segs.xn(1,:)/10^3,segs.yn(1,:)/10^3,'-b','linewidth',1,'handlevisibility','off');
        end

    end

    % loop through each segment
    segs.xpts = segs.xn; segs.ypts = segs.yn;
    for i=1:length(segs.x(:,1))
        % If segment crossed through polygon
        if any(inpolygon(segs.xn(i,:),segs.yn(i,:),ol.x,ol.y))
            % Find the intersection points 
            [segs.intx(i,:),segs.inty(i,:)] = polyxpoly(ol.x,ol.y,segs.xn(i,:),segs.yn(i,:));
            % Calculate width from these segments
            segs.W(i) = sqrt((segs.intx(i,1)-segs.intx(i,2))^2+(segs.inty(i,1)-segs.inty(i,2))^2);
            % Make points outside polygon equal NaN
            segs.xpts(i,1:dsearchn([segs.xn(i,:)',segs.yn(i,:)'],[segs.intx(i,1),segs.inty(i,1)])) = NaN;
                segs.xpts(i,dsearchn([segs.xn(i,:)',segs.yn(i,:)'],[segs.intx(i,2),segs.inty(i,2)])+1:end) = NaN;
            segs.ypts(i,1:dsearchn([segs.xn(i,:)',segs.yn(i,:)'],[segs.intx(i,1),segs.inty(i,1)])) = NaN;
                segs.ypts(i,dsearchn([segs.xn(i,:)',segs.yn(i,:)'],[segs.intx(i,2),segs.inty(i,2)])+1:end) = NaN;
        else 
            segs.intx(i,:) = NaN;
            segs.inty(i,:) = NaN;
            segs.W(i) = NaN;
        end

        % Plot resulting segments onto figure above
        figure(1); hold on;
        if i==1
            plot(segs.intx(i,:)/10^3,segs.inty(i,:)/10^3,'LineWidth',2,'color',col(3,:),...
                'DisplayName','calculated width segments');
        else
            plot(segs.intx(i,:)/10^3,segs.inty(i,:)/10^3,'LineWidth',2,...
                'color',col(3,:),'HandleVisibility','off');        
        end
    end 

% load and interpolate velocities along width segments
cd([homepath,'data/velocities/']);
    % load ITS_LIVE speeds
    ILfiles = dir('ANT*.nc');
    % Loop through all files, interpolate along centerline
    for i=1:length(ILfiles)

        [X,Y] = meshgrid(ncread(ILfiles(i).name,'x'),ncread(ILfiles(i).name,'y'));
        u = ncread(ILfiles(i).name,'v')'; % m/y
        u(u==-32767) = NaN; %Replace no data values with NaN
        U(i).date = str2double(ILfiles(i).name(11:14));                    % observation date
        uw=NaN*zeros(length(cl.X),length(ii)); % initialize
        for j=1:length(cl.X)
            ii = find(~isnan(segs.xpts(j,:))); % find where segment coordinates are real
            for k=1:2:length(ii)
                uw(j,k) = interp2(X,Y,u,segs.xpts(j,ii(k)),segs.ypts(j,ii(k))); % interpolate speed
            end
            U(i).speed(j) = nanmean(uw)./3.1536e7;                          % width-averaged speed (m/s)
            disp([i,j]);
        end
        U(i).units = "m/s";                                                 % speed units
        U(i).source = "ITS-LIVE";                                           % speed data source
        U(i).numPts = length(U(i).speed(~isnan(U(i).speed)));               % number of data points
        % display mean error
        %disp([num2str(nanmean(U(i).date)),' mean error = ',num2str(nanmean(U(i).speed_err)),...
        %    ' m/s']);

        % Grab cropped 2018 .mat file on last iteration
        if i==length(ILfiles)
            load('ANT_G0240_2018_crop.mat');
            U(i+1).date = 2018;
            for j=1:length(cl.X)
                ii = find(~isnan(segs.xpts(j,:)));
                uw=zeros(1,length(ii));
                for k=1:length(ii)
                    uw(k) = interp2(X,Y,u,segs.xpts(j,ii(k)),segs.ypts(j,ii(k)));
                end
                U(i+1).speed(j) = nanmean(uw)./3.1536e7;                  % interpolated speed (m/s)
                disp([i,j]);
            end
            U(i+1).units = "m/s";                                                 % speed units
            U(i+1).source = "ITS-LIVE";                                           % speed data source
            U(i+1).numPts = length(U(i+1).speed(~isnan(U(i+1).speed)));           % number of data points
        end

    end

    % load TSX speeds
    TSXDates = [2007:2017];
    V = load('LarsenB_TSX_velocities.mat').V;
    for i=1:length(TSXDates)

        %Load velocity & coordinates
        x_v = {V(i).x}; x_v = x_v{:,:}; y_v = {V(i).y}; y_v = y_v{:,:}; % m/d
        v = {V(i).speed}; v = v{:,:}; 

        % Interpolate speed across width at each point
        int = i+length(ILfiles)+1;
        for j=1:length(cl.X)
            ii = find(~isnan(segs.xpts(j,:)));
            uw=zeros(1,length(ii));
            for k=1:length(ii)
                uw(k) = interp2(x_v,y_v,v,segs.xpts(j,ii(k)),segs.ypts(j,ii(k))); 
            end
            U(int).speed(j) = nanmean(uw)./(3600*24);   % interpolated speed (m/s)
            disp([i,j]);
        end

        %Save interpolated speed info as structure
        U(int).date = V(i).decidates(1);                % observation date
        U(int).units = 'm/s';                           % speed units
        U(int).source = 'TSX';                          % speed data source
        U(int).numPts = ...
            length(U(int).speed(~isnan(U(int).speed))); % number of data points

    end

    % sort and linearly extrapolate speeds to fill in gaps
    U = sortStruct(U,'date',1);
    % mean of first velocities where data exist
    U1 = nanmean([U(7).speed(1) U(11).speed(1) U(15).speed(1) ...
        U(18).speed(1) U(19).speed(1)],'all'); % m/s
    for k=1:length(U)
        % use the average speed at first point if NaN
        if isnan(U(k).speed(1))
            U(k).speed(1) = U1;
        end
        U(k).speed_linearextrap = feval(fit(x(~isnan(U(k).speed))',...
            U(k).speed(~isnan(U(k).speed)),'pchip'),x(1:135));
        U(k).speed_linearextrap(136:186) = U(k).speed_linearextrap(135);
    end

% plot results
figure(3); clf; hold on;
set(gca,'FontSize',18,'FontName','Arial','linewidth',2);
grid on; xlabel('Distance Along Centerline (km)'); ylabel('Speed (m a^{-1})');
legend('Location','northwest');
n = [5 8 13 15:20];
col = parula(length(n)+2);
for i=1:length(n)
    plot(x(1:135)/10^3,movmean(U(n(i)).speed(1:135),10)*3.1536e7,'linewidth',2,'color',col(i,:),'displayname',string(U(n(i)).date));
end

% save results
if save_U_widthavg
    cd([homepath,'inputs-outputs/']);
    U_widthavg = U;
    save('centerlineSpeedsWidthAveraged_2007-2018.mat','U_widthavg');
    disp('width-averaged speed saved');
end



