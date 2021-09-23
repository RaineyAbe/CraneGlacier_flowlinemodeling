% Rainey Aberle
% Summer/Fall 2020
% Script to create a timeseries of Crane Glacier geometry, ice speed, and 
% terminus position along glacier centerline, 
% used to tune glacier flowline model

% Order of operations:
%   1. Initial setup
%   2. Ice surface elevation (PGC & OIB), bed elevation profile (OIB), and surface speeds
%   (ITS_LIVE & TSX)
%   3. Terminus positions
%   4. OIB bed picks

%% 1. Initial setup

clear all; close all;

% Define homepath
homepath = '/Users/raineyaberle/Desktop/Research/CraneModeling/CraneGlacier_flowlinemodeling/';

% Add path to required functions
addpath([homepath,'matlabFunctions/hugheylab-nestedSortStruct']);
addpath([homepath,'matlabFunctions/']);

bed_save = 0;           % = 1 to save bed
surface_save = 0;       % = 1 to save surface
velocity_save = 0;      % = 1 to save velocity
terminus_save = 0;      % = 1 to save terminus positions
dHdt_save = 0;          % = 1 to save dHdt
figures_save = 0;       % = 1 to save figures 

% Load Crane centerline
cd([homepath,'inputs-outputs/']);
cl.X = load('Crane_centerline.mat').x; cl.Y = load('Crane_centerline.mat').y;

% Define x as distance along centerline
x = zeros(1,length(cl.X));
for i=2:(length(cl.X))
    x(i)=sqrt((cl.X(i)-cl.X(i-1))^2+(cl.Y(i)-cl.Y(i-1))^2)+x(i-1);
end

% Convert to lon/lat, load geoidheight at each pt along centerline
[cl.Lon,cl.Lat] = ps2wgs(cl.X,cl.Y,'StandardParallel',-71,'StandardMeridian',0);
h_geoid = geoidheight(cl.Lat,cl.Lon);

%% 2. Glacier ice surface elevation (PGC & OIB), bed (OIB), and surface speeds (ITS_LIVE & TSX)

%Note: OIB columns
%Elevation = 5
%Surface elevation = surface(5) - elevation(7)
%Thickness = 4
%Bed = elevation(5) - thickness(4)
%Convert elevations to geoid
%No data values = -9999

% If bed not previously saved, load and interpolate bed observations
if bed_save
    % Grab Tate's observed bed & surface, convert coordinates,
    % interpolate to centerline
    
    % Load Tate's OIB file
    obs = load('TateBed_IRMCR1B_20181016_01_005.mat');
    obs=obs.mdata_WGS84;
    hb0 = obs.bedElevation;
    h0 = obs.Surface_Elev;
    lat_hb = obs.Latitude;
    lon_hb = obs.Longitude;
    
    % Convert coordinates to Antarctic polar stereographic
    [x_hb0,y_hb0] = wgs2ps(lon_hb,lat_hb,'StandardParallel',-71,'StandardMeridian',0);
    
    % Interpolate bed elevation for each point along the centerline
    % within a certain distance
    maxDist = 5e3; % m
    hb = zeros(1,length(cl.X));
    for i=1:length(cl.X)
        Ihb = dsearchn([x_hb0 y_hb0],[cl.X(i) cl.Y(i)]);
        dist = sqrt((x_hb0(Ihb)-cl.X(i))^2+(y_hb0(Ihb)-cl.Y(i))^2);
        if dist<=maxDist
            hb(i) = hb0(Ihb);
        else
            hb(i)=NaN;
        end
    end
    
    hb(1) = hb(2);
    
    % save hb and coordinates in structure
    HB.hb0 = hb; HB.X = cl.X'; HB.Y = cl.Y';
    clear hb; hb=HB; clear HB;
    
    % Save resulting bed profile as variable
    cd([homepath,'inputs-outputs/']);
    save('observedBed_Tate.mat','hb');
    disp('bed variable saved');
    
else
    
    cd([homepath,'inputs-outputs/']);
    hb = load('observedBed_Tate.mat').hb.hb0;
    
end

% If surface not previously saved, load and interpolate surface
% observations
if surface_save
    % Grab OIB observation dates
    cd([homepath,'data/OIB/L3']);
    files = dir('Converted*');
    for i=1:length(files)
        OIBDates(i,:) = [num2str(files(i).name(18:21)),'-',num2str(files(i).name(22:23)),'-',num2str(files(i).name(24:25))]; % OIB temporal coverage
    end
    
    % Load PGC DEM surface data
    cd([homepath,'data/']);
    h = load('Crane_SurfaceObservations_PGCDEMs.mat').h;
    k = length(h); % Counter for number of fields in h
    
    % Set up figure, plot PGC observations
    col = parula(length(OIBDates(:,1))+k); % Color scheme for plot
    figure(1); hold on; set(gcf,'Units','centimeters','Position',[0 0 30 30]);
    grid on; xlabel('Distance Along Centerline (m)'); ylabel('Elevation (m)');
    xlim([0 7.5e4]); title('Crane Glacier Ice Surface Observations');
    set(gca,'FontName','Calibri','FontSize',13); legend;
    
    % Loop through all OIB files, interpolate observations along centerline,
    % save in h variable, plot all observations
    cd([homepath,'data/OIB/L3/']);
    for i=1:length(OIBDates(:,1))
        
        ib = readmatrix(files(i).name);
        ib_lat = ib(:,1); ib_lon = ib(:,2);
        
        %Create new matrix w polar stereographic (ps) coordinates
        P = wgs2ps(ib_lon,ib_lat,'StandardParallel',-71,'StandardMeridian',0);
        ps_ib = ib; ps_ib(:,1) = P(:,1); ps_ib(:,2) = P(:,2);
        ps_ib(ps_ib==-9999)=NaN; %Replace no data values with NaN
        
        % Find OIB points closest to each Crane centerline point
        maxDist = 1500;
        
        surf = [];
        
        for j=1:length(cl.X)
            n = dsearchn([ps_ib(:,1) ps_ib(:,2)],[cl.X(j) cl.Y(j)]);
            if pdist([ps_ib(n,1) ps_ib(n,2);cl.X(j) cl.Y(j)],'euclidean')<=maxDist
                surf(j) = ps_ib(n,5)-ps_ib(n,7)-h_geoid(j);
            else
                surf(j)= NaN;
            end
        end
        
        % Save in h structure
        k = k+1;
        h(k).surface = surf; %surf_m;
        h(k).date = OIBDates(i,:);
        h(k).x = x;
        h(k).source = 'OIB';
        
    end
    
    % Sort h by observation date, save variable
    [h,index] = nestedSortStruct(h,'date');
    cd([homepath,'inputs-outputs/']);
    save('surfaceElevationObs.mat','h');
    disp('surface variable saved');
    
else
    cd([homepath,'inputs-outputs/']);
    h = load('surfaceElevationObs.mat').h;
end

% Load bathymetry observations
bathym = load('bathymetryData.mat').cl_trough;
bathym = -1.*bathym; % Data reported in depth beneath the surface (convert to negative values for plotting)

% Load terminus positions
cd([homepath,'inputs-outputs/']);
termX = load('LarsenB_centerline.mat').centerline.termx;
termY = load('LarsenB_centerline.mat').centerline.termy;
termx = x(dsearchn([cl.X cl.Y],[termX' termY']));
termdate = load('LarsenB_centerline.mat').centerline.termdate;
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

% Plot surface
n = [1 2 3 5 10 13 15 18 29 36]; % indices of surfaces to use for annual tracking
figure(2); clf
subplot(2,1,1);
    ylabel('Elevation (m)'); xlabel('Distance Along Centerline (km)');
    set(gca,'fontname','Arial','fontsize',12,'linewidth',2,'fontweight','bold');
    xlim([0 70]); title('a) Observed Glacier Geometry');
    set(gcf,'Units','centimeters','position',[5 5 20 25]); grid on; legend;
    col = parula(length(h)+3); % Color scale for plotting surface lines
for i=1:length(h)
    if i==4 || i==16
        
    elseif any(n==i)
        hold on; plot(x(1:h(i).Iterm)./10^3,h(i).surface(1:h(i).Iterm),...
            '-','Color',col(i,:),'linewidth',1.5,'markersize',10,...
            'displayname',h(i).date(1:4));
        plot([x(h(i).Iterm)./10^3 x(h(i).Iterm)./10^3],...
            [hb(h(i).Iterm) h(i).surface(h(i).Iterm)],...
            'linewidth',1.5,'color',col(i,:),'handlevisibility','off');
    else
        hold on; plot(x(1:h(i).Iterm)./10^3,h(i).surface(1:h(i).Iterm),...
            '-','Color',col(i,:),'linewidth',1.5,'markersize',10,...
            'handlevisibility','off');
        plot([x(h(i).Iterm)./10^3 x(h(i).Iterm)./10^3],...
            [hb(h(i).Iterm) h(i).surface(h(i).Iterm)],...
            'linewidth',1.5,'color',col(i,:),'handlevisibility','off');
    end
end

% Plot observed bed profile and bathymetry observations
plot(x./10^3,hb,'-k','linewidth',2,'displayname','Bed (2018)');
plot(x(1:length(bathym))./10^3,bathym,'--c','linewidth',2,'displayname','Bathymetry (2006)');
hold off;

% calculate & plot dHdt
n=[1 2 3 9 11 21 22 25 28 30 35 36]; % use most continuous profiles (and at least one per year)
dHdt.dHdt = zeros(length(h),length(x)); dt = zeros(1,length(n));
for i=1:length(n)
    if i==1
        dHdt.dHdt(i,:) = dHdt.dHdt(i,:); % maintain zeros for year 0
    else
        dHdt.dt(i) = (datenum(h(n(i)).date)-datenum(h(n(i-1)).date)).*8.64e3; % (day*8.64e3s/day = s)
        dHdt.dHdt(i,:) = (h(n(i)).surface - h(n(i-1)).surface)./(dHdt.dt(i)); % m/s
    end
end
% dHdt total
dHdt.h1 = h(1).surface; dHdt.h2 = h(36).surface;
dHdt.dt_total = (datenum(h(36).date)-datenum(h(1).date)).*8.64e3; % s
dHdt.dH_total = dHdt.h1-dHdt.h2;
dHdt.dHdt_total = (dHdt.h1-dHdt.h2)./(dHdt.dt_total); % m/s
figure(3); clf; set(gcf,'Position',[386 285 1028 520]);
    hold on; grid on;
    yyaxis left; ylabel('ice surface (m)'); title('dH_{tot}');
    plot(x./10^3,h(1).surface,x./10^3,h(36).surface,'-b','linewidth',2);
    yyaxis right; plot(x./10^3,dHdt.dH_total,'--m','linewidth',2);
    set(gca,'fontsize',11,'linewidth',2); ylabel('dH (m)');

% save dHdt
if dHdt_save
    cd([homepath,'inputs-outputs/']);
    save('dHdt_2009-2018.mat','dHdt');
    disp('dHdt saved');
end

% If bed hasn't been saved, add the bed with bathymetry observations as
% final elevations
if bed_save
    cd([homepath,'inputs-outputs/']);
    save('observedBed.mat','hb_bathymobs','-append');
    disp('bed saved');
end

% Ice speeds (TSX, & ITS_LIVE)

% Set up subplot
figure(2);
subplot(2,1,2); hold on;
    set(gca,'FontSize',12,'FontName','Arial','linewidth',2,'fontweight','bold');
    grid on; xlabel('Distance Along Centerline (km)'); ylabel('Speed (m a^{-1})');
    title('b) Observed Glacier Speed'); legend('Location','northwest');

% Load and save velocities if not saved
if velocity_save
    
    cd([homepath,'data/velocities/']);
    
    % load ANT speeds
    ANTfiles = dir('ANT*.nc');
    % Loop through all files, interpolate along centerline
    for i=1:length(ANTfiles)
        
        [X,Y] = meshgrid(ncread(ANTfiles(i).name,'x'),ncread(ANTfiles(i).name,'y'));
        u = ncread(ANTfiles(i).name,'v')'; % m/y
        u_err = transpose(ncread(ANTfiles(i).name,'v_err'));
        u(u==-32767) = NaN; %Replace no data values with NaN
        
        U(i).date = str2double(ANTfiles(i).name(11:14));                % observation date
        U(i).speed = interp2(X,Y,u,cl.X,cl.Y)./3.1536e7;                % interpolated speed (m/s)
        U(i).speed_err = interp2(X,Y,u_err,cl.X,cl.Y)./3.1536e7;        % interpolated speed error (m/s)
        U(i).units = "m/s";                                             % speed units
        U(i).source = "ITS-LIVE";                                       % speed data source
        U(i).numPts = length(U(i).speed(~isnan(U(i).speed)));           % number of data points
        
        % display mean error
        disp([num2str(nanmean(U(i).date)),' mean error = ',num2str(nanmean(U(i).speed_err)),...
            ' m/s']);
    end
    
    % load TSX speeds
    TSXDates = [2007:2017];
    V = load('LarsenB_TSX_velocities.mat').V;
    for i=1:length(TSXDates)
        
        %Load velocity & coordinates
        x_v = {V(i).x}; x_v = x_v{:,:}; y_v = {V(i).y}; y_v = y_v{:,:};
        v = {V(i).speed}; v = v{:,:}; v_err = {V(i).speed_err}; v_err = v_err{:,:};
        
        %Interpolate speed and error along centerline
        u_cl = (interp2(x_v,y_v,v,cl.X,cl.Y))./(3600*24); %m s^-1
        u_cl_err = (interp2(x_v,y_v,v_err,cl.X,cl.Y))./(3600*24); %m s^-1
        
        %Save interpolated speed info as structure
        int = i+length(ANTfiles);
        U(int).date = V(i).decidates(1);                % observation date
        U(int).speed = u_cl;                            % interpolated speed (m/s)
        U(int).speed_err = u_cl_err;                    % interpolated speed error (m/s)
        U(int).units = 'm/s';                           % speed units
        U(int).source = 'TSX';                          % speed data source
        U(int).numPts = ...
            length(U(int).speed(~isnan(U(int).speed))); % number of data points
        
    end
    
    % sort and linearly extrapolate speeds to fill in gaps
    U = sortStruct(U,'date',1);
    % mean of first velocities where data exist
    U1 = nanmean([U(15).speed(1) U(16).speed(1) U(17).speed(1) ...
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
    
    % save u variable as structure
    cd([homepath,'inputs-outputs/']);
    save('centerlineSpeeds_2007-2018.mat','U');
    disp('velocity variable saved');
    
else
    cd([homepath,'inputs-outputs/']);
    U = load('centerlineSpeeds_2007-2018.mat').U;
end

% Plot velocities
n = [2 4 6 8 9 14 15:19];
col = parula(length(n)+2); % Color scheme for plotting u profiles
figure(2); subplot(2,1,2);
for i=1:length(n)
    plot(x./10^3,U(n(i)).speed.*3.1536e7,'-','color',col(i,:),'linewidth',2,...
        'displayname',num2str(round(U(n(i)).date)));
end

if figures_save
    cd([homepath,'figures/']);
    figure(2);
    saveas(gcf,'Crane_centerlineObs.png');
    disp('figure 2 saved');
end

hold off;

%% 3. Terminus Position

close all;

% Load LarsenB_centerline variable (Dryak & Enderlin)
cd([homepath,'inputs-outputs/']);
centerline = load('LarsenB_centerline.mat').centerline;

% Set up figures
figure(4); clf % Terminus positions along centerline
set(gcf,'Position',[50 100 1200 600]);
subplot(1,2,1);
    hold on; plot(cl.X,cl.Y,'-m','linewidth',2,'displayname','Centerline');
    title('Terminus Position Coordinates');
    set(gca,'FontSize',14,'FontName','Arial','linewidth',2);
    xlabel('Easting (m)'); ylabel('Northing (m)');
    xlim([-2.41e6 -2.405e6]); ylim([1.262e6 1.273e6]);
    legend('Location','northwest');
    grid on; %axis equal
subplot(1,2,2); % Terminus positions over time
    hold on;
    title('Terminus Position Over Time');
    set(gca,'FontSize',14,'FontName','Arial','LineWidth',2);
    xlabel('Year'); ylabel('Distance Along Centerline (km)');
    xlim([2000 2020]);
    grid on;

% Loop through all terminus positions/dates, plot
col = parula(length(centerline.termx)+10); % Color scheme for plotting term positions
for i=1:length(centerline.termx)
    
    termDate = centerline.termdate(i); % Date
    termx = centerline.termx(i); % Terminus x coord
    termy = centerline.termy(i); % Terminus y coord
    xn = dsearchn(cl.Y,termy);
    if mod(i-1,10) == 0 % Only include every 10 in legend
        subplot(1,2,1); hold on; plot(termx,termy,'.','color',col(i,:),'markersize',30,'displayname',num2str(termDate));
        subplot(1,2,2); hold on; plot(termDate,x(xn)/10^3,'.','color',col(i,:),'markersize',30);
    else
        subplot(1,2,1); hold on; plot(termx,termy,'.','color',col(i,:),'markersize',30,'handlevisibility','off');
        subplot(1,2,2); hold on; plot(termDate,x(xn)/10^3,'.','color',col(i,:),'markersize',30);
    end
    
    x_term(i,:) = ([termDate x(xn)]);
    
    term(i).decidate = termDate;
    term(i).X = termx;
    term(i).Y = termy;
    term(i).x = x(xn);
    
end

% Save term as .mat variable
if terminus_save
    cd([homepath,'inputs-outputs']);
    save('terminusPositions_2002-2019.mat','term');
end

% Save figure as image
if figures_save
    cd([homepath,'figures']);
    figure(4);
    saveas(gcf,'Crane_terminusPositions_2002-2019.png');
    disp('figure 4 saved');
end

%% 4. OIB bed picks

close all;

cd([homepath,'inputs-outputs/']);

% load OIB picks from previous years
OIB16_1 = load('OIBPicks_2016_005.mat').CraneOIBPicks_2016;
OIB16_2 = load('OIBPicks_2016_2.mat').CraneOIBPicks_2016;
OIB17 = load('OIBPicks_2017_006.mat').CraneOIBPicks_2017;
% interpolate coordinates to centerline
OIB16_1.x = x(dsearchn([cl.X cl.Y],[OIB16_1.Easting OIB16_1.Northing]))';
OIB16_2.x = x(dsearchn([cl.X cl.Y],[OIB16_2.Easting OIB16_2.Northing]))';
OIB17.x = x(dsearchn([cl.X cl.Y],[OIB17.Easting OIB17.Northing]))';

% plot with 2018 bed
figure(5); clf; hold on; grid on;
    set(gca,'fontname','Arial','fontsize',12,'linewidth',2,'fontweight','bold');
    ylabel('Elevation (m)'); xlabel('Distance Along Centerline (km)');
    xlim([0 70]); title('OIB Bed Picks');
    set(gcf,'Units','centimeters','position',[5 7.8 30 20]); grid on; legend;
    col = winter(4); % color scale for plotting bed elevation lines
    plot(x/10^3,hb,'color',col(1,:),'linewidth',2,'displayname','2018');
    plot(OIB17.x/10^3,OIB17.hb_geoid,'color',col(2,:),'linewidth',2,'displayname','2017');
    plot(OIB16_1.x/10^3,OIB16_1.hb_geoid,'color',col(3,:),'linewidth',2,'displayname','2016 (1)');
    plot(OIB16_2.x/10^3,OIB16_2.hb_geoid,'color',col(4,:),'linewidth',2,'displayname','2016 (2)');
    
%% 5. Width-averaged thickness and bed elevation profile
% Assume a parabolic bed and take mean thickness at each width segment

close all;

save_hb_adj = 1; % = 1 to save adjusted hb

cd([homepath,'inputs-outputs/']);

% load 2009 surface
h_2009 = load('surfaceElevationObs.mat').h(1).surface;
h_2009(find(isnan(h_2009),1,'first'):end) = 0;

% load width
W = load('calculatedWidth.mat').width.W;

% use bathymetry observations for centerline bed where they exist
hb(~isnan(bathym)) = bathym(~isnan(bathym));

% calculate thickness
H = h_2009-hb;

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
hb_adj = h_2009-H_adj;

% plot results
figure(6); clf; hold on; legend;
set(gcf,'position',[200 200 1000 800]); title('Width-Averaged Bed');
xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)');
set(gca,'linewidth',2,'fontsize',18); grid on;
plot(x/10^3,hb,'--k','linewidth',2,'DisplayName','hb_{cl}');
plot(x/10^3,hb_adj,'-k','linewidth',2,'DisplayName','hb_{adj}');
plot(x/10^3,h_2009,'-b','linewidth',2,'DisplayName','h');

% save results
if save_hb_adj
    save('delineatedBedWidthAveraged.mat','hb_adj','x');
    disp('hb_adj and hb_adj2 saved');
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
    ANTfiles = dir('ANT*.nc');
    % Loop through all files, interpolate along centerline
    for i=1:length(ANTfiles)

        [X,Y] = meshgrid(ncread(ANTfiles(i).name,'x'),ncread(ANTfiles(i).name,'y'));
        u = ncread(ANTfiles(i).name,'v')'; % m/y
        u(u==-32767) = NaN; %Replace no data values with NaN
        U(i).date = str2double(ANTfiles(i).name(11:14));                    % observation date
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
        if i==length(ANTfiles)
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
        int = i+length(ANTfiles)+1;
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



