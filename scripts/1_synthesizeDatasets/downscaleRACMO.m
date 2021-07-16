%% Script to load RACMO2.3 climate variables and statistically downscale SMB for elevation-dependence
% using methods adapted from Noel et al. (2016) 
% along the Crane Glacier centerline
% Rainey Aberle
% Spring 2020

% Walk through SMB downscaling method for one year 
%   (skip to #6 if computing for multiple years)
%   1. Load Crane centerline & mean 2011 RACMO variables 
%       (sf, sm, ro, smb, height), plot
%   2. Load RACMO height, interpolate along Crane centerline
%   3. Adjust sf, sm, and smb for elevation-dependence 
%       (Method adapted from Noel et al. (2016))
%   4. Use linear trendlines of sf and sm to calculate smb along centerline
%   5. Adjust air temperature T along Crane centerline using a 
%       dry adiabatic lapse rate
% Compute mean adjusted annual SMB and air temperature for 2009-2019
%   6. Reconduct steps 1-5 above for 2009-2019

close all; clear all;  

% Define homepath
homepath = '/Users/raineyaberle/Desktop/Research/CraneModeling/CraneGlacier_flowlinemodeling/';
cd(homepath);

% Add path with necessary functions, data, inputs/outputs
addpath([homepath,'matlabFunctions/']);
addpath([homepath,'matlabFunctions/cmocean_v2.0/cmocean/']);
addpath([homepath,'data/RACMO2.3/']); 
addpath([homepath,'inputs-outputs/']);

save_smb = 0; % = 1 to save resulting smb

%% 1. Load centerline, Annual Snowfall (sf), Snowmelt (sm), Runoff (ro),
%Surface Mass Balance (SMB), Air Temperature (airtemp) along centerline
    
% Load centerline
cl.X = load('Crane_centerline.mat','x').x; 
cl.Y = load('Crane_centerline.mat','y').y;
[cl.Lon,cl.Lat] = ps2wgs(cl.X,cl.Y,'StandardParallel',-71,'StandardMeridian',0);

% Densities of ice and water
rho_i = 917;    % kg/m^3
rho_w = 1000;   % kg/m^3

% Plot full gridded RACMO variables for year and Crane centerline
% Note: RACMO day # = (year - 1950)*365 + (day # in year)
% Ex: Feb 1, 2000 = (2000-1950)*365 + 1
yr = 2009;
day_start = (yr-1950)*365 + 1;
day_end = (yr-1950)*365 + 365;
        
% SNOWFALL (sf)
    sf.Lat = ncread('RACMO2.3p2_XPEN055_snowfall_monthly_1979_2016.nc','lat'); % degrees north
    sf.Lon = ncread('RACMO2.3p2_XPEN055_snowfall_monthly_1979_2016.nc','lon'); % degrees east
    sf.h = ncread('RACMO2.3p2_XPEN055_snowfall_monthly_1979_2016.nc','height'); % m above surface   
    sf.sf = ncread('RACMO2.3p2_XPEN055_snowfall_monthly_1979_2016.nc','snowfall'); % kg m-2 mo-1 
    time = ncread('RACMO2.3p2_ANT27_snowfall_monthly_1979_2016.nc','time'); % days since 1950/01/01
    Iday_start = dsearchn(time,day_start); % index of day start in RACMO time
    Iday_end = dsearchn(time,day_end); % index of day start in RACMO time
    
    % Grab mean sf for every month in year (RACMO days day_start:day_end)
    sf_yr = zeros(length(sf.sf(:,1,1)),length(sf.sf(1,:,1)),1); % initialize
    for i=1:length(sf.sf(:,1,1))
        for j=1:length(sf.sf(1,:,1))
            sf_yr(i,j,1) = nanmean(sf.sf(i,j,Iday_start:Iday_end))*12./rho_i; % m/yr
        end     
    end 
        
    % Extrapolate mean annual sf along Crane centerline (RACMO days day_start:day_end)
    % Grab RACMO grid
    RACMOy=zeros(1,length(cl.Lat)); RACMOx=zeros(1,length(cl.Lat)); % initialize
    for i=1:length(cl.Lat)
        lat_diff = abs(cl.Lat(i)*ones(size(sf.Lat)) - sf.Lat);
        lon_diff = abs(cl.Lon(i)*ones(size(sf.Lon)) - sf.Lon);
        diff_map = sqrt(lat_diff.^2+lon_diff.^2);
        RACMO_ref = find(diff_map==min(min(diff_map)));
        [RACMOy(i),RACMOx(i)] = ind2sub(size(squeeze(nanmean(sf.sf(:,:,1,270:450),4))),RACMO_ref);
    end 
    sf.cl=[];
    for i=1:length(RACMOy)
        for j=Iday_start:Iday_end
            sf.cl(i,j) = sf.sf(RACMOy(i),RACMOx(i),1,j);
        end    
    end  
    % Take the average at each point along centerline
    sf.cl(:,1:Iday_start-1)=[]; % Started adding points in column mo_start
    sf.cl(:,1) = nanmean(sf.cl,2);
    sf.cl(:,2:end) = []; % kg m-2 yr-1
    sf.cl = sf.cl./rho_i; % m yr-1
        
    figure(1); clf; hold on; 
        colormap(cmocean('algae'));
        set(gcf,'Units','centimeters','Position',[5 20 18 13]);
        surf(sf.Lon,sf.Lat,sf_yr); view(2);
        plot3(cl.Lon,cl.Lat,ones(length(cl.Lon))*3000,'-m','LineWidth',4,'DisplayName','Crane Centerline');
        set(gca,'FontName','Arial','FontSize',14,'linewidth',2);
        xlabel('Lon'); ylabel('Lat'); title(['RACMO Mean ',num2str(yr),' Snowfall']); 
        h = colorbar; set(get(h,'title'),'string','m yr^-^1');
        ylim([-66 -65]); xlim([-63.5 -61.5]); caxis([0 0.9.*max(sf_yr(:))]);
        hold off;
    
% SNOWMELT (sm)
    sm.Lat = ncread('RACMO2.3p2_XPEN055_snowmelt_monthly_1979_2016.nc','lat'); % degrees north
    sm.Lon = ncread('RACMO2.3p2_XPEN055_snowmelt_monthly_1979_2016.nc','lon'); % degrees east
    sm.h = ncread('RACMO2.3p2_XPEN055_snowmelt_monthly_1979_2016.nc','height'); % m above surface   
    sm.sm = squeeze(ncread('RACMO2.3p2_XPEN055_snowmelt_monthly_1979_2016.nc','snowmelt')); % kg m-2 mo-1
    
     %Grab mean sm for every month in year (RACMO days day_start:day_end)
    sm_yr = zeros(length(sm.sm(:,1,1)),length(sm.sm(1,:,1))); % initialize
    for i = 1:length(sm.sm(:,1,1))
        for j=1:length(sm.sm(1,:,1))
            sm_yr(i,j,1) = nanmean(sm.sm(i,j,Iday_start:Iday_end))*12; % kg/m2/yr
        end     
    end 
    sm_yr = sm_yr./rho_i; % m yr^-1

    % Interpolate mean annual sm along Crane centerline (RACMO days day_start:day_end)
    sm.cl=zeros(1,length(RACMOy));
    for i=1:length(RACMOy)
        for j=Iday_start:Iday_end
            sm.cl(i,j) = sm.sm(RACMOy(i),RACMOx(i),j);
        end    
    end  

   % Take the average at each point along centerline
    sm.cl(:,1:Iday_start-1)=[]; % Started adding points in column Iday_start
    sm.cl(:,1) = nanmean(sm.cl,2);
    sm.cl(:,2:end) = []; % kg m-2 yr-1
    sm.cl = sm.cl./rho_i; % m yr-1
        
    figure(2); clf; hold on; set(gcf,'Units','centimeters','Position',[5 0 18 13]);
        colormap(cmocean('algae'));
        surf(sm.Lon,sm.Lat,sm_yr); view(2);
        plot3(cl.Lon,cl.Lat,ones(length(cl.Lon))*3000,'-m','LineWidth',4,'DisplayName','Crane Centerline');
        set(gca,'FontName','Arial','FontSize',14,'linewidth',2);
        xlabel('Lon'); ylabel('Lat'); title(['RACMO Mean ',num2str(yr),' Snowmelt']);  
        h = colorbar; set(get(h,'title'),'string','m yr^-^1');
        ylim([-66 -65]); xlim([-63.5 -61.5]); caxis([0 0.9.*max(sm_yr(:))]);
        hold off;
    
% RUNOFF (ro)
    ro.Lat = ncread('RACMO2.3p2_XPEN055_runoff_monthly_1979_2016.nc','lat'); % degrees north
    ro.Lon = ncread('RACMO2.3p2_XPEN055_runoff_monthly_1979_2016.nc','lon'); % degrees east
    ro.h = ncread('RACMO2.3p2_XPEN055_runoff_monthly_1979_2016.nc','height'); % m above surface   
    ro.ro = squeeze(ncread('RACMO2.3p2_XPEN055_runoff_monthly_1979_2016.nc','runoff')); % kg m-2 mo-1
    
    % Grab mean annual ro for full RACMO grid (RACMO days day_start:day_end)
    ro_yr = zeros(length(ro.ro(:,1,1)),length(ro.ro(1,:,1))); % initialize
    for i = 1:length(ro.ro(:,1,1))
        for j=1:length(ro.ro(1,:,1))
            ro_yr(i,j,1) = nanmean(ro.ro(i,j,Iday_start:Iday_end))*12; % kg/m^2/yr
        end     
    end 
    ro_yr = ro_yr./rho_i; % m yr^-1

    % Interpolate mean annual ro along Crane centerline (RACMO days day_start:day_end)
    ro.cl=zeros(length(RACMOy),length(Iday_start:Iday_end)); % initialize
    for i=1:length(RACMOy)
        for j=Iday_start:Iday_end
            ro.cl(i,j) = ro.ro(RACMOy(i),RACMOx(i),j)*12; % kg/m^2/yr
        end    
    end  

    % Take the average at each point along centerline
    ro.cl(:,1:Iday_start-1)=[]; % Started adding points in column mo_start
    ro.cl(:,1) = nanmean(ro.cl,2);
    ro.cl(:,2:end) = []; % kg m-2 yr-1
    ro.cl = ro.cl./rho_i; % m yr-1

    figure(3); clf; hold on; set(gcf,'Units','centimeters','Position',[25 20 18 13]);
        colormap(cmocean('algae'));        
        surf(ro.Lon,ro.Lat,ro_yr); view(2);
        plot3(cl.Lon,cl.Lat,ones(length(cl.Lon))*3000,'-m','LineWidth',4,'DisplayName','Crane Centerline');
        set(gca,'FontName','Arial','FontSize',14,'linewidth',2);
        xlabel('Lon'); ylabel('Lat'); title(['RACMO Mean ',num2str(yr),' Runoff']);  
        h = colorbar; set(get(h,'title'),'string','m yr^-^1');
        ylim([-66 -65]); xlim([-63.5 -61.5]); caxis([0 0.9.*max(ro_yr(:))]);
        hold off;    

% SURFACE MASS BALANCE (SMB)
    smb.Lat = ncread('SMB_RACMO2.3p2_monthly_XPEN055_197901_201908.nc','lat'); %degrees north
    smb.Lon = ncread('SMB_RACMO2.3p2_monthly_XPEN055_197901_201908.nc','lon'); %degrees east
    smb.h = ncread('SMB_RACMO2.3p2_monthly_XPEN055_197901_201908.nc','height'); %m above the surface   
    smb.smb = squeeze(ncread('SMB_RACMO2.3p2_monthly_XPEN055_197901_201908.nc','smb')); %kg m-2 mo-1
    
    % Grab mean smb for every month in year (RACMO days Iday_start:Iday_end)
    smb_yr = zeros(length(smb.smb(:,1,1)),length(smb.smb(1,:,1))); % initialize
    for i = 1:length(smb.smb(:,1,1))
        for j=1:length(smb.smb(1,:,1))
            smb_yr(i,j,1) = nanmean(smb.smb(i,j,Iday_start:Iday_end))*12; % kg/m^2/yr
        end     
    end 
    smb_yr = smb_yr./rho_i; %m yr^-1
    smb_yr_cl = zeros(1,length(RACMOx));
    for i=1:length(RACMOx)
        smb.cl_yr(i) = smb_yr(RACMOy(i),RACMOx(i));
    end

    % Interpolate mean annual SMB along Crane centerline (RACMO days day_start:day_end)
    smb.cl=zeros(length(RACMOy),length(Iday_start:Iday_end)); % initialize
    for i=1:length(RACMOy)
        for j=Iday_start:Iday_end
            smb.cl(i,j) = smb.smb(RACMOy(i),RACMOx(i),j);
        end    
    end  

    % Take the average at each point along centerline
    smb.cl(:,1:Iday_start-1)=[]; % Started adding points in column Iday_start
    smb.cl(:,1) = nanmean(smb.cl,2);
    smb.cl(:,2:end) = []; % kg m-2 yr-1
    smb.cl = smb.cl./rho_i; % m yr-1

    figure(4); clf; hold on; set(gcf,'Units','centimeters','Position',[25 0 18 13]);
        colormap(cmocean('algae'));
        surf(smb.Lon,smb.Lat,smb_yr); view(2);
        plot3(cl.Lon,cl.Lat,ones(length(cl.Lon))*3000,'-m','LineWidth',4,'DisplayName','Crane Centerline');
        set(gca,'FontName','Arial','FontSize',14);
        xlabel('Lon'); ylabel('Lat'); title(['RACMO Mean ',num2str(yr),' SMB']);  
        h = colorbar; set(get(h,'title'),'string','m/yr');
        ylim([-66 -65]); xlim([-63.5 -61.5]); 
        caxis([0 0.9.*max(smb_yr(:))]);
        hold off;

%% 2. Load RACMO height, linearly interpolate along centerline
    
close all; clear h

h.Lat = ncread('Height_latlon_XPEN055.nc','lat'); %degrees north
h.Lon = ncread('Height_latlon_XPEN055.nc','lon'); %degrees east
h.h = ncread('Height_latlon_XPEN055.nc','height'); %m above the surface
h_cl_RACMO = griddata(h.Lon,h.Lat,h.h,cl.Lon,cl.Lat);

% Define x as distance along centerline
x=zeros(1,length(cl.X)); % initialize
for i=2:(length(cl.X))
    x(i)=sqrt((cl.X(i)-cl.X(i-1))^2+(cl.Y(i)-cl.Y(i-1))^2)+x(i-1);
end

%Plot RACMO gridded height variable and Crane centerline
figure(1); clf; hold on; 
    colormap(cmocean('algae'));
    set(gcf,'Units','centimeters','Position',[3 8 21 18]);
    surf(h.Lon,h.Lat,h.h); view(2);
    plot3(cl.Lon,cl.Lat,h_cl_RACMO,'-m','LineWidth',4,'DisplayName','Centerline');
    set(gca,'FontName','Arial','FontSize',14,'linewidth',2);
    xlabel('Lon'); ylabel('Lat'); title('RACMO Height'); 
    c = colorbar; set(get(c,'title'),'string','(m)');
    ylim([-66 -65]); xlim([-63.5 -61.5]);
    hold off;

figure(2); clf; hold on; 
    set(gcf,'Units','centimeters','Position',[28 8 21 18]);
    plot(x/10^3,h_cl_RACMO,'-k','LineWidth',3,'displayname','RACMO height'); grid on; 
    set(gca,'FontName','Arial','FontSize',14,'LineWidth',2); 
    title('RACMO height at Crane Centerline');
    xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)'); 
    hold off;
        
%% 3. Downscale SF, SM, & SMB  along centerline 
%   for elevation dependence, adapted from Noel et al.(2016)

close all;

%Convert RACMO grid to polar stereographic
[h.X,h.Y] = wgs2ps(h.Lon,h.Lat,'StandardParallel',-71,'StandardMeridian',0);
[sf.X,sf.Y] = wgs2ps(sf.Lon,sf.Lat,'StandardParallel',-71,'StandardMeridian',0);
[sm.X,sm.Y] = wgs2ps(sm.Lon,sm.Lat,'StandardParallel',-71,'StandardMeridian',0);
[smb.X,smb.Y] = wgs2ps(smb.Lon,smb.Lat,'StandardParallel',-71,'StandardMeridian',0);

%Load Crane ice surface
h_cl = load("surfaceElevationObs.mat").h(28).surface';

maxDist = 10e3; 
   
%Adjust snowfall (sf) 
    
    %Loop through each point along centerline
    sf.interp = ones(length(h_cl),1);
    for i=1:length(x)
        
        %Grab points within a certain distance from centerline, interpolate
        nearbyPts = zeros(length(h.X(:,1)),length(h.X(1,:))); %Hold points within a certain distance
        numPts = 0; %Total number of points found
    
        xi = cl.X(i); yi = cl.Y(i);
        for j=1:length(h.X(:,1))
            for k=1:length(h.X(1,:))
                dist = sqrt((xi-h.X(j,k))^2 + (yi-h.Y(j,k))^2);
                if dist<=maxDist
                    numPts = numPts+1;
                    nearbyPts(numPts,1:2) = ([h.h(j,k) sf_yr(j,k)]); 
                end 
            end 
        end 
        
        %Calculate a linear trendline for nearby points
        P = polyfit(nearbyPts(:,1),nearbyPts(:,2),1);
        int = 0:1600; 
        yfit = P(1)*int+P(2);
    
        %Apply regression slope to current grid cell, solve for eqn
        b = sf.cl(i)-P(1)*h_cl_RACMO(i);
        sf.interp(i) = P(1)*h_cl(i)+b; 
        
    end 
    
%Ajust snowmelt (sm)
    
    %Loop through each point along centerline
    sm.interp = ones(length(h_cl),1);
    for i=1:length(x)
    
        %Grab points within a certain distance from centerline, interpolate
        nearbyPts = []; %Hold points within a certain distance
        numPts = 0; %Total number of points found
    
        xi = cl.X(i); yi = cl.Y(i);
        for j=1:length(h.X(:,1))
            for k=1:length(h.X(1,:))
                dist = sqrt((xi-h.X(j,k))^2 + (yi-h.Y(j,k))^2);
                if dist<=maxDist
                    numPts = numPts+1;
                    nearbyPts(numPts,1:2) = ([h.h(j,k) sm_yr(j,k)]); 
                end 
            end 
        end 
        
        %Calculate a linear trendline for nearby points
        P = polyfit(nearbyPts(:,1),nearbyPts(:,2),1);
        int = 0:1600; 
        yfit = P(1)*int+P(2);
    
        %Apply regression slope to current grid cell, solve for eqn
        m = P(1); b = sm.cl(i)-m*h_cl_RACMO(i);
        smfit = [h_cl_RACMO m.*h_cl_RACMO+b];
        sm.interp(i) = m*h_cl(i)+b; 
        
    end 

%Ajust surface mass balance (smb)

    %Loop through each point along centerline
    smb.interp = ones(length(h_cl),1);
    for i=1:length(x)
    
        %Grab points within a certain distance from centerline, interpolate
        nearbyPts = []; %Hold points within a certain distance
        numPts = 0; %Total number of points found
    
        xi = cl.X(i); yi = cl.Y(i);
        for j=1:length(h.X(:,1))
            for k=1:length(h.X(1,:))
                dist = sqrt((xi-h.X(j,k))^2 + (yi-h.Y(j,k))^2);
                if dist<=maxDist
                    numPts = numPts+1;
                    nearbyPts(numPts,1:2) = ([h.h(j,k) smb_yr(j,k)]); 
                end 
            end 
        end 
        
        %Calculate a linear trendline for nearby points
        P = polyfit(nearbyPts(:,1),nearbyPts(:,2),1);
        int = -10:1800; 
        yfit = P(1)*int+P(2);
    
        %Apply regression slope to current grid cell, solve for eqn
        m = P(1); b = smb.cl(i)-m*h_cl_RACMO(i);
        smbfit = [h_cl_RACMO m.*h_cl_RACMO+b];
        smb.interp(i) = m*h_cl(i)+b; 
        
        if i==1
            figure(1); clf; hold on; grid on; set(gcf,'Units','centimeters','Position',[3 8 21 18]);
            set(gca,'FontSize',14,'FontName','Arial'); 
            title('Regression Calculation for SMB, i=1');
            plot(smbfit(:,1),smbfit(:,2),'-b','LineWidth',2,'DisplayName','a');            
            plot(yfit,'-r','LineWidth',2,'DisplayName','b'); 
            plot(h_cl_RACMO(i),smb.cl(i),'.b','MarkerSize',30,'DisplayName','Current grid cell');
            plot(nearbyPts(:,1),nearbyPts(:,2),'.r','MarkerSize',30,'DisplayName','Adjacent grid cells');
            plot(h_cl(i),smb.interp(i),'*b','MarkerSize',25,'DisplayName','Resulting SMB');
            xlabel('Height (m)'); ylabel('SMB (m yr^-^1)');
            legend('Location','southeast');
            hold off; 
        end 
    end 

%Plot adjusted climate variables
figure(2); clf; hold on; set(gcf,'Units','centimeters','Position',[28 8 21 18]);
    plot(h_cl,sf.interp,'.b','MarkerSize',15,'DisplayName','Snowfall'); 
    plot(h_cl,sm.interp,'.m','MarkerSize',15,'DisplayName','Snowmelt');
    plot(h_cl,smb.interp,'.g','MarkerSize',15,'DisplayName','SMB');
    grid on; title('Adjusted Climate Variables');
    set(gca,'FontSize',14,'FontName','Arial');
    xlabel('Height (m)'); legend; ylabel('(m yr^-^1)');
    hold off;
    
%% 4. Use linear trendline to calculate smb slope

close all;

%Linear trend in adjusted climate variables
sf.linear = fit(h_cl(~isnan(h_cl)),sf.interp(~isnan(h_cl)),'poly1'); 
    sf.linear = feval(sf.linear,h_cl);
sm.linear = fit(h_cl(~isnan(h_cl)),sm.interp(~isnan(h_cl)),'poly1');
    sm.linear_eval = sm.linear.p1*h_cl+sm.linear.p2; 
    sm.linear = feval(sm.linear,h_cl);
smb.linear = fit(h_cl(~isnan(h_cl)),smb.interp(~isnan(h_cl)),'poly1');
    smb.linear = feval(smb.linear,h_cl);
smb.linear2 = sf.linear-sm.linear;

%Plot Results
figure(1); clf; hold on; 
    grid on; set(gcf,'Units','centimeters','Position',[3 8 21 18]);
    set(gca,'FontName','Arial','FontSize',14,'linewidth',2); 
    plot(h_cl,sf.linear,'.c','MarkerSize',12,'DisplayName','Snowfall');
    plot(h_cl,sm.linear,'.m','MarkerSize',12,'DisplayName','Snowmelt');
    plot(h_cl,smb.linear,'.b','MarkerSize',12,'DisplayName','SMB');
    plot(h_cl,smb.linear2,'.-r','MarkerSize',10,'DisplayName','SMB (SF-SM)');
    xlabel('Height (m)'); ylabel('m yr^-^1'); title('Linear Trends in Adjusted Climate Variables');
    legend('Location','east');
    hold off;

figure(2); clf; hold on; 
    grid on; set(gcf,'Units','centimeters','Position',[28 8 21 18]);
    set(gca,'FontName','Arial','FontSize',14,'linewidth',2); 
    plot(x,smb.interp,'x-b','MarkerSize',14,'DisplayName','Adjusted');
    plot(x,smb.linear,'--b','LineWidth',2,'DisplayName','Adjusted Linear Trend');
    plot(x,smb.cl,'-b','LineWidth',2,'DisplayName','UN-Adjusted');
    plot(x,smb.linear2,'.-r','MarkerSize',10,'DisplayName','SMB (SF-SM)');    
    title('SMB Along Crane Centerline'); legend('Location','southwest');
    xlabel('Distance Along Centerline (m)'); ylabel('SMB (m yr^-^1)');
    hold off;  

%% 5. Adjust air temperature at ice surface using a dry adiabatic lapse rate

close all;

% Load RACMO air temperature
cd([homepath,'data/RACMO2.3']);
T.Lat = ncread('RACMO2.3p2_XPEN055_T2m_monthly_1979_2016.nc','lat'); % degrees north
T.Lon = ncread('RACMO2.3p2_XPEN055_T2m_monthly_1979_2016.nc','lon'); % degrees east
T.h = ncread('RACMO2.3p2_XPEN055_T2m_monthly_1979_2016.nc','height'); % m above height variable   
T.T = squeeze(ncread('RACMO2.3p2_XPEN055_T2m_monthly_1979_2016.nc','t2m')); % degrees K

T.T_yr = zeros(length(T.T(:,1,1)),length(T.T(1,:,1))); % initialize
for i = 1:length(T.T(:,1,1))
    for j=1:length(T.T(1,:,1))
        T.T_yr(i,j,1) = nanmean(T.T(i,j,Iday_start:Iday_end));
    end     
end 
T.T_yr = T.T_yr - 273.15; %degrees C

%Interpolate mean annual T along centerline (RACMO days day_start:day_end)
T.cl=zeros(length(RACMOy),length(Iday_start:Iday_end));
for i=1:length(RACMOy)
    for j=Iday_start:Iday_end
        T.cl(i,j) = T.T(RACMOy(i),RACMOx(i),j);
    end    
end  

% Take the average at each point along centerline
T.cl(:,1:Iday_start-1)=[]; % started adding points in column Iday_start
T.cl(:,1) = nanmean(T.cl,2);
T.cl(:,2:end) = []; % deg K
T.cl = T.cl - 273.15; % degrees C

% Loop through each point along centerline
T.interp = ones(length(h_cl),1);
for i=1:length(x)

    % Grab points within a certain distance from centerline, interpolate
    maxDist = 9e3; 
    nearbyPts = []; % Hold points within a certain distance
    numPts = 0; % Total number of points found

    xi = cl.X(i); yi = cl.Y(i);
    for j=1:length(h.X(:,1))
        for k=1:length(h.X(1,:))
            dist = sqrt((xi-h.X(j,k))^2 + (yi-h.Y(j,k))^2);
            if dist<=maxDist
                numPts = numPts+1;
                nearbyPts(numPts,1:2) = ([h.h(j,k) T.T_yr(j,k)]); 
            end 
        end 
    end 

    %Calculate a linear trendline for nearby points
    P = polyfit(nearbyPts(:,1),nearbyPts(:,2),1);
    int = -10:1800; 
    yfit = P(1)*int+P(2);

    %Apply regression slope to current grid cell, solve for eqn
    m = P(1); b = T.cl(i)-m*h_cl_RACMO(i);
    T.fit = [h_cl_RACMO m.*h_cl_RACMO+b];
    T.interp(i) = m*h_cl(i)+b; 
end 

%Dry adiabatic lapse rate    
lr = 9.8e-3; %lapse rate (degrees C m^-1)

%Calculate new air temperature using lapse rate and RACMO reference height
T.cl_m = movmean(T.cl,10,'omitnan');
h_cl_m = movmean(h_cl,10,'omitnan');
T.cl_yr = (h_cl_RACMO-h_cl)*lr+T.cl_m; 

%Plot height profiles
figure(1); clf; hold on; 
    set(gcf,'Units','centimeters','Position',[2 10 22 20]);
    plot(x/10^3,h_cl_RACMO,'-b','LineWidth',3,'DisplayName','RACMO height'); 
    plot(x/10^3,h_cl,'-m','LineWidth',3,'DisplayName','Crane Surface');
    grid on; xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)');
    set(gca,'FontSize',14,'FontName','Arial'); legend('Location','south'); 
    title('Height Profiles');
    hold off;

%Plot temperature profile at each surface
figure(2); clf; hold on; 
    set(gcf,'Units','centimeters','Position',[26 10 22 20]);
    plot(x/10^3,T.interp,'-b','LineWidth',3,'DisplayName','RACMO height'); 
    plot(x/10^3,T.cl_yr,'-m','LineWidth',3,'DisplayName','Crane Surface');
    grid on; xlabel('Distance Along Centerline (km)'); ylabel('Temperature (^oC)');
    set(gca,'FontSize',14,'FontName','Arial'); legend('Location','south');
    title('RACMO Air Temperature');
    hold off;

if save_smb
    cd([homepath,'data/RACMO2.3/']);
    save_smb('downscaledSMB_2011.mat','x','h_cl_11','smb.interp','smb.linear');
    save_smb('adjustedAirTemp_2011.mat','T.cl_11');
    disp('downscaled SMB and air temperature saved for one year.');
end 

%% 6. Repeat above steps to load smb for 2009-2019

close all; 

figure_save = 0;    % = 1 to save figure
save_variables = 0;       % = 1 to save final downscaled SMB

years = 2009:2019; % Define years
maxDist = 8e3;  % Define distance from centerline over which values will 
                % be used in the downscaling computation 

col = parula(length(years)+1); % Define color scheme for plotting

% Load and initialize variables
loop=1;
while loop==1

    % Load centerline
    cl.X = load('Crane_centerline.mat','x').x; 
    cl.Y = load('Crane_centerline.mat','y').y;
    cd([homepath,'inputs-outputs/']);
        
    % Load most advanced terminus position (2019)
    term = dsearchn([cl.X cl.Y],[load('terminusPositions_2002-2019.mat').term(61).X ...
    load('terminusPositions_2002-2019.mat').term(61).Y]); 
    % clip centerline at this point
    %cl.X = cl.X(1:term); cl.Y = cl.Y(1:term);
    
    % Define x as distance along centerline
    x = zeros(1,length(cl.X));
    for i=2:length(cl.X)
        x(i)=sqrt((cl.X(i)-cl.X(i-1))^2+(cl.Y(i)-cl.Y(i-1))^2)+x(i-1);
    end 
    
    % Convert coords to lon lat
    [cl.Lon,cl.Lat] = ps2wgs(cl.X,cl.Y,'StandardParallel',-71,'StandardMeridian',0);
    
    rho_i = 917; % kg/m^3
    
    % Load RACMO variables
    SF.Lat = ncread('RACMO2.3p2_XPEN055_snowfall_monthly_1979_2016.nc','lat'); % degrees north
    SF.Lon = ncread('RACMO2.3p2_XPEN055_snowfall_monthly_1979_2016.nc','lon'); % degrees east
    SF.sf = squeeze(ncread('RACMO2.3p2_XPEN055_snowfall_monthly_1979_2016.nc','snowfall')); % kg m-2 mo-1 
    time = ncread('RACMO2.3p2_XPEN055_snowfall_monthly_1979_2016.nc','time'); % days since 1950/01/01
    
    SM.Lat = ncread('RACMO2.3p2_XPEN055_snowmelt_monthly_1979_2016.nc','lat'); % degrees north
    SM.Lon = ncread('RACMO2.3p2_XPEN055_snowmelt_monthly_1979_2016.nc','lon'); % degrees east
    SM.sm = squeeze(ncread('RACMO2.3p2_XPEN055_snowmelt_monthly_1979_2016.nc','snowmelt')); % kg m-2 mo-1
    
    RO.Lat = ncread('RACMO2.3p2_XPEN055_runoff_monthly_1979_2016.nc','lat'); % degrees north
    RO.Lon = ncread('RACMO2.3p2_XPEN055_runoff_monthly_1979_2016.nc','lon'); % degrees east
    RO.ro = squeeze(ncread('RACMO2.3p2_XPEN055_runoff_monthly_1979_2016.nc','runoff')); % kg m-2 mo-1
    
    SMB.Lat = ncread('SMB_RACMO2.3p2_monthly_XPEN055_197901_201908.nc','lat'); % degrees north
    SMB.Lon = ncread('SMB_RACMO2.3p2_monthly_XPEN055_197901_201908.nc','lon'); % degrees east
    SMB.smb = squeeze(ncread('SMB_RACMO2.3p2_monthly_XPEN055_197901_201908.nc','smb')); % kg m-2 mo-1
    
    h.Lat = ncread('Height_latlon_XPEN055.nc','lat'); % degrees north
    h.Lon = ncread('Height_latlon_XPEN055.nc','lon'); % degrees east
    h.h = ncread('Height_latlon_XPEN055.nc','height'); % m above the surface
    h.cl = griddata(h.Lon,h.Lat,h.h,cl.Lon,cl.Lat); 
    
    T.Lat = ncread('RACMO2.3p2_XPEN055_T2m_monthly_1979_2016.nc','lat'); % degrees north
    T.Lon = ncread('RACMO2.3p2_XPEN055_T2m_monthly_1979_2016.nc','lon'); % degrees east
    T.T = squeeze(ncread('RACMO2.3p2_XPEN055_T2m_monthly_1979_2016.nc','t2m')); % degrees K
    
    % Grab RACMO grid
        % Note: RACMO day # = (year - 1950)*365 + (day # in year)
        % Ex: 1 Feb/2000 = (2000-1950)*365 + 1 = 18251
        for i=1:length(cl.Lon)
            lat_diff = abs(cl.Lat(i)*ones(size(SF.Lat)) - SF.Lat);
            lon_diff = abs(cl.Lat(i,1)*ones(size(SF.Lon)) - SF.Lon);
            diff_map = sqrt(lat_diff.^2+lon_diff.^2);
            RACMO_ref = find(diff_map==min(min(diff_map)));
            [RACMOy(i),RACMOx(i)] = ind2sub(size(nanmean(SF.sf(:,:,1),4)),RACMO_ref);
        end     
        
    % Convert RACMO grid to polar stereographic
    [h.X,h.Y] = wgs2ps(h.Lon,h.Lat,'StandardParallel',-71,'StandardMeridian',0);
    [SF.X,SF.Y] = wgs2ps(SF.Lon,SF.Lat,'StandardParallel',-71,'StandardMeridian',0);
    [SM.X,SM.Y] = wgs2ps(SM.Lon,SM.Lat,'StandardParallel',-71,'StandardMeridian',0);
    [SMB.X,SMB.Y] = wgs2ps(SMB.Lon,SMB.Lat,'StandardParallel',-71,'StandardMeridian',0);
    
    % Load 2016 Crane ice surface
    cd([homepath,'inputs-outputs/']);
    h_cl = load("surfaceElevationObs.mat").h(28).surface';   
        
    % Grab points within a certain distance from centerline, interpolate
    nearbyPts = []; % Hold points within a certain distance
    numPts = 0; % Total number of points found    
    
    % Initialize variables
    SF.extrap = NaN*zeros(length(years),length(x)); SF.downscaled = NaN*zeros(length(years),length(x)); SF.linear = NaN*zeros(length(years),length(x)); SF.fullGrid = NaN*zeros(length(RACMOx),length(RACMOy));
    SM.extrap = NaN*zeros(length(years),length(x)); SM.downscaled = NaN*zeros(length(years),length(x)); SM.linear = NaN*zeros(length(years),length(x)); SM.fullGrid = NaN*zeros(length(RACMOx),length(RACMOy));
    SMB.extrap = NaN*zeros(length(years),length(x)); SMB.downscaled = NaN*zeros(length(years),length(x)); SMB.linear = NaN*zeros(length(years),length(x)); SMB.fullGrid = NaN*zeros(length(RACMOx),length(RACMOy));

    loop=loop+1;

end

% Loop through years to calculate Adjusted Snowfall (SF), Snowmelt (SM), 
% and Surface Mass Balance (SMB) along centerline
for i=1:length(years)

    day_start = (years(i)-1950)*365+1; % first day in year in RACMO days
        Iday_start = dsearchn(time,day_start); % index of first day in RACMO days
    day_end = (years(i)-1950)*365+365; % last day in year in RACMO days
        Iday_end = dsearchn(time,day_end); % index of last day in RACMO days
    
    % set up figure
    if i==1
        figure(1); clf
        set(gcf,'units','centimeters','position',[5 0 40 30]);
        subplot(2,3,1); hold on;  % adjusted snowfall
            set(gca,'FontName','Arial','FontSize',14);
            xlabel('Distance Along Centerline (km)'); ylabel('(m a^{-1})'); 
            title('Downscaled Snowfall');grid on;
        subplot(2,3,2); hold on; % adjusted snowmelt
            set(gca,'FontName','Arial','FontSize',14);
            xlabel('Distance Along Centerline (km)'); ylabel('(m a^{-1})'); 
            title('Downscaled Snowmelt');grid on; 
        subplot(2,3,3); hold on; % adjusted SMB
            set(gca,'FontName','Arial','FontSize',14);
            xlabel('Distance Along Centerline (km)'); ylabel('(m a^{-1})'); 
            title('Downscaled SMB');grid on;         
        subplot(2,3,[4,4.5]); hold on; % linearized SMB
            set(gca,'FontName','Arial','FontSize',14);
            xlabel('Distance Along Centerline (km)'); ylabel('(m a^{-1})'); 
            title('Estimated Runoff');grid on; 
        subplot(2,3,[5.5,6]); hold on; % adjusted air temp
            set(gca,'FontName','Arial','FontSize',14);
            xlabel('Distance Along Centerline (km)'); ylabel('(^oC)'); 
            title('Adjusted Air Temp'); grid on; 
        c = colorbar('position',[0.93 0.4 0.02 0.3],'fontsize',14);
        c.Ticks = [0 0.5 1.0]; c.TickLabels = {'2002','2010','2019'};
    end 

    % SNOWFALL (SF)
    if years(i)<=2016

        % Grab mean SF for every month for full RACMO grid
        for j=1:length(SF.sf(:,1,1))
            for k=1:length(SF.sf(1,:,1))
                SF.fullGrid(j,k) = nanmean(SF.sf(j,k,Iday_start:Iday_end))./rho_i*12; % m a^-1
            end  
        end 
        % Extrapolate mean annual SF along Crane centerline
        for j=1:length(x)
            SF.extrap(i,j) = SF.fullGrid(RACMOy(j),RACMOx(j)); 
        end

        % Adjust snowfall (SF) for elevation-dependence

        % Loop through each point along centerline
        for j=1:length(x)

            % Grab points within a certain distance from centerline, interpolate
            nearbyPts = []; % Hold points within a certain distance
            numPts = 0; % Total number of points found

            for k=1:length(h.X(:,1))
                for l=1:length(h.X(1,:))
                    dist = sqrt((cl.X(j)-h.X(k,l))^2 + (cl.Y(j)-h.Y(k,l))^2);
                    if dist<=maxDist
                        numPts = numPts+1;
                        nearbyPts(numPts,1:2) = ([h.h(k,l) SF.fullGrid(k,l)]); 
                    end 
                end 
            end 

            % Grab the mean of nearby points
            SF.mean(i,j) = nanmean(nearbyPts(:,2));
            
            %Calculate a linear trendline for nearby points
            P = polyfit(nearbyPts(:,1),nearbyPts(:,2),1);
            int = -500:2000; 
            yfit = P(1)*int+P(2);
    
            SF.downscaled(i,j) = interp1(int,yfit,h_cl(j));
    
        end     
    
        % Fit a linear trendline to snowfall
        SF.linear(i,~isnan(h_cl)) = feval(fit(h_cl(~isnan(h_cl)),SF.downscaled(i,~isnan(h_cl))','poly1'),h_cl(~isnan(h_cl)))';
        SF.linear(i,find(isnan(h_cl),1,'first'):length(x)) = SF.linear(i,find(isnan(h_cl),1,'first')-1);
    
        % Plot
        figure(1); subplot(2,3,1);
        plot(x/10^3,SF.downscaled(i,:),'--','color',col(i,:),'linewidth',1,'displayname',num2str(years(i))); 
        plot(x/10^3,SF.linear(i,:),'-','color',col(i,:),'linewidth',2,'HandleVisibility','off');     
        drawnow
        
        % Calculate and plot the mean SF at each point 
        if years(i)==2016
            SF.downscaled_average = nanmean(SF.downscaled,1);
            SF.downscaled_average_linear = feval(fit(h_cl(~isnan(h_cl)),SF.downscaled_average(~isnan(h_cl))','poly1'),h_cl);
            plot(x/10^3,SF.downscaled_average_linear,'-k','linewidth',2);
        end

    end   

    % SNOWMELT (SM)
    if years(i)<=2016

        % Grab mean SM for every month for full RACMO grid
        for j=1:length(SM.sm(:,1,1))
            for k=1:length(SM.sm(1,:,1))
                SM.fullGrid(j,k) = nanmean(SM.sm(j,k,Iday_start:Iday_end))./rho_i.*12; % m a-1
            end   
        end 
        % Extrapolate mean annual SM along Crane centerline
        for j=1:length(x)
            SM.extrap(i,j) = SM.fullGrid(RACMOy(j),RACMOx(j));
        end

        % Adjust SM for elevation-dependence

        % Loop through each point along centerline
        for j=1:length(x)

            % Grab points within a certain distance from centerline, interpolate
            nearbyPts = []; % Hold points within a certain distance
            numPts = 0; % Total number of points found

            for k=1:length(h.X(:,1))
                for l=1:length(h.X(1,:))
                    dist = sqrt((cl.X(j)-h.X(k,l))^2 + (cl.Y(j)-h.Y(k,l))^2);
                    if dist<=maxDist
                        numPts = numPts+1;
                        nearbyPts(numPts,1:2) = ([h.h(k,l) SM.fullGrid(k,l)]); 
                    end 
                end 
            end 

            % Grab the mean of nearby points
            SM.mean(i,j) = nanmean(nearbyPts(:,2));
            
            %Calculate a linear trendline for nearby points
            P = polyfit(nearbyPts(:,1),nearbyPts(:,2),1);
            int = -500:2000; 
            yfit = P(1)*int+P(2);
    
            SM.downscaled(i,j) = interp1(int,yfit,h_cl(j));
    
        end     
    
        % Fit a linear trendline to the downscaled SM
        SM.linear(i,~isnan(h_cl)) = feval(fit(h_cl(~isnan(h_cl)),SM.downscaled(i,~isnan(h_cl))','poly1'),h_cl(~isnan(h_cl)))';
        SM.linear(i,find(isnan(h_cl),1,'first'):length(x)) = SM.linear(i,find(isnan(h_cl),1,'first')-1);
    
        % Plot
        figure(1); subplot(2,3,2);
        plot(x/10^3,SM.downscaled(i,:),'--','color',col(i,:),'linewidth',1,'displayname',num2str(years(i))); 
        plot(x/10^3,SM.linear(i,:),'-','color',col(i,:),'linewidth',2,'HandleVisibility','off');     
        drawnow

        % Calculate and plot the mean SM at each point 
        if years(i)==2016
            SM.downscaled_average = nanmean(SM.downscaled,1);
            SM.downscaled_average_linear = feval(fit(h_cl(~isnan(h_cl)),SM.downscaled_average(~isnan(h_cl))','poly1'),h_cl);
            plot(x/10^3,SM.downscaled_average_linear,'-k','linewidth',2);
        end
        
    end   

    % SURFACE MASS BALANCE (SMB)
    if years(i)<=2019

        % Grab mean SMB for every month for full RACMO grid
        for j = 1:length(SMB.smb(:,1,1))
            for k=1:length(SMB.smb(1,:,1))
                SMB.fullGrid(j,k) = nanmean(SMB.smb(j,k,Iday_start:Iday_end))./rho_i*12; % m a^-1
            end    
        end  
        % Extrapolate mean annual SF along Crane centerline
        for j=1:length(x)
            SMB.extrap(i,j) = SMB.fullGrid(RACMOy(j),RACMOx(j));
        end        

        % Adjust SMB for elevation-dependence

        % Loop through each point along centerline
        for j=1:length(x)

            % Grab points within a certain distance from centerline, interpolate
            nearbyPts = []; % Hold points within a certain distance
            numPts = 0; % Total number of points found

            for k=1:length(h.X(:,1))
                for l=1:length(h.X(1,:))
                    dist = sqrt((cl.X(j)-h.X(k,l))^2 + (cl.Y(j)-h.Y(k,l))^2);
                    if dist<=maxDist
                        numPts = numPts+1;
                        nearbyPts(numPts,1:2) = ([h.h(k,l) SMB.fullGrid(k,l)]); 
                    end 
                end 
            end 

            % Grab the mean of nearby points
            SMB.mean(i,j) = nanmean(nearbyPts(:,2));
            
            %Calculate a linear trendline for nearby points
            P = polyfit(nearbyPts(:,1),nearbyPts(:,2),1);
            int = -500:2000; 
            yfit = P(1)*int+P(2);
    
            SMB.downscaled(i,j) = interp1(int,yfit,h_cl(j));
    
        end     
    
        % Fit a linear trendline to SMB
        SMB.linear(i,~isnan(h_cl)) = feval(fit(h_cl(~isnan(h_cl)),SMB.downscaled(i,~isnan(h_cl))','poly1'),h_cl(~isnan(h_cl)))';
        SMB.linear(i,isnan(h_cl)) = SMB.linear(i,find(isnan(h_cl),1,'first')-1);
    
        % Plot
        figure(1); subplot(2,3,3);
        plot(x/10^3,SMB.downscaled(i,:),'--','color',col(i,:),'linewidth',1,'displayname',num2str(years(i))); 
        plot(x/10^3,SMB.linear(i,:),'-','color',col(i,:),'linewidth',2,'HandleVisibility','off');     
        drawnow

        % Calculate and plot the mean SM at each point 
        if years(i)==2016
            SMB.downscaled_average = nanmean(SMB.downscaled,1);
            SMB.downscaled_average_linear = feval(fit(h_cl(~isnan(h_cl)),SMB.downscaled_average(~isnan(h_cl))','poly1'),h_cl);
            plot(x/10^3,SMB.downscaled_average_linear,'-k','linewidth',2);
        end

    end   

    % RUNOFF (RO)
    if years(i)<=2019
    
        % estimate using the difference between the native 
        % and downscaled resolution snowmelt
        if years(i)<2016
            RO.diff(i,:) = SM.downscaled(i,:)-SM.extrap(i,:);
        else
            RO.diff(i,:) = SMB.downscaled(i,:)-SMB.extrap(i,:);
        end
        
        % calculate linear trendline
        RO.linear(i,:) = feval(fit(x(~isnan(RO.diff(i,:)))',RO.diff(i,~isnan(RO.diff(i,:)))','poly1'),x')';
        
        % Plot
        figure(1);
        subplot(2,3,[4,4.5]);    
        plot(x/10^3,RO.diff(i,:),'--','color',col(i,:),'linewidth',1);
        plot(x/10^3,RO.linear(i,:),'-','color',col(i,:),'linewidth',2);

        % Calculate and plot the mean RO at each point 
        if years(i)==2016
            RO.downscaled_average = nanmean(RO.diff,1);
            RO.downscaled_average_linear = feval(fit(h_cl(~isnan(h_cl)),RO.downscaled_average(~isnan(h_cl))','poly1'),h_cl);
            plot(x/10^3,RO.downscaled_average_linear,'-k','linewidth',2);
        end

    end

    % AIR TEMPERATURE (T)
    T.extrap = NaN*zeros(length(years),length(x)); 
    T.mean = NaN*zeros(length(years),length(x));
    T.downscaled = NaN*zeros(length(years),length(x));
    if years(i)<=2019
 
        T.fullGrid = NaN*zeros(length(T.T(:,1,1)),length(T.T(1,:,1))); % initialize
        for j=1:length(T.T(:,1,1))
            for k=1:length(T.T(1,:,1))
                T.fullGrid(j,k) = nanmean(T.T(j,k,Iday_start:Iday_end))-273.15; % ^oC
            end     
        end 
        % Extrapolate mean annual T along centerline (RACMO days day_start:day_end)
        for j=1:length(RACMOy)
            T.extrap(i,j) = T.fullGrid(RACMOy(j),RACMOx(j));
        end  
   
        % Loop through each point along centerline
        for j=1:length(x)
        
            % Grab points within a certain distance from centerline, interpolate
            maxDist = 9e3; 
            nearbyPts = []; % Hold points within a certain distance
            numPts = 0; % Total number of points found
        
            for k=1:length(h.X(:,1))
                for l=1:length(h.X(1,:))
                    dist = sqrt((cl.X(j)-h.X(k,l))^2 + (cl.Y(j)-h.Y(k,l))^2);
                    if dist<=maxDist
                        numPts = numPts+1;
                        nearbyPts(numPts,1:2) = ([h.h(k,l) T.fullGrid(k,l)]); 
                    end 
                end 
            end 
        
            % Grab the mean of nearby points
            T.mean(i,j) = nanmean(nearbyPts(:,2));
           
            % Calculate a linear trendline for nearby points
            P = polyfit(nearbyPts(:,1),nearbyPts(:,2),1);
            int = -10:1800; 
            yfit = P(1)*int+P(2);
        
            % Apply regression slope to current grid cell, solve for eqn
            m = P(1); b = T.extrap(i,j)-m*h.cl(j);
            T_fit = [h.cl m.*h.cl+b];
            T.downscaled(i,j) = m*h_cl(j)+b; 

        end 
        
        % Dry adiabatic lapse rate    
        lr = 9.8e-3; % lapse rate (degrees C m^-1)
        
        % Calculate new air temperature using lapse rate and RACMO reference height
        T_cl_m = movmean(T.downscaled(i,:),10,'omitnan');
        h_cl_m = movmean(h_cl,10,'omitnan');
        T.adjusted(i,:) = (h.cl-h_cl)'*lr+T.extrap(i,:); 
    
        % Plot
        figure(1);
        subplot(2,3,[5.5,6]);
        plot(x/10^3,T.adjusted(i,:),'color',col(i,:),'linewidth',2);

        % Calculate and plot the mean SM at each point 
        if years(i)==2016
            T.adjusted_average = nanmean(T.adjusted,1);
            T.adjusted_average_linear = feval(fit(h_cl(~isnan(h_cl)),T.adjusted_average(~isnan(h_cl))','poly1'),h_cl);
            plot(x/10^3,T.adjusted_average_linear,'-k','linewidth',2);
        end

    end

end 
    
% save figure
if figure_save
    cd([homepath,'figures/']);
    saveas(gcf,'adjustedRACMOVariables_2009-2019.png','png');
    disp('Figure 1 saved');
end 

% save downscaled SMB and Runoff
if save_variables
    cd([homepath,'inputs-outputs/']);
    save('downscaledClimateVariables_2009-2019.mat','RO','SMB','T');
    disp('Climate variables saved');
end 



