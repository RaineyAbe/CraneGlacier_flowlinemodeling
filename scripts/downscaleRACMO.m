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
homepath = '/Users/raineyaberle/Desktop/Research/CraneGlacier_flowlinemodeling/';
cd(homepath);

% Add path with necessary functions, data, inputs/outputs
addpath([homepath,'matlabFunctions/']);
addpath([homepath,'../matlabFunctions/cmocean_v2.0/cmocean/']);
addpath([homepath,'data/RACMO2.3/']); 
addpath([homepath,'inputs-outputs/']);

save_smb = 0; % = 1 to save resulting SMB
save_temp = 0; % = 1 to save resulting air temperature

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
yr = 2016;
day_start = (yr-1950)*365 + 1;
day_end = (yr-1950)*365 + 365;
        
% 2011 SNOWFALL (sf)
    sf.Lat = ncread('RACMO2.3p2_XPEN055_snowfall_monthly_1979_2016.nc','lat'); % degrees north
    sf.Lon = ncread('RACMO2.3p2_XPEN055_snowfall_monthly_1979_2016.nc','lon'); % degrees east
    sf.h = ncread('RACMO2.3p2_XPEN055_snowfall_monthly_1979_2016.nc','height'); % m above surface   
    sf.sf = ncread('RACMO2.3p2_XPEN055_snowfall_monthly_1979_2016.nc','snowfall'); % kg m-2 d-1 
    time = ncread('RACMO2.3p2_ANT27_snowfall_monthly_1979_2016.nc','time'); % days since 1950/01/01
    Iday_start = dsearchn(time,day_start); % index of day start in RACMO time
    Iday_end = dsearchn(time,day_end); % index of day start in RACMO time
    
    % Grab mean sf for every month in year (RACMO days mo_start:mo_end)
    sf_yr = zeros(length(sf.sf(:,1,1)),length(sf.sf(1,:,1)),1); % initialize
    for i = 1:length(sf.sf(:,1,1))
        for j=1:length(sf.sf(1,:,1))
            sf_yr(i,j,1) = nanmean(sf.sf(i,j,Iday_start:Iday_end)); % kg/m^2/d
        end     
    end 
    % kg/m^2/d / (917 kg/m^3) * (365 d/yr) = m/yr
    sf_yr = sf_yr./rho_i.*365; % m/yr
        
    % Interpolate mean annual sf along Crane centerline (RACMO days day_start:day_end)
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
    sf.cl(:,2:end) = []; % kg m-2 day-1
    sf.cl = sf.cl.*365./rho_i; % m yr-1
        
    figure(1); clf; hold on; 
        colormap(cmocean('algae'));
        set(gcf,'Units','centimeters','Position',[5 20 18 13]);
        surf(sf.Lon,sf.Lat,sf_yr); view(2);
        plot3(cl.Lon,cl.Lat,ones(length(cl.Lon))*3000,'-m','LineWidth',4,'DisplayName','Crane Centerline');
        set(gca,'FontName','Arial','FontSize',14,'linewidth',2);
        xlabel('Lon'); ylabel('Lat'); title(['RACMO Mean ',num2str(yr),' Snowfall']); 
        h = colorbar; set(get(h,'title'),'string','m yr^-^1');
        ylim([-66 -65]); xlim([-63.5 -61.5]); caxis([0 100]);
        hold off;
    
% 2011 SNOWMELT (sm)
    sm.Lat = ncread('RACMO2.3p2_XPEN055_snowmelt_monthly_1979_2016.nc','lat'); % degrees north
    sm.Lon = ncread('RACMO2.3p2_XPEN055_snowmelt_monthly_1979_2016.nc','lon'); % degrees east
    sm.h = ncread('RACMO2.3p2_XPEN055_snowmelt_monthly_1979_2016.nc','height'); % m above surface   
    sm.sm = squeeze(ncread('RACMO2.3p2_XPEN055_snowmelt_monthly_1979_2016.nc','snowmelt')); % kg m-2 d-1
    
     %Grab mean sm for every month in year (RACMO days day_start:day_end)
    sm_yr = zeros(length(sm.sm(:,1,1)),length(sm.sm(1,:,1))); % initialize
    for i = 1:length(sm.sm(:,1,1))
        for j=1:length(sm.sm(1,:,1))
            sm_yr(i,j,1) = nanmean(sm.sm(i,j,Iday_start:Iday_end));
        end     
    end 
    sm_yr = sm_yr./rho_i.*365; %m yr^-1

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
    sm.cl(:,2:end) = []; % kg m-2 day-1
    sm.cl = sm.cl.*365./rho_i; % m yr-1
        
    figure(2); clf; hold on; set(gcf,'Units','centimeters','Position',[5 0 18 13]);
        colormap(cmocean('algae'));
        surf(sm.Lon,sm.Lat,sm_yr); view(2);
        plot3(cl.Lon,cl.Lat,ones(length(cl.Lon))*3000,'-m','LineWidth',4,'DisplayName','Crane Centerline');
        set(gca,'FontName','Arial','FontSize',14,'linewidth',2);
        xlabel('Lon'); ylabel('Lat'); title(['RACMO Mean ',num2str(yr),' Snowmelt']);  
        h = colorbar; set(get(h,'title'),'string','m yr^-^1');
        ylim([-66 -65]); xlim([-63.5 -61.5]); caxis([0 100]);
        hold off;
    
% 2011 RUNOFF (ro)
    ro.Lat = ncread('RACMO2.3p2_XPEN055_runoff_monthly_1979_2016.nc','lat'); % degrees north
    ro.Lon = ncread('RACMO2.3p2_XPEN055_runoff_monthly_1979_2016.nc','lon'); % degrees east
    ro.h = ncread('RACMO2.3p2_XPEN055_runoff_monthly_1979_2016.nc','height'); % m above surface   
    ro.ro = squeeze(ncread('RACMO2.3p2_XPEN055_runoff_monthly_1979_2016.nc','runoff')); % kg m-2 d-1
    
    % Grab mean annual ro for full RACMO grid (RACMO days day_start:day_end)
    ro_yr = zeros(length(ro.ro(:,1,1)),length(ro.ro(1,:,1))); % initialize
    for i = 1:length(ro.ro(:,1,1))
        for j=1:length(ro.ro(1,:,1))
            ro_yr(i,j,1) = nanmean(ro.ro(i,j,Iday_start:Iday_end));
        end     
    end 
    ro_yr = ro_yr./rho_i.*365; %m yr^-1

    % Interpolate mean annual ro along Crane centerline (RACMO days day_start:day_end)
    ro.cl=zeros(length(RACMOy),length(Iday_start:Iday_end)); % initialize
    for i=1:length(RACMOy)
        for j=Iday_start:Iday_end
            ro.cl(i,j) = ro.ro(RACMOy(i),RACMOx(i),j);
        end    
    end  

    % Take the average at each point along centerline
    ro.cl(:,1:Iday_start-1)=[]; % Started adding points in column mo_start
    ro.cl(:,1) = nanmean(ro.cl,2);
    ro.cl(:,2:end) = []; % kg m-2 day-1
    ro.cl = ro.cl.*365./rho_i; % m yr-1

    figure(3); clf; hold on; set(gcf,'Units','centimeters','Position',[25 20 18 13]);
        colormap(cmocean('algae'));        
        surf(ro.Lon,ro.Lat,ro_yr); view(2);
        plot3(cl.Lon,cl.Lat,ones(length(cl.Lon))*3000,'-m','LineWidth',4,'DisplayName','Crane Centerline');
        set(gca,'FontName','Arial','FontSize',14,'linewidth',2);
        xlabel('Lon'); ylabel('Lat'); title(['RACMO Mean ',num2str(yr),' Runoff']);  
        h = colorbar; set(get(h,'title'),'string','m yr^-^1');
        ylim([-66 -65]); xlim([-63.5 -61.5]); caxis([0 100]);
        hold off;    

% SURFACE MASS BALANCE (SMB)
    smb.Lat = ncread('SMB_RACMO2.3p2_monthly_XPEN055_197901_201908.nc','lat'); %degrees north
    smb.Lon = ncread('SMB_RACMO2.3p2_monthly_XPEN055_197901_201908.nc','lon'); %degrees east
    smb.h = ncread('SMB_RACMO2.3p2_monthly_XPEN055_197901_201908.nc','height'); %m above the surface   
    smb.smb = squeeze(ncread('SMB_RACMO2.3p2_monthly_XPEN055_197901_201908.nc','smb')); %kg m-2 d-1
    
    % Grab mean smb for every month in year (RACMO days Iday_start:Iday_end)
    smb_yr = zeros(length(smb.smb(:,1,1)),length(smb.smb(1,:,1))); % initialize
    for i = 1:length(smb.smb(:,1,1))
        for j=1:length(smb.smb(1,:,1))
            smb_yr(i,j,1) = nanmean(smb.smb(i,j,Iday_start:Iday_end));
        end     
    end 
    smb_yr = smb_yr./rho_i.*365; %m yr^-1

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
    smb.cl(:,2:end) = []; % kg m-2 day-1
    smb.cl = smb.cl.*365./rho_i; % m yr-1

    figure(4); clf; hold on; set(gcf,'Units','centimeters','Position',[25 0 18 13]);
        colormap(cmocean('algae'));
        surf(smb.Lon,smb.Lat,smb_yr); view(2);
        plot3(cl.Lon,cl.Lat,ones(length(cl.Lon))*3000,'-m','LineWidth',4,'DisplayName','Crane Centerline');
        set(gca,'FontName','Arial','FontSize',14);
        xlabel('Lon'); ylabel('Lat'); title(['RACMO Mean ',num2str(yr),' SMB']);  
        h = colorbar; set(get(h,'title'),'string','m/yr');
        ylim([-66 -65]); xlim([-63.5 -61.5]); caxis([0 100]);
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
        
%% 3. Adjust SF, SM, & SMB  along centerline 
%   for elevation dependence, adapted from Noel et al.(2016)

close all;

%Convert RACMO grid to polar stereographic
[h.X,h.Y] = wgs2ps(h.Lon,h.Lat,'StandardParallel',-71,'StandardMeridian',0);
[sf.X,sf.Y] = wgs2ps(sf.Lon,sf.Lat,'StandardParallel',-71,'StandardMeridian',0);
[sm.X,sm.Y] = wgs2ps(sm.Lon,sm.Lat,'StandardParallel',-71,'StandardMeridian',0);
[smb.X,smb.Y] = wgs2ps(smb.Lon,smb.Lat,'StandardParallel',-71,'StandardMeridian',0);

%Load Crane ice surface
h_cl = load("Crane_surfaceElevationObs.mat").h(28).surface';

maxDist = 9e3; 
   
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
    save_smb('Crane_downscaledSMB_2011.mat','x','h_cl_11','smb.interp','smb.linear');
    save_smb('Crane_adjustedAirTemp_2011.mat','T.cl_11');
    disp('downscaled SMB and air temperature saved for one year.');
end 

%% 6. Repeat above steps to load smb for 2009-2019

close all; 

figure_save = 0;    % = 1 to save figure
save_smb = 0;       % = 1 to save final downscaled SMB

years = 2009:2019; % Define years
col = parula(length(years)+1); % color scheme for plotting

% Load centerline
cl.X = load('Crane_centerline.mat','x').x; 
cl.Y = load('Crane_centerline.mat','y').y;
cd([homepath,'inputs-outputs/']);
    
    % Load most advanced terminus position (2019)
    term = dsearchn([cl.X cl.Y],[load('Crane_terminusPositions_2002-2019.mat').term(61).X ...
    load('Crane_terminusPositions_2002-2019.mat').term(61).Y]); 
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
    sf.Lat = ncread('RACMO2.3p2_XPEN055_snowfall_monthly_1979_2016.nc','lat'); % degrees north
    sf.Lon = ncread('RACMO2.3p2_XPEN055_snowfall_monthly_1979_2016.nc','lon'); % degrees east
    sf.sf = squeeze(ncread('RACMO2.3p2_XPEN055_snowfall_monthly_1979_2016.nc','snowfall')); % kg m-2 d-1 
    time = ncread('RACMO2.3p2_XPEN055_snowfall_monthly_1979_2016.nc','time'); % days since 1950/01/01

    sm.Lat = ncread('RACMO2.3p2_XPEN055_snowmelt_monthly_1979_2016.nc','lat'); % degrees north
    sm.Lon = ncread('RACMO2.3p2_XPEN055_snowmelt_monthly_1979_2016.nc','lon'); % degrees east
    sm.sm = squeeze(ncread('RACMO2.3p2_XPEN055_snowmelt_monthly_1979_2016.nc','snowmelt')); % kg m-2 d-1

    ro.Lat = ncread('RACMO2.3p2_XPEN055_runoff_monthly_1979_2016.nc','lat'); % degrees north
    ro.Lon = ncread('RACMO2.3p2_XPEN055_runoff_monthly_1979_2016.nc','lon'); % degrees east
    ro.ro = squeeze(ncread('RACMO2.3p2_XPEN055_runoff_monthly_1979_2016.nc','runoff')); % kg m-2 d-1
   
    smb.Lat = ncread('SMB_RACMO2.3p2_monthly_XPEN055_197901_201908.nc','lat'); % degrees north
    smb.Lon = ncread('SMB_RACMO2.3p2_monthly_XPEN055_197901_201908.nc','lon'); % degrees east
    smb.smb = squeeze(ncread('SMB_RACMO2.3p2_monthly_XPEN055_197901_201908.nc','smb')); % kg m-2 d-1
    
    h.Lat = ncread('Height_latlon_XPEN055.nc','lat'); % degrees north
    h.Lon = ncread('Height_latlon_XPEN055.nc','lon'); % degrees east
    h.h = ncread('Height_latlon_XPEN055.nc','height'); % m above the surface
    h_cl_RACMO = griddata(h.Lon,h.Lat,h.h,cl.Lon,cl.Lat); 
    
    T_Lat = ncread('RACMO2.3p2_XPEN055_T2m_monthly_1979_2016.nc','lat'); % degrees north
    T_Lon = ncread('RACMO2.3p2_XPEN055_T2m_monthly_1979_2016.nc','lon'); % degrees east
    T_T = squeeze(ncread('RACMO2.3p2_XPEN055_T2m_monthly_1979_2016.nc','t2m')); % degrees K

    % Grab RACMO grid
        % Note: RACMO day # = (year - 1950)*365 + (day # in year)
        % Ex: 1 Feb/2000 = (2000-1950)*365 + 1 = 18251
        for i=1:length(cl.Lon)
            lat_diff = abs(cl.Lat(i)*ones(size(sf.Lat)) - sf.Lat);
            lon_diff = abs(cl.Lat(i,1)*ones(size(sf.Lon)) - sf.Lon);
            diff_map = sqrt(lat_diff.^2+lon_diff.^2);
            RACMO_ref = find(diff_map==min(min(diff_map)));
            [RACMOy(i),RACMOx(i)] = ind2sub(size(nanmean(sf.sf(:,:,1),4)),RACMO_ref);
        end     
        
    % Convert RACMO grid to polar stereographic
    [h.X,h.Y] = wgs2ps(h.Lon,h.Lat,'StandardParallel',-71,'StandardMeridian',0);
    [sf.X,sf.Y] = wgs2ps(sf.Lon,sf.Lat,'StandardParallel',-71,'StandardMeridian',0);
    [sm.X,sm.Y] = wgs2ps(sm.Lon,sm.Lat,'StandardParallel',-71,'StandardMeridian',0);
    [smb.X,smb.Y] = wgs2ps(smb.Lon,smb.Lat,'StandardParallel',-71,'StandardMeridian',0);
    
    % Load 2016 Crane ice surface
    cd([homepath,'inputs-outputs/']);
    h_cl = load("Crane_surfaceElevationObs.mat").h(28).surface';
    
    maxDist = 8e3;     
    
        % Grab points within a certain distance from centerline, interpolate
        nearbyPts = []; % Hold points within a certain distance
        numPts = 0; % Total number of points found    

for i=1:length(years)

% Calculate Adjusted Snowfall (sf), Snowmelt (sm), 
% and Surface Mass Balance (SMB) along centerline

    day_start = (years(i)-1950)*365+1; Iday_start = dsearchn(time,day_start);
    day_end = (years(i)-1950)*365+365; Iday_end = dsearchn(time,day_end);
    
    % set up figures
    if i==1
        figure(1); clf
        set(gcf,'units','centimeters','position',[5 0 40 30]);
        subplot(2,3,1); hold on;  % adjusted snowfall
            set(gca,'FontName','Arial','FontSize',14);
            xlabel('Elevation (m)'); ylabel('(m a^{-1})'); 
            title('Adjusted Snowfall');grid on;
        subplot(2,3,2); hold on; % adjusted snowmelt
            set(gca,'FontName','Arial','FontSize',14);
            xlabel('Elevation (m)'); ylabel('(m a^{-1})'); 
            title('Adjusted Snowmelt');grid on; 
        subplot(2,3,3); hold on; % adjusted SMB
            set(gca,'FontName','Arial','FontSize',14);
            xlabel('Elevation (m)'); ylabel('(m a^{-1})'); 
            title('Adjusted SMB');grid on;         
        subplot(2,3,[4,4.5]); hold on; % linearized SMB
            set(gca,'FontName','Arial','FontSize',14);
            xlabel('Distance Along Centerline (km)'); ylabel('(m a^{-1})'); 
            title('Linearized SMB');grid on; 
        subplot(2,3,[5.5,6]); hold on; % adjusted air temp
            set(gca,'FontName','Arial','FontSize',14);
            xlabel('Distance Along Centerline (km)'); ylabel('(^oC)'); 
            title('Adjusted Air Temp'); grid on; 
        c = colorbar('position',[0.93 0.4 0.02 0.3],'fontsize',14);
        c.Ticks = [0 0.5 1.0]; c.TickLabels = {'2009','2014','2019'};
    end 

% SNOWFALL (sf)
if years(i)<=2016
    % Grab mean sf for every month for full RACMO grid
        sf.Yr = zeros(length(sf.sf(:,1,1)),length((sf.sf(1,:,1))));
        for j = 1:length(sf.sf(:,1,1))
            for k=1:length(sf.sf(1,:,1))
                sf.Yr(j,k) = nanmean(sf.sf(j,k,Iday_start:Iday_end)); % kg/m^2/day
            end     
        end 
        sf.Yr = sf.Yr./rho_i.*365; % m/yr

    % Interpolate mean annual sf along Crane centerline
        sf.cl=[];
        for j=1:length(RACMOy)
            for k=Iday_start:Iday_end
                sf.cl(j,k) = sf.sf(RACMOy(j),RACMOx(j),k);
            end    
        end  

       % Take the average at each point along centerline
        sf.cl(:,1:(Iday_start-1))=[]; % Delete empty columns
        sf.cl(:,1) = nanmean(sf.cl,2);
        sf.cl(:,2:end) = []; % kg m-2 day-1
        sf.cl = sf.cl.*365./rho_i; % m a-1

   % Adjust snowfall (sf) for elevation-dependence

    % Loop through each point along centerline
    sf.interp = ones(length(h_cl),1);
    for j=1:length(x)

        % Grab points within a certain distance from centerline, interpolate
        nearbyPts = []; % Hold points within a certain distance
        numPts = 0; % Total number of points found

        xi = cl.X(j); yi = cl.Y(j);
        for k=1:length(h.X(:,1))
            for l=1:length(h.X(1,:))
                dist = sqrt((xi-h.X(k,l))^2 + (yi-h.Y(k,l))^2);
                if dist<=maxDist
                    numPts = numPts+1;
                    nearbyPts(numPts,1:2) = ([h.h(k,l) sf.Yr(k,l)]); 
                end 
            end 
        end 

        %Calculate a linear trendline for nearby points
        P = polyfit(nearbyPts(:,1),nearbyPts(:,2),1);
        int = -500:2000; 
        yfit = P(1)*int+P(2);

        sf.interp(j) = interp1(int,yfit,h_cl(j));

    end     
    
    % Fit a linear trendline to snowfall
    sf.linear(i,1:length(find(~isnan(h_cl)))) = feval(fit(h_cl(~isnan(h_cl)),sf.interp(~isnan(h_cl)),'poly1'),h_cl(~isnan(h_cl)))';
    sf.linear(i,length(find(~isnan(h_cl)))+1:length(x)) = sf.linear(end);

    figure(1); subplot(2,3,1);
    plot(h_cl,sf.interp,'--','color',col(i,:),'linewidth',1,'displayname',num2str(years(i))); 
    plot(h_cl,sf.linear(i,:),'-','color',col(i,:),'linewidth',2,'HandleVisibility','off');     
    drawnow
    
end   

%SNOWMELT (sm)
if years(i)<=2016
    %Grab mean sm for every month for full RACMO grid
        sm.Yr = zeros(length(sm.sm(:,1,1)),length(sm.sm(1,:,1)));
        for j = 1:length(sm.sm(:,1,1))
            for k=1:length(sm.sm(1,:,1))
                sm.Yr(j,k,1) = nanmean(sm.sm(j,k,Iday_start:Iday_end)); % kg m^2 / day
            end     
        end 
        sm.Yr = sm.Yr./rho_i.*365; %m a^-1
        
    %Interpolate mean annual sm along Crane centerline
        sm.cl=zeros(length(RACMOy),length(Iday_start:Iday_end));
        for j=1:length(RACMOy)
            for k=Iday_start:Iday_end
                sm.cl(j,k) = sm.sm(RACMOy(j),RACMOx(j),k);
            end    
        end  
        
       %Take the average at each point along centerline
        sm.cl(:,1:(Iday_start-1))=[]; %Delete empty columns
        sm.cl(:,1) = nanmean(sm.cl,2);
        sm.cl(:,2:end) = []; %kg m-2 day-1
        sm.cl = sm.cl.*365./rho_i; %m a-1
        
    %Ajust snowmelt for elevation-dependence
    
        %Loop through each point along centerline
        sm.interp = ones(length(h_cl),1);
        for j=1:length(x)
    
            %Grab points within a certain distance from centerline, interpolate
            nearbyPts = []; %Hold points within a certain distance
            numPts = 0; %Total number of points found
    
            xi = cl.X(j); yi = cl.Y(j);
            for k=1:length(h.X(:,1))
                for l=1:length(h.X(1,:))
                    dist = sqrt((xi-h.X(k,l))^2 + (yi-h.Y(k,l))^2);
                    if dist<=maxDist
                        numPts = numPts+1;
                        nearbyPts(numPts,1:2) = ([h.h(k,l) sm.Yr(k,l)]); 
                    end 
                end 
            end 
        
            %Calculate a linear trendline for nearby points
            P = polyfit(nearbyPts(:,1),nearbyPts(:,2),1);
            int = -20:1600; 
            yfit = P(1)*int+P(2);
            sm.interp(j) = interp1(int,yfit,h_cl(j));
        
        end 
    
    % Fit a linear trendline to snowmelt
    sm.linear(i,1:length(find(~isnan(h_cl)))) = feval(fit(h_cl(~isnan(h_cl)),sm.interp(~isnan(h_cl)),'poly1'),h_cl(~isnan(h_cl)))';
    sm.linear(i,length(find(~isnan(h_cl)))+1:length(x)) = sm.linear(end);

    figure(1); subplot(2,3,2);
    plot(h_cl,sm.interp,'--','color',col(i,:),'linewidth',1,'displayname',num2str(years(i))); 
    plot(h_cl,sm.linear(i,:),'-','color',col(i,:),'linewidth',2,'HandleVisibility','off');     
    drawnow
    
end

%SURFACE MASS BALANCE (smb)

    %Grab mean smb for every month in year
        smb.Yr = zeros(length(smb.smb(:,1,1)),length(smb.smb(1,:,1)));
        for j = 1:length(smb.smb(:,1,1))
            for k=1:length(smb.smb(1,:,1))
                smb.Yr(j,k,1) = nanmean(smb.smb(j,k,Iday_start:Iday_end)); % kg a^-1
            end     
        end 
        smb.Yr = smb.Yr./rho_i.*365;  % m a^-1
            
        %Interpolate mean annual smb along Crane centerline (RACMO days day_start:day_end)
        smb.cl=zeros(length(RACMOy),length(Iday_start:Iday_end));
        for j=1:length(RACMOy)
            for k=Iday_start:Iday_end
                smb.cl(j,k) = smb.smb(RACMOy(j),RACMOx(j),k);
            end    
        end  
        
       %Take the average at each point along centerline
        smb.cl(:,1:(Iday_start-1))=[]; %Delete empty columns
        smb.cl(:,1) = nanmean(smb.cl,2);
        smb.cl(:,2:end) = []; %kg m-2 day-1
        smb.cl = smb.cl.*365./rho_i; %m a-1    

    %Ajust surface mass balance for elevation-dependence

    %Loop through each point along centerline
    smb.interp = ones(length(h_cl),1);
    for j=1:length(x)
    
        %Grab points within a certain distance from centerline, interpolate
        nearbyPts = []; %Hold points within a certain distance
        numPts = 0; %Total number of points found
    
        xi = cl.X(j); yi = cl.Y(j);
        for k=1:length(h.X(:,1))
            for l=1:length(h.X(1,:))
                dist = sqrt((xi-h.X(k,l))^2 + (yi-h.Y(k,l))^2);
                if dist<=maxDist
                    numPts = numPts+1;
                    nearbyPts(numPts,1:2) = ([h.h(k,l) smb.Yr(k,l)]); 
                end 
            end 
        end 
        
        %Calculate a linear trendline for nearby points
        smb.interp(j) = polyval(polyfit(nearbyPts(:,1),...
            nearbyPts(:,2),1),h_cl(j));
       
    end 
    
    % Fit a linear trendline to snowfall
    smb.linear(i,1:length(find(~isnan(h_cl)))) = feval(fit(h_cl(~isnan(h_cl)),smb.interp(~isnan(h_cl)),'poly1'),h_cl(~isnan(h_cl)))';
    smb.linear(i,length(find(~isnan(h_cl)))+1:length(x)) = smb.linear(end);

    figure(1); subplot(2,3,3);
    plot(h_cl,smb.interp,'--','color',col(i,:),'linewidth',1,'displayname',num2str(years(i))); 
    plot(h_cl,smb.linear(i,:),'-','color',col(i,:),'linewidth',2,'HandleVisibility','off');     
    drawnow
    
    % Save in structure
    SMB(i).year = years(i); % Year
    sigma_h_cl_11 = 25; % uncertainty in ice surface (m) - two points
    SMB(i).smb_interp = smb.interp; % interpolated smb (m/a)
    
    % Calculate uncertainty in final SMB vector: 
    % sigma_smb: [10% SMB raw product] &
    %           [0.49 mWE/a from mean bias for downscaled product vs.
    %           smb observed for eight transects in GrIS (Noel et al., 2016)] 
    %           & [uncertainty in ice surface elevation]
    %           = 0.1.*smb_interp + 0.49*(rho_w/rho_i*sigma_smbWE)
    % mWE/a -> m ice/a: rho_w*H_w = rho_i*H_i
    %                   ->  H_i = rho_w/rho_i/H_w
    H_w = 0.49; % mWE/a mean bias vs. observations (Noel et al., 2016)
    rho_i = 917; rho_w = 1000; % kg/m^3
    for j=1:length(x)
        SMB(i).sigma_smb(j) = SMB(i).smb_interp(j).*sqrt((0.1.*smb.interp(j)/smb.interp(j))^2 ...
        +(rho_w/rho_i*H_w/H_w)^2+(sigma_h_cl_11/h_cl(j))^2);
    end 
    
    % Calculate final SMB as sf-sm evaluated at the ice surface elevation
    if years(i)<2016
        P = fit(h_cl(~isnan(h_cl)),sf.linear(i,~isnan(h_cl))'-sm.linear(i,~isnan(h_cl))','poly1');
        SMB(i).smb_adj = P.p1*h_cl+P.p2;
    else
        P = fit(h_cl(~isnan(h_cl)),smb.linear(i,~isnan(h_cl))','poly1');        
        SMB(i).adj = P.p1*h_cl+P.p2;
        SMB(i).adj(find(isnan(SMB(i).adj),1,'first')+1:end) = SMB(i).adj(find(isnan(SMB(i).adj),1,'first')-1);
    end
    subplot(2,3,[4,4.5]);
    plot(x/10^3,SMB(i).smb_adj,'linewidth',2,'color',col(i,:));
    drawnow
    
% AIR TEMPERATURE (T)
    T(i).T_yr = zeros(length(T_T(:,1,1)),length(T_T(1,:,1))); % initialize
    for j=1:length(T_T(:,1,1))
        for k=1:length(T_T(1,:,1))
            T(i).T_yr(j,k,1) = nanmean(T_T(j,k,Iday_start:Iday_end));
        end     
    end 
    T(i).T_yr = T(i).T_yr - 273.15; % degrees C
    % Interpolate mean annual T along centerline (RACMO days day_start:day_end)
    T_cl=zeros(length(RACMOy),length(Iday_start:Iday_end));
    for j=1:length(RACMOy)
        for k=Iday_start:Iday_end
            T_cl(j,k) = T_T(RACMOy(j),RACMOx(j),k);
        end    
    end  
    % Take the average at each point along centerline
    T_cl(:,1:Iday_start-1)=[]; % started adding points in column Iday_start
    T_cl(:,1) = nanmean(T_cl,2);
    T_cl(:,2:end) = []; % deg K
    T_cl = T_cl - 273.15; % degrees C

    % Loop through each point along centerline
    T(i).interp = ones(length(h_cl),1);
    for j=1:length(x)

        % Grab points within a certain distance from centerline, interpolate
        maxDist = 9e3; 
        nearbyPts = []; % Hold points within a certain distance
        numPts = 0; % Total number of points found

        xi = cl.X(i); yi = cl.Y(i);
        for k=1:length(h.X(:,1))
            for l=1:length(h.X(1,:))
                dist = sqrt((xi-h.X(k,l))^2 + (yi-h.Y(k,l))^2);
                if dist<=maxDist
                    numPts = numPts+1;
                    nearbyPts(numPts,1:2) = ([h.h(k,l) T(i).T_yr(k,l)]); 
                end 
            end 
        end 

        % Calculate a linear trendline for nearby points
        P = polyfit(nearbyPts(:,1),nearbyPts(:,2),1);
        int = -10:1800; 
        yfit = P(1)*int+P(2);

        % Apply regression slope to current grid cell, solve for eqn
        m = P(1); b = T_cl(i)-m*h_cl_RACMO(i);
        T_fit = [h_cl_RACMO m.*h_cl_RACMO+b];
        T(i).interp(i) = m*h_cl(i)+b; 
    end 

    % Dry adiabatic lapse rate    
    lr = 9.8e-3; % lapse rate (degrees C m^-1)

    % Calculate new air temperature using lapse rate and RACMO reference height
    T_cl_m = movmean(T_cl,10,'omitnan');
    h_cl_m = movmean(h_cl,10,'omitnan');
    T(i).cl_yr = (h_cl_RACMO-h_cl)*lr+T_cl; 
    
    % Plot
    figure(1);
    subplot(2,3,[5.5,6]);
    plot(x/10^3,T(i).cl_yr,'color',col(i,:),'linewidth',2);

end 
    
% save figure
if figure_save
    cd([homepath,'figures/']);
    saveas(gcf,'Crane_adjustedRACMOVariables_2011-2019.png','png');
    disp('Figure 1 saved');
end 

% save downscaled SMB
if save_smb
    cd([homepath,'inputs-outputs/']);
    save('Crane_downscaledSMB_2009-2019.mat','SMB');
    disp('SMB saved');
end 

% save air temperature
if save_temp
    cd([homepath,'inputs-outputs/']);
    save('Crane_downscaledAnnualAirTemp.mat','T');
    disp('T saved');
end


