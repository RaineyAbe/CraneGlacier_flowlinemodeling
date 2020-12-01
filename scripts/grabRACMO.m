%Rainey Aberle
%Spring 2020
%Grab RACMO2.3 climate variables and downscale SMB for elevation-dependence
%at Crane Glacier, Antarctic Peninsula

%Order of operations:
%   1. Load Crane centerline & mean 2011 RACMO variables 
%       (sf, sm, ro, smb, height), plot
%   2. Load RACMO height, interpolate along Crane centerline
%   3. Adjust sf, sm, and smb for elevation-dependence 
%       (Method adapted from Noel et al. (2016))
%   4. Use linear trendlines of sf and sm to calculate smb along centerline
%   5. Adjust air temperature T along Crane centerline using a 
%       dry adiabatic lapse rate
%   6. Complete the above steps for 2011-2016

close all; clear all; warning off; 

%1. Load centerline, 2011 Snowfall (sf), Snowmelt (sm), Runoff (ro),
%Surface Mass Balance (SMB), Air Temperature (airtemp) along centerline

    homepath = '/Users/raineyaberle/Desktop/Research/CraneGlacier_modeling/';
    cd(homepath);
    
    % Add path with MATLAB functions (ps2wgs, etc.), data, inputs/outputs
    addpath('/Users/raineyaberle/Desktop/Research/general_matlabcodes/');
    addpath([homepath,'data/RACMO2.3']); addpath([homepath,'inputs-outputs']);
    
    %Load centerline
    cl.X = load('Crane_centerline.mat','x').x; 
    cl.Y = load('Crane_centerline.mat','y').y;

    %Grab 2011 RACMO values

    [cl.Lon,cl.Lat] = ps2wgs(cl.X,cl.Y,'StandardParallel',-71,'StandardMeridian',0);
    rho_i = 917;    %kg/m^3
    rho_w = 1000;   % kg/m^3
        
    %Plot full gridded RACMO variables for  2011 and Crane centerline
        %Note: RACMO month # = (year - 1979)*12 + (mo. # in year)
        %Ex: Feb/2000 = (2000-1979)*12 + 2 = 254
        
%2011 SNOWFALL (sf)
        sf.Lat = ncread('RACMO2.3p2_XPEN055_snowfall_monthly_1979_2016.nc','lat'); %degrees north
        sf.Lon = ncread('RACMO2.3p2_XPEN055_snowfall_monthly_1979_2016.nc','lon'); %degrees east
        sf.h = ncread('RACMO2.3p2_XPEN055_snowfall_monthly_1979_2016.nc','height'); %m above surface   
        sf.sf = ncread('RACMO2.3p2_XPEN055_snowfall_monthly_1979_2016.nc','snowfall'); %kg m-2 d-1 
    
    %Grab mean sf for every month in 2011 (RACMO Months 385:396)
        sf_11 = zeros(length(sf.sf(:,1,1)),length(sf.sf(1,:,1)),1); % initialize
        for i = 1:length(sf.sf(:,1,1))
            for j=1:length(sf.sf(1,:,1))
                sf_11(i,j,1) = nanmean(sf.sf(i,j,385:396)); % m/s
            end     
        end 
        sf_11 = sf_11./rho_i.*365; %m yr^-1
        
        
    %Interpolate 2011 mean annual sf along Crane centerline (RACMO Months 385:396)
        %Grab RACMO grid
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
            for j=385:396
                sf.cl(i,j) = sf.sf(RACMOy(i),RACMOx(i),1,j);
            end    
        end  
        
       %Take the average at each point along centerline
        sf.cl(:,1:384)=[]; %Started adding points in column 385
        sf.cl(:,1) = nanmean(sf.cl,2);
        sf.cl(:,2:end) = []; %kg m-2 day-1
        sf.cl = sf.cl.*365./rho_i; %m yr-1
        
    figure; hold on; 
        set(gcf,'Units','centimeters','Position',[5 20 18 13]);
        surf(sf.Lon,sf.Lat,sf_11); view(2);
        plot3(cl.Lon,cl.Lat,ones(length(cl.Lon))*3000,'-m','LineWidth',4,'DisplayName','Crane Centerline');
        set(gca,'FontName','Arial','FontSize',14,'linewidth',2);
        xlabel('Lon'); ylabel('Lat'); title('RACMO Mean 2011 Snowfall'); 
        h = colorbar; set(get(h,'title'),'string','m yr^-^1');
        ylim([-66 -65]); xlim([-63.5 -61.5]);
        hold off;
    
%2011 SNOWMELT (sm)
    sm.Lat = ncread('RACMO2.3p2_XPEN055_snowmelt_monthly_1979_2016.nc','lat'); %degrees north
    sm.Lon = ncread('RACMO2.3p2_XPEN055_snowmelt_monthly_1979_2016.nc','lon'); %degrees east
    sm.h = ncread('RACMO2.3p2_XPEN055_snowmelt_monthly_1979_2016.nc','height'); %m above surface   
    sm.sm = squeeze(ncread('RACMO2.3p2_XPEN055_snowmelt_monthly_1979_2016.nc','snowmelt')); %kg m-2 d-1
    
    %Grab mean sm for every month in 2011 (RACMO Months 385:396)
        sm_11 = zeros(length(sm.sm(:,1,1)),length(sm.sm(1,:,1))); % initialize
        for i = 1:length(sm.sm(:,1,1))
            for j=1:length(sm.sm(1,:,1))
                sm_11(i,j,1) = nanmean(sm.sm(i,j,385:396));
            end     
        end 
        sm_11 = sm_11./rho_i.*365; %m yr^-1
        
        %Interpolate 2011 mean annual sm along Crane centerline (RACMO Months 385:396)
        sm.cl=zeros(1,length(RACMOy));
        for i=1:length(RACMOy)
            for j=385:396
                sm.cl(i,j) = sm.sm(RACMOy(i),RACMOx(i),j);
            end    
        end  
        
       %Take the average at each point along centerline
        sm.cl(:,1:384)=[]; %Started adding points in column 385
        sm.cl(:,1) = nanmean(sm.cl,2);
        sm.cl(:,2:end) = []; %kg m-2 day-1
        sm.cl = sm.cl.*365./rho_i; %m yr-1
        
    figure; hold on; set(gcf,'Units','centimeters','Position',[5 0 18 13]);
        surf(sm.Lon,sm.Lat,sm_11); view(2);
        plot3(cl.Lon,cl.Lat,ones(length(cl.Lon))*3000,'-m','LineWidth',4,'DisplayName','Crane Centerline');
        set(gca,'FontName','Arial','FontSize',14,'linewidth',2);
        xlabel('Lon'); ylabel('Lat'); title('RACMO Mean 2011 Snowmelt'); 
        h = colorbar; set(get(h,'title'),'string','m yr^-^1');
        ylim([-66 -65]); xlim([-63.5 -61.5]);
        hold off;
    
%2011 RUNOFF (ro)
    ro.Lat = ncread('RACMO2.3p2_XPEN055_runoff_monthly_1979_2016.nc','lat'); %degrees north
    ro.Lon = ncread('RACMO2.3p2_XPEN055_runoff_monthly_1979_2016.nc','lon'); %degrees east
    ro.h = ncread('RACMO2.3p2_XPEN055_runoff_monthly_1979_2016.nc','height'); %m above surface   
    ro.ro = squeeze(ncread('RACMO2.3p2_XPEN055_runoff_monthly_1979_2016.nc','runoff')); %kg m-2 d-1
    
    %Grab 2011 mean annual ro for full RACMO grid (RACMO Months 385:396)
        ro_11 = zeros(length(ro.ro(:,1,1)),length(ro.ro(1,:,1))); % initialize
        for i = 1:length(ro.ro(:,1,1))
            for j=1:length(ro.ro(1,:,1))
                ro_11(i,j,1) = nanmean(ro.ro(i,j,385:396));
            end     
        end 
        ro_11 = ro_11./rho_i.*365; %m yr^-1
        
        %Interpolate 2011 mean annual ro along Crane centerline (RACMO Months 385:396)
        ro.cl=zeros(length(RACMOy),length(385:396)); % initialize
        for i=1:length(RACMOy)
            for j=385:396
                ro.cl(i,j) = ro.ro(RACMOy(i),RACMOx(i),j);
            end    
        end  
        
       %Take the average at each point along centerline
        ro.cl(:,1:384)=[]; %Started adding points in column 385
        ro.cl(:,1) = nanmean(ro.cl,2);
        ro.cl(:,2:end) = []; %kg m-2 day-1
        ro.cl = ro.cl.*365./rho_i; %m yr-1
    
    figure; hold on; set(gcf,'Units','centimeters','Position',[25 20 18 13]);
        surf(ro.Lon,ro.Lat,ro_11); view(2);
        plot3(cl.Lon,cl.Lat,ones(length(cl.Lon))*3000,'-m','LineWidth',4,'DisplayName','Crane Centerline');
        set(gca,'FontName','Arial','FontSize',14,'linewidth',2);
        xlabel('Lon'); ylabel('Lat'); title('RACMO Mean 2011 Runoff'); 
        h = colorbar; set(get(h,'title'),'string','m yr^-^1');
        ylim([-66 -65]); xlim([-63.5 -61.5]); 
        hold off;    

%2011 SURFACE MASS BALANCE (smb)
    smb.Lat = ncread('SMB_RACMO2.3p2_monthly_XPEN055_197901_201908.nc','lat'); %degrees north
    smb.Lon = ncread('SMB_RACMO2.3p2_monthly_XPEN055_197901_201908.nc','lon'); %degrees east
    smb.h = ncread('SMB_RACMO2.3p2_monthly_XPEN055_197901_201908.nc','height'); %m above the surface   
    smb.smb = squeeze(ncread('SMB_RACMO2.3p2_monthly_XPEN055_197901_201908.nc','smb')); %kg m-2 d-1
    
    %Grab mean smb for every month in 2011 (RACMO Months 385:396)
        smb_11 = zeros(length(smb.smb(:,1,1)),length(smb.smb(1,:,1))); % initialize
        for i = 1:length(smb.smb(:,1,1))
            for j=1:length(smb.smb(1,:,1))
                smb_11(i,j,1) = nanmean(smb.smb(i,j,385:396));
            end     
        end 
        smb_11 = smb_11./rho_i.*365; %m yr^-1
            
        %Interpolate 2011 mean annual smb along Crane centerline (RACMO Months 385:396)
        smb.cl=zeros(length(RACMOy),length(385:396)); % initialize
        for i=1:length(RACMOy)
            for j=385:396
                smb.cl(i,j) = smb.smb(RACMOy(i),RACMOx(i),j);
            end    
        end  
        
       %Take the average at each point along centerline
        smb.cl(:,1:384)=[]; %Started adding points in column 385
        smb.cl(:,1) = nanmean(smb.cl,2);
        smb.cl(:,2:end) = []; %kg m-2 day-1
        smb.cl = smb.cl.*365./rho_i; %m yr-1
    
    figure; hold on; set(gcf,'Units','centimeters','Position',[25 0 18 13]);
        surf(smb.Lon,smb.Lat,smb_11); view(2);
        plot3(cl.Lon,cl.Lat,ones(length(cl.Lon))*3000,'-m','LineWidth',4,'DisplayName','Crane Centerline');
        set(gca,'FontName','Arial','FontSize',14);
        xlabel('Lon'); ylabel('Lat'); title('RACMO Mean 2011 SMB'); 
        h = colorbar; set(get(h,'title'),'string','m yr^-^1');
        ylim([-66 -65]); xlim([-63.5 -61.5]);
        hold off;

%% 2. Load RACMO height, linearly interpolate along Crane centerline
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
    figure; hold on; 
        set(gcf,'Units','centimeters','Position',[3 8 21 18]);
        surf(h.Lon,h.Lat,h.h); view(2);
        plot3(cl.Lon(1:76),cl.Lat(1:76),h_cl_RACMO(1:76),'-m','LineWidth',4,'DisplayName','Centerline');
        set(gca,'FontName','Arial','FontSize',14,'linewidth',2);
        xlabel('Lon'); ylabel('Lat'); title('RACMO Height'); 
        c = colorbar; set(get(c,'title'),'string','(m)');
        ylim([-66 -65]); xlim([-64 -62]);
        hold off;
    
    figure; hold on; 
        set(gcf,'Units','centimeters','Position',[28 8 21 18]);
        plot(x,h_cl_RACMO,'-b','LineWidth',3); grid on; 
        set(gca,'FontName','Arial','FontSize',14,'LineWidth',2); 
        title('RACMO height at Crane Centerline');
        xlabel('Distance Along Centerline (m)'); ylabel('RACMO height (m)'); 
        hold off;
        
%% 3. Adjust SF, SM, & SMB  along the Crane centerline 
%   for elevation dependence, adapted from Noel et al.(2016)

close all;

    %Convert RACMO grid to polar stereographic
    [h.X,h.Y] = wgs2ps(h.Lon,h.Lat,'StandardParallel',-71,'StandardMeridian',0);
    [sf.X,sf.Y] = wgs2ps(sf.Lon,sf.Lat,'StandardParallel',-71,'StandardMeridian',0);
    [sm.X,sm.Y] = wgs2ps(sm.Lon,sm.Lat,'StandardParallel',-71,'StandardMeridian',0);
    [smb.X,smb.Y] = wgs2ps(smb.Lon,smb.Lat,'StandardParallel',-71,'StandardMeridian',0);
    
    %Load 2011 Crane ice surface
    h_cl_11 = load("Crane_SurfaceObservations_2009-2018.mat").h(3).surface';
    
    maxDist = 9e3; 
   
%Adjust snowfall (sf) 
    
    %Loop through each point along centerline
    sf.interp = ones(length(h_cl_11),1);
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
                    nearbyPts(numPts,1:2) = ([h.h(j,k) sf_11(j,k)]); 
                end 
            end 
        end 
        
        %Calculate a linear trendline for nearby points
        P = polyfit(nearbyPts(:,1),nearbyPts(:,2),1);
        int = 0:1600; 
        yfit = P(1)*int+P(2);
    
        %Apply regression slope to current grid cell, solve for eqn
        b = sf.cl(i)-P(1)*h_cl_RACMO(i);
        sf.interp(i) = P(1)*h_cl_11(i)+b; 
        
    end 
    
%Ajust snowmelt
    
    %Loop through each point along centerline
    sm.interp = ones(length(h_cl_11),1);
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
                    nearbyPts(numPts,1:2) = ([h.h(j,k) sm_11(j,k)]); 
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
        sm.interp(i) = m*h_cl_11(i)+b; 
        
    end 

%Ajust surface mass balance

    %Loop through each point along centerline
    smb.interp = ones(length(h_cl_11),1);
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
                    nearbyPts(numPts,1:2) = ([h.h(j,k) smb_11(j,k)]); 
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
        smb.interp(i) = m*h_cl_11(i)+b; 
        
        if i==1
            figure; hold on; grid on; set(gcf,'Units','centimeters','Position',[3 8 21 18]);
            set(gca,'FontSize',14,'FontName','Arial'); 
            title('Regression Calculation for SMB, i=1');
            plot(smbfit(:,1),smbfit(:,2),'-b','LineWidth',2,'DisplayName','a');            
            plot(yfit,'-r','LineWidth',2,'DisplayName','b'); 
            plot(h_cl_RACMO(i),smb.cl(i),'.b','MarkerSize',30,'DisplayName','Current grid cell');
            plot(nearbyPts(:,1),nearbyPts(:,2),'.r','MarkerSize',30,'DisplayName','Adjacent grid cells');
            plot(h_cl_11(i),smb.interp(i),'*b','MarkerSize',25,'DisplayName','Resulting SMB');
            xlabel('Height (m)'); ylabel('SMB (m yr^-^1)');
            legend('Location','southeast');
            hold off; 
        end 
    end 

    %Plot adjusted climate variables
    figure; hold on; set(gcf,'Units','centimeters','Position',[28 8 21 18]);
    plot(h_cl_11,sf.interp,'.b','MarkerSize',15,'DisplayName','Snowfall'); 
    plot(h_cl_11,sm.interp,'.m','MarkerSize',15,'DisplayName','Snowmelt');
    plot(h_cl_11,smb.interp,'.g','MarkerSize',15,'DisplayName','SMB');
    grid on; title('Adjusted Climate Variables');
    set(gca,'FontSize',14,'FontName','Arial');
    xlabel('Height (m)'); legend; ylabel('(m yr^-^1)');
    hold off;
    
%% 4. Use linear trendline to calculate smb slope
    % Previously saved 'Crane_downscaledSMB_2011.mat' using smb.linear2

    close all;

    %Linear trend in adjusted climate variables
    sf.linear = fit(h_cl_11(~isnan(h_cl_11)),sf.interp(~isnan(h_cl_11)),'poly1'); 
        sf.linear = feval(sf.linear,h_cl_11);
    sm.linear = fit(h_cl_11(~isnan(h_cl_11)),sm.interp(~isnan(h_cl_11)),'poly1');
        sm.linear = feval(sm.linear,h_cl_11);
    smb.linear = fit(h_cl_11(~isnan(h_cl_11)),smb.interp(~isnan(h_cl_11)),'poly1');
        smb.linear = feval(smb.linear,h_cl_11);
    smb.linear2 = sf.linear-sm.linear;
    
    %Plot Results
    figure; hold on; 
        grid on; set(gcf,'Units','centimeters','Position',[3 8 21 18]);
        set(gca,'FontName','Arial','FontSize',14,'linewidth',2); 
        plot(h_cl_11,sf.linear,'.c','MarkerSize',12,'DisplayName','Snowfall');
        plot(h_cl_11,sm.linear,'.m','MarkerSize',12,'DisplayName','Snowmelt');
        plot(h_cl_11,smb.linear,'.b','MarkerSize',12,'DisplayName','SMB');
        plot(h_cl_11,smb.linear2,'.-r','MarkerSize',10,'DisplayName','SMB (SF-SM)');
        xlabel('Height (m)'); ylabel('m yr^-^1'); title('Linear Trends in Adjusted Climate Variables');
        legend('Location','east');
        hold off;
    
    figure; hold on; 
        grid on; set(gcf,'Units','centimeters','Position',[28 8 21 18]);
        set(gca,'FontName','Arial','FontSize',14,'linewidth',2); 
        plot(x,smb.interp,'x-b','MarkerSize',14,'DisplayName','Adjusted');
        plot(x,smb.linear,'--b','LineWidth',2,'DisplayName','Adjusted Linear Trend');
        plot(x,smb.cl,'-b','LineWidth',2,'DisplayName','UN-Adjusted');
        plot(x,smb.linear2,'.-r','MarkerSize',10,'DisplayName','SMB (SF-SM)');    
        title('SMB Along Crane Centerline'); legend('Location','southwest');
        xlabel('Distance Along Centerline (m)'); ylabel('SMB (m yr^-^1)');
        hold off;  
    
%% Adjust air temperature at ice surface using a dry adiabatic lapse rate
%close all;

    %Load RACMO air temperature
    cd /Users/raineyaberle/Desktop/Research/Crane_modeling/RACMO2.3
        T.Lat = ncread('RACMO2.3p2_XPEN055_T2m_monthly_1979_2016.nc','lat'); %degrees north
        T.Lon = ncread('RACMO2.3p2_XPEN055_T2m_monthly_1979_2016.nc','lon'); %degrees east
        T.h = ncread('RACMO2.3p2_XPEN055_T2m_monthly_1979_2016.nc','height'); %m above height variable   
        T.T = squeeze(ncread('RACMO2.3p2_XPEN055_T2m_monthly_1979_2016.nc','t2m')); % degrees K
        
        T.T_11 = zeros(length(T.T(:,1,1)),length(T.T(1,:,1))); % initialize
        for i = 1:length(T.T(:,1,1))
            for j=1:length(T.T(1,:,1))
                T.T_11(i,j,1) = nanmean(T.T(i,j,385:396));
            end     
        end 
        T.T_11 = T.T_11 - 273.15; %degrees C
        
        %Interpolate 2011 mean annual T along Crane centerline (RACMO Months 385:396)
        T.cl=zeros(length(RACMOy),length(385:396));
        for i=1:length(RACMOy)
            for j=385:396
                T.cl(i,j) = T.T(RACMOy(i),RACMOx(i),j);
            end    
        end  
        
       %Take the average at each point along centerline
        T.cl(:,1:384)=[]; %Started adding points in column 385
        T.cl(:,1) = nanmean(T.cl,2);
        T.cl(:,2:end) = []; %deg K
        T.cl = T.cl - 273.15; %degrees C
        
    %Loop through each point along centerline
    T.interp = ones(length(h_cl_11),1);
    for i=1:length(x)
    
        %Grab points within a certain distance from centerline, interpolate
        maxDist = 7e3; %~3 RACMO grid points
        nearbyPts = []; %Hold points within a certain distance
        numPts = 0; %Total number of points found
    
        xi = cl.X(i); yi = cl.Y(i);
        for j=1:length(h.X(:,1))
            for k=1:length(h.X(1,:))
                dist = sqrt((xi-h.X(j,k))^2 + (yi-h.Y(j,k))^2);
                if dist<=maxDist
                    numPts = numPts+1;
                    nearbyPts(numPts,1:2) = ([h.h(j,k) T.T_11(j,k)]); 
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
        T.interp(i) = m*h_cl_11(i)+b; 
    end 
        
    %Dry adiabatic lapse rate    
    lr = 9.8e-3; %lapse rate (degrees C m^-1)

    %Calculate new air temperature using lapse rate and RACMO reference height
    T.cl_m = movmean(T.cl_11,10,'omitnan');
    h_cl_m = movmean(h_cl_11,10,'omitnan');
    T.cl_11 = (h_cl_RACMO-h_cl_m)*lr+T.cl_m; 

    %Plot height profiles
    figure; hold on; 
        set(gcf,'Units','centimeters','Position',[2 10 22 20]);
        plot(x,h_cl_RACMO,'-b','LineWidth',3,'DisplayName','RACMO height'); 
        plot(x,h_cl_11,'-m','LineWidth',3,'DisplayName','2011 Crane Surface');
        grid on; xlabel('Distance Along Centerline (m)'); ylabel('Elevation (m)');
        set(gca,'FontSize',14,'FontName','Arial'); legend('Location','south'); 
        title('Height Profiles');
        hold off;

    %Plot temperature profile at each surface
    figure; hold on; 
        set(gcf,'Units','centimeters','Position',[26 10 22 20]);
        plot(x,T.interp,'-b','LineWidth',3,'DisplayName','RACMO height'); 
        plot(x,T.cl_11,'-m','LineWidth',3,'DisplayName','2011 Crane Surface');
        grid on; xlabel('Distance Along Centerline (m)'); ylabel('Temperature (^oC)');
        set(gca,'FontSize',14,'FontName','Arial'); legend('Location','south');
        title('RACMO Air Temperature');
        hold off;
    
%% Save adjusted variables
cd /Users/raineyaberle/Desktop/Research/RACMO2.3
save('Crane_downscaledSMB_2011.mat','x','h_cl_11','smb.interp','smb.linear');

save('Crane_adjustedAirTemp_2011.mat','T.cl_11');

%% Repeat above steps to load smb for 2009-2019

close all; clear all;

figure_save = 1; % =1 to save figure
smb_save = 1; % =1 to save final downscaled SMB

years = 2009:2019;
col = parula(length(years)+1);

%Load centerline
    homepath='/Users/raineyaberle/Desktop/Research/';
    cd([homepath,'CraneGlacier_modeling']);
    addpath([homepath,'matlabFunctions']);
    addpath([homepath,'CraneGlacier_modeling/data/RACMO2.3']);
    cl.X = load('Crane_centerline.mat','x').x; 
    cl.Y = load('Crane_centerline.mat','y').y;
    
    % Load most advanced terminus position (2019)
    term = dsearchn([cl.X cl.Y],[load('CraneTerminusPosition_2002-2019.mat').term(61).X ...
    load('CraneTerminusPosition_2002-2019.mat').term(61).Y]); 
    % clip centerline at this point
    %cl.X = cl.X(1:term); cl.Y = cl.Y(1:term);

    %Define x as distance along centerline
    x = zeros(1,length(cl.X));
    for i=2:length(cl.X)
        x(i)=sqrt((cl.X(i)-cl.X(i-1))^2+(cl.Y(i)-cl.Y(i-1))^2)+x(i-1);
    end 
    
    %Convert coords to lon lat
    [cl.Lon,cl.Lat] = ps2wgs(cl.X,cl.Y,'StandardParallel',-71,'StandardMeridian',0);
    
    rho_i = 917; %kg/m^3

%Load RACMO variables
    sf.Lat = ncread('RACMO2.3p2_XPEN055_snowfall_monthly_1979_2016.nc','lat'); %degrees north
    sf.Lon = ncread('RACMO2.3p2_XPEN055_snowfall_monthly_1979_2016.nc','lon'); %degrees east
    sf.sf = squeeze(ncread('RACMO2.3p2_XPEN055_snowfall_monthly_1979_2016.nc','snowfall')); %kg m-2 d-1 

    sm.Lat = ncread('RACMO2.3p2_XPEN055_snowmelt_monthly_1979_2016.nc','lat'); %degrees north
    sm.Lon = ncread('RACMO2.3p2_XPEN055_snowmelt_monthly_1979_2016.nc','lon'); %degrees east
    sm.sm = squeeze(ncread('RACMO2.3p2_XPEN055_snowmelt_monthly_1979_2016.nc','snowmelt')); %kg m-2 d-1

    ro.Lat = ncread('RACMO2.3p2_XPEN055_runoff_monthly_1979_2016.nc','lat'); %degrees north
    ro.Lon = ncread('RACMO2.3p2_XPEN055_runoff_monthly_1979_2016.nc','lon'); %degrees east
    ro.ro = squeeze(ncread('RACMO2.3p2_XPEN055_runoff_monthly_1979_2016.nc','runoff')); %kg m-2 d-1
   
    smb.Lat = ncread('SMB_RACMO2.3p2_monthly_XPEN055_197901_201908.nc','lat'); %degrees north
    smb.Lon = ncread('SMB_RACMO2.3p2_monthly_XPEN055_197901_201908.nc','lon'); %degrees east
    smb.smb = squeeze(ncread('SMB_RACMO2.3p2_monthly_XPEN055_197901_201908.nc','smb')); %kg m-2 d-1
    
    h.Lat = ncread('Height_latlon_XPEN055.nc','lat'); %degrees north
    h.Lon = ncread('Height_latlon_XPEN055.nc','lon'); %degrees east
    h.h = ncread('Height_latlon_XPEN055.nc','height'); %m above the surface
    h_cl_RACMO = griddata(h.Lon,h.Lat,h.h,cl.Lon,cl.Lat); 
    
    %Grab RACMO grid
        %Note: RACMO month # = (year - 1979)*12 + (mo. # in year)
        %Ex: Feb/2000 = (2000-1979)*12 + 2 = 254
        for i=1:length(cl.Lon)
            lat_diff = abs(cl.Lat(i)*ones(size(sf.Lat)) - sf.Lat);
            lon_diff = abs(cl.Lat(i,1)*ones(size(sf.Lon)) - sf.Lon);
            diff_map = sqrt(lat_diff.^2+lon_diff.^2);
            RACMO_ref = find(diff_map==min(min(diff_map)));
            [RACMOy(i),RACMOx(i)] = ind2sub(size(nanmean(sf.sf(:,:,1),4)),RACMO_ref);
        end     
        
    %Convert RACMO grid to polar stereographic
    [h.X,h.Y] = wgs2ps(h.Lon,h.Lat,'StandardParallel',-71,'StandardMeridian',0);
    [sf.X,sf.Y] = wgs2ps(sf.Lon,sf.Lat,'StandardParallel',-71,'StandardMeridian',0);
    [sm.X,sm.Y] = wgs2ps(sm.Lon,sm.Lat,'StandardParallel',-71,'StandardMeridian',0);
    [smb.X,smb.Y] = wgs2ps(smb.Lon,smb.Lat,'StandardParallel',-71,'StandardMeridian',0);
    
    %Load 2011 Crane ice surface
     cd([homepath,'CraneGlacier_modeling/inputs-outputs']);
     h_cl_11 = load("Crane_SurfaceObservations_2009-2018.mat").h(3).surface';
    
    maxDist = 200e3;     
    
        %Grab points within a certain distance from centerline, interpolate
        nearbyPts = []; %Hold points within a certain distance
        numPts = 0; %Total number of points found    

for i=1:length(years)

% Calculate Adjusted Snowfall (sf), Snowmelt (sm), 
% and Surface Mass Balance (SMB) along centerline

    mo_start = (years(i)-1979)*12+1;
    mo_end = (years(i)-1979)*12+12;
    if years(i)==2019
        mo_end = (years(i)-1979)*12+8;
    end 
    
    % set up figures
    if i==1
        fig1=figure;
        set(gcf,'units','centimeters','position',[5 0 40 30]);
        subplot(2,2,1); hold on;  % adjusted snowfall
            set(gca,'FontName','Arial','FontSize',14);
            xlabel('Distance Along Centerline (m)'); ylabel('(m yr^{-1})'); 
            title('Adjusted Snowfall');grid on; legend;
        subplot(2,2,2); hold on; % adjusted snowmelt
            set(gca,'FontName','Arial','FontSize',14);
            xlabel('Distance Along Centerline (m)'); ylabel('(m yr^{-1})'); 
            title('Adjusted Snowmelt');grid on; legend;
        subplot(2,2,3); hold on; % adjusted SMB
            set(gca,'FontName','Arial','FontSize',14);
            xlabel('Distance Along Centerline (m)'); ylabel('(m yr^{-1})'); 
            title('Adjusted SMB');grid on; legend;        
        subplot(2,2,4); hold on; % linearized SMB
            set(gca,'FontName','Arial','FontSize',14);
            xlabel('Distance Along Centerline (m)'); ylabel('(m yr^{-1})'); 
            title('Linearized SMB');grid on; legend; 
    end 

%SNOWFALL (sf)
if years(i)<=2016
    %Grab mean sf for every month for full RACMO grid
        sf.Yr = zeros(length(sf.sf(:,1,1)),length((sf.sf(1,:,1))));
        for j = 1:length(sf.sf(:,1,1))
            for k=1:length(sf.sf(1,:,1))
                sf.Yr(j,k) = nanmean(sf.sf(j,k,mo_start:mo_end));
            end     
        end 
        sf.Yr = sf.Yr./rho_i.*365; %m yr^-1

    %Interpolate mean annual sf along Crane centerline
        sf.cl=[];
        for j=1:length(RACMOy)
            for k=mo_start:mo_end
                sf.cl(j,k) = sf.sf(RACMOy(j),RACMOx(j),k);
            end    
        end  

       %Take the average at each point along centerline
        sf.cl(:,1:(mo_start-1))=[]; %Delete empty columns
        sf.cl(:,1) = nanmean(sf.cl,2);
        sf.cl(:,2:end) = []; %kg m-2 day-1
        sf.cl = sf.cl.*365./rho_i; %m yr-1

   %Adjust snowfall (sf) for elevation-dependence

    %Loop through each point along centerline
    sf.interp = ones(length(h_cl_11),1);
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
                    nearbyPts(numPts,1:2) = ([h.h(k,l) sf.Yr(k,l)]); 
                end 
            end 
        end 

        %Calculate a linear trendline for nearby points
        P = polyfit(nearbyPts(:,1),nearbyPts(:,2),1);
        int = -500:2000; 
        yfit = P(1)*int+P(2);

        sf.interp(j) = interp1(int,yfit,h_cl_11(j));

    end        

    figure(1); subplot(2,2,1);
    plot(x,sf.interp,'-','color',col(i,:),'linewidth',2,'displayname',num2str(years(i))); 
    drawnow
end   

%SNOWMELT (sm)
if years(i)<=2016
    %Grab mean sm for every month for full RACMO grid
        sm.Yr = zeros(length(sm.sm(:,1,1)),length(sm.sm(1,:,1)));
        for j = 1:length(sm.sm(:,1,1))
            for k=1:length(sm.sm(1,:,1))
                sm.Yr(j,k,1) = nanmean(sm.sm(j,k,mo_start:mo_end)); % kg m^2 / day
            end     
        end 
        sm.Yr = sm.Yr./rho_i.*365; %m yr^-1
        
    %Interpolate mean annual sm along Crane centerline
        sm.cl=zeros(length(RACMOy),length(mo_start:mo_end));
        for j=1:length(RACMOy)
            for k=mo_start:mo_end
                sm.cl(j,k) = sm.sm(RACMOy(j),RACMOx(j),k);
            end    
        end  
        
       %Take the average at each point along centerline
        sm.cl(:,1:(mo_start-1))=[]; %Delete empty columns
        sm.cl(:,1) = nanmean(sm.cl,2);
        sm.cl(:,2:end) = []; %kg m-2 day-1
        sm.cl = sm.cl.*365./rho_i; %m yr-1
        
    %Ajust snowmelt for elevation-dependence
    
        %Loop through each point along centerline
        sm.interp = ones(length(h_cl_11),1);
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
            sm.interp(j) = interp1(int,yfit,h_cl_11(j));
        
        end 
    
    figure(1); subplot(2,2,2);
    plot(x,sm.interp,'-','color',col(i,:),'linewidth',2,'displayname',num2str(years(i)));
    drawnow
end

%SURFACE MASS BALANCE (smb)

    %Grab mean smb for every month in year
        smb.Yr = zeros(length(smb.smb(:,1,1)),length(smb.smb(1,:,1)));
        for j = 1:length(smb.smb(:,1,1))
            for k=1:length(smb.smb(1,:,1))
                smb.Yr(j,k,1) = nanmean(smb.smb(j,k,mo_start:mo_end)); % kg /yr
            end     
        end 
        smb.Yr = smb.Yr./rho_i.*365;  % m yr^-1
            
        %Interpolate mean annual smb along Crane centerline (RACMO Months 385:396)
        smb.cl=zeros(length(RACMOy),length(mo_start:mo_end));
        for j=1:length(RACMOy)
            for k=mo_start:mo_end
                smb.cl(j,k) = smb.smb(RACMOy(j),RACMOx(j),k);
            end    
        end  
        
       %Take the average at each point along centerline
        smb.cl(:,1:(mo_start-1))=[]; %Delete empty columns
        smb.cl(:,1) = nanmean(smb.cl,2);
        smb.cl(:,2:end) = []; %kg m-2 day-1
        smb.cl = smb.cl.*365./rho_i; %m yr-1    

    %Ajust surface mass balance for elevation-dependence

    %Loop through each point along centerline
    smb.interp = ones(length(h_cl_11),1);
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
            nearbyPts(:,2),1),h_cl_11(j));
       
    end 
    
    figure(1); subplot(2,2,3);
    plot(x,smb.interp,'-','color',col(i,:),'linewidth',2,'displayname',num2str(years(i)));
    drawnow 
    
    %Calculate smb linear trendline    
    n = find(~isnan(h_cl_11) & ~isnan(smb.interp));
    smb.P = polyval(polyfit(h_cl_11(n),...
        smb.interp(n),1),h_cl_11(n));
        smb.lin(i,1:2) = polyfit(transpose(x(n)),smb.P,1);
        smb.linear(:,i) = polyval(smb.lin(i,:),x);
        
    figure(1); subplot(2,2,4);
    plot(x,smb.linear(:,i),'-','color',col(i,:),'linewidth',2,'displayname',num2str(years(i)));
    
    % Save in structure
    SMB(i).year = years(i); % Year
    SMB(i).smb_linear = smb.linear(:,i); % linear trendline smb (m/yr)
        sigma_h_cl_11 = 25; % uncertainty in ice surface (m) - two points

    SMB(i).smb_interp = smb.interp; % interpolated smb (m/yr)
    
    % Calculate uncertainty in final SMB vector: 
    % sigma_smb: [10% SMB raw product] &
    %           [0.49 mWE/yr from mean bias for downscaled product vs.
    %           smb observed for eight transects in GrIS (Noel et al., 2016)] 
    %           & [uncertainty in ice surface elevation]
    %           = 0.1.*smb_interp + 0.49*(rho_w/rho_i*sigma_smbWE)
    % mWE/yr -> m ice/yr: rho_w*H_w = rho_i*H_i
    %                   ->  H_i = rho_w/rho_i/H_w
    H_w = 0.49; % mWE/yr mean bias vs. observations (Noel et al., 2016)
    rho_i = 917; rho_w = 1000; % kg/m^3
    for j=1:length(x)
        SMB(i).sigma_smb(j) = SMB(i).smb_interp(j).*sqrt((0.1.*smb.interp(j)/smb.interp(j))^2 ...
        +(rho_w/rho_i*H_w/H_w)^2+(sigma_h_cl_11/h_cl_11(j))^2);
    end 
    
end 
        
if figure_save
    cd([homepath,'figures']);
    saveas(fig1,'Crane_AdjustedRACMOVariables_2011-2019.png','png');
    disp('Figure 1 saved');
end 
if smb_save
    cd([homepath,'CraneGlacier_modeling/inputs-outputs']);
    save('Crane_downscaledSMB_2009-2019.mat','SMB');
    disp('SMB saved');
end 


    
    
    
    
    
    
    
    
    
    

