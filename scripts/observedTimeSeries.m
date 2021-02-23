% Rainey Aberle
% Summer/Fall 2020
% Script to create a timeseries of Crane Glacier geometry, ice speed, and 
% terminus position along glacier centerline, 
% used to tune glacier flowline model

% Order of operations:
%   1. Ice surface (PGC & OIB) and bed elevation profile
%   2. Ice speeds (ITS_LIVE & TSX)
%       - Calculate full speed magnitude using surface and bed slopes
%   3. Terminus positions

clear all; close all;

bed_save = 'n';        % Switch 'y'/'n' to save/not save bed
surface_save = 'n';    % Switch 'y'/'n' to save/not save surface
velocity_save = 'y';   % Switch 'y'/'n' to save/not save velocity
dHdt_save = 'n';       % Switch 'y'/'n' to save/not save dHdt
figures_save = 'y';    % Switch 'y'/'n' to save/not save figures 

%% 1. Glacier ice surface (PGC & OIB)
    %Note: OIB columns
    %Elevation = 5
    %Surface elevation = surface(5) - elevation(7)
    %Thickness = 4
    %Bed = elevation(5) - thickness(4)
    %Convert elevations to geoid
    %No data values = -9999
    
    % Add path to folders holding required functions 
    homepath = '/Users/raineyaberle/Desktop/Research/';
    addpath([homepath,'matlabFunctions']);
    addpath([homepath,'matlabFunctions/hugheylab-nestedSortStruct']);
    
    % Load Crane centerline
    cd([homepath,'CraneGlacier_modeling/inputs-outputs']);
    cl.X = load('Crane_centerline.mat').x; cl.Y = load('Crane_centerline.mat').y;
    
        % Define x as distance along centerline
        x = zeros(1,length(cl.X));
        for i=2:(length(cl.X))
            x(i)=sqrt((cl.X(i)-cl.X(i-1))^2+(cl.Y(i)-cl.Y(i-1))^2)+x(i-1);
        end
        
        % Convert to lon/lat, load geoidheight at each pt along centerline
        [cl.Lon,cl.Lat] = ps2wgs(cl.X,cl.Y);
        h_geoid = geoidheight(cl.Lat,cl.Lon);
        
    % If bed not previously saved, load and interpolate bed observations
    if strcmp(bed_save,'y')
        % Grab Tate's observed bed & surface, convert coordinates, 
        % interpolate to centerline

        % Load Tate's OIB file
        cd([homepath,'Crane_modelingOLD/BedandSurface']);
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
        cd([homepath,'CraneGlacier_modeling/inputs-outputs']);
        save('Crane_ObservedBed_Tate.mat','hb');
        disp('bed variable saved');
        
    else
        
        cd([homepath,'CraneGlacier_modeling/inputs-outputs']);
        hb = load('Crane_ObservedBed_Tate.mat').hb.hb0; 
        
    end 
    
    % If surface not previously saved, load and interpolate surface
    % observations
    if strcmp(surface_save,'y')
        % Grab OIB observation dates
        cd([homepath,'Crane_modelingOLD/BedandSurface/OIB/ConvertedOIB']);
        files = dir('Converted*');
        for i=1:length(files)
            OIBDates(i,:) = [num2str(files(i).name(18:21)),'-',num2str(files(i).name(22:23)),'-',num2str(files(i).name(24:25))]; % OIB temporal coverage
        end     

        % Load PGC DEM surface data
        cd([homepath,'Crane_modelingOLD/BedandSurface']);
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

        cd([homepath,'Crane_modelingOLD/BedandSurface/OIB/ConvertedOIB']);
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
        cd([homepath,'CraneGlacier_modeling/inputs-outputs']);
        save('CraneSurfaceObservations_2009-2018.mat','h');
        disp('surface variable saved');
    
    else 
        cd([homepath,'CraneGlacier_modeling/inputs-outputs']);
        h = load('Crane_SurfaceObservations_2009-2018.mat').h;
    end 
    
    % Load bathymetry observations
    %    hb.bathym = load('Crane_BathymetryData.mat').cl_trough;
    %    hb.bathym = -1.*hb.bathym; % Data reported in depth beneath the surface (convert to negative values for plotting)
    
    % Load terminus positions
    cd([homepath,'CraneGlacier_modeling/inputs-outputs']);
    termX = load('LarsenB_centerline.mat').centerline.termx;
    termY = load('LarsenB_centerline.mat').centerline.termy;
    termx = x(dsearchn([cl.X cl.Y],[termX' termY']));
    termdate = load('LarsenB_centerline.mat').centerline.termdate;
        termx = feval(fit(termdate',termx','poly2'),termdate');
    % reformat surface dates to match terminus dates
    hdates = zeros(1,length(h));
    % set up figure
    figure; plot(termdate,termx,'.-b','displayname','term','markersize',15);
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
        %plot(x(1:length(hb.bathym))./10^3,hb.bathym,'--c','linewidth',2,'displayname','Bathymetry (2006)');
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
    
    % If bed hasn't been saved, add the bed with bathymetry observations as
    % final elevations
    if strcmp(bed_save,'y')
        cd([homepath,'Crane_modelingOLD/BedandSurface']);
        save('Crane_ObservedBed.mat','hb_bathymobs','-append');
    end 
    
% Save figure as image if desired
    if strcmp(figures_save,'y')
        cd([homepath,'figures']);
        saveas(gcf,'Crane_ObservedTimeSeries.png');
        disp('Figure 1 saved');
    end
    
% Save dHdt as .mat file
    if strcmp(dHdt_save,'y')
        cd([homepath,'CraneGlacier_modeling/inputs-outputs']);
        save('dHdt_2009-2018.mat','dHdt'); % save dHdt structure
        disp('dHdt variable saved'); 
    end
        
% Ice speeds (TSX, & ITS_LIVE)
        
    % Set up subplot
    figure(2); 
    subplot(2,1,2);
        hold on;
        set(gca,'FontSize',12,'FontName','Arial','linewidth',2,'fontweight','bold');
        grid on; xlabel('Distance Along Centerline (km)'); ylabel('Speed (m/yr)');
        title('b) Observed Glacier Speed'); 
        legend('Location','northwest');
    
    % Load and save velocities if not saved
    if strcmp(velocity_save,'y')
        
        cd([homepath,'CraneGlacier_modeling/data/velocities']);
        
        % load ANT speeds
        ANTfiles = dir('ANT*');
        % Loop through all files, interpolate along centerline
        for i=1:length(ANTfiles)

            [X,Y] = meshgrid(ncread(ANTfiles(i).name,'x'),ncread(ANTfiles(i).name,'y'));
            u = transpose(ncread(ANTfiles(i).name,'v')); % m/y
            u_err = transpose(ncread(ANTfiles(i).name,'v_err'));         
            u(u==-32767) = NaN; %Replace no data values with NaN
        
            U(i).date = str2double(ANTfiles(i).name(11:14));                % observation date
            U(i).speed = interp2(X,Y,u,cl.X,cl.Y)./3.1536e7;                % interpolated speed (m/s)
            U(i).speed_err = interp2(X,Y,u_err,cl.X,cl.Y)./3.1536e7;        % interpolated speed error (m/s)
            U(i).units = "m/s";                                             % speed units
            U(i).source = "ITS-LIVE";                                       % speed data source
            U(i).numPts = length(U(i).speed(~isnan(U(i).speed)));           % number of data points 

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
        cd([homepath,'CraneGlacier_modeling/inputs-outputs']);
        save('Crane_CenterlineSpeeds_2007-2017.mat','U');
        disp('velocity variable saved');
        
    else
        cd([homepath,'CraneGlacier_modeling/inputs-outputs/']);
        U = load('Crane_CenterlineSpeeds_2007-2017.mat').U; 
    end     
    
    % Plot velocities
    col = parula(length(n)+2); % Color scheme for plotting u profiles
    n = [2 4 6 8 9 14 15:19];
    figure(2); subplot(2,1,2);
       for i=1:length(n)
           plot(x./10^3,U(n(i)).speed.*3.1536e7,'-','color',col(i,:),'linewidth',2,...
               'displayname',num2str(round(U(n(i)).date)));
       end 
                      
    if strcmp(figures_save,'y')
        cd([homepath,'figures']);
        saveas(gcf,'CraneCenterlineSpeeds_2007-2017.png');
        disp('Figure 3 saved');
    end 
    
    hold off;
    
%%
% Solve for the magnitude of the velocity using the surface slope and 
% velocity map (which only prescribes the horizontal velocity component)

    figure(4); clf
    col = parula(length(h));

% Loop through all surfaces
for i=1:length(h)
        
    % Create a staggered grid to calculate slope of the ice bed and surface
        hm = (h(i).surface(2:end)+h(i).surface(1:end-1))./2; % 2011 OIB surface (m)
        hm = [hm 0]; hm = movmean(hm,10);
        hbm = (hb.hb(2:end)+hb.hb(1:end-1))./2; % bed elevation (m)
        hbm = [hbm; 0]; hbm=movmean(hbm,10);
        xm = (x(2:end)+x(1:end-1))./2;
        xm = [xm,0];
   % hm = movmean(h(i).surface,20);

     % Match surface year to velocity year
     year = h(i).date(1:4);
        for j=1:length(u)
            if strcmp(u(j).year,year)
                ux = u(j).speed;
            end 
        end 
     
        % Pre-allocate vectors
        theta_h = zeros(1,length(x)-1);
        theta_hb = zeros(1,length(x)-1);
        theta = zeros(1,length(x)-1);
        
        % Linearly interpolate speed where no data exist
        %ux = fit(transpose(x(find(~isnan(ux)))),ux(find(~isnan(ux))),'linearinterp');
        %    ux = feval(ux,x);
        uy = zeros(1,length(x)-1);
        U = zeros(1,length(x)-1); 

    % Loop through all points on the centerline
    for j=2:length(cl.X)-1

        % Calculate surface angle from the horizontal
        theta_h(j) = atan((hm(j+1)-hm(j-1))/(xm(j)-xm(j-1)));
        theta_hb(j) = atan((hbm(j+1)-hbm(j-1))/(xm(j)-xm(j-1)));
        theta(j)=nanmean([theta_h(j) theta_hb(j)]);
        
        % Calculate vertical component of velocity
        uy(j) = ux(j)*tan(theta(j)); 

        % Calculate full velocity
        U(j) = sqrt(ux(j)^2 + uy(j)^2); 

    end 

    % Plot resulting velocity profile
    figure(3); hold on; 
        plot(x(2:end-1)./10^3,U(2:end).*3.1536e7,'-','color',col(i,:),...
            'linewidth',2,'displayname',year);
        set(gca,'fontname','arial','fontsize',14); grid on;
        xlabel('distance along centerline (km)'); ylabel('speed (m/yr)');
        
end 

hold off;

%% 2. Terminus Position

    % Load LarsenB_centerline variable (Dryak & Enderlin)
    cd([homepath,'CraneGlacier_modeling/inputs-outputs']);
    centerline = load('LarsenB_centerline.mat').centerline;
    
    % Set up figures
    fig_termMap = figure; % Terminus positions along centerline
        hold on; plot(cl.X,cl.Y,'-m','linewidth',2,'displayname','Centerline');
        title('Crane Terminus Position Along Centerline');
        set(gca,'FontSize',14,'FontName','Arial','linewidth',2);
        set(gcf,'Units','centimeters','Position',[2 5 25 20]);
        xlabel('Easting (m)'); ylabel('Northing (m)');
        xlim([-2.41e6 -2.405e6]); ylim([1.262e6 1.273e6]);
        legend('Location','northwest');
        grid on; %axis equal
        
    fig_termTime = figure; % Terminus positions over time
        hold on; 
        title('Crane Terminus Position Over Time');
        set(gca,'FontSize',14,'FontName','Arial','LineWidth',2);
        set(gcf,'Units','centimeters','Position',[28 5 25 20]);        
        xlabel('Year'); ylabel('Distance Along Centerline (m)');
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
            figure(fig_termMap); hold on; plot(termx,termy,'.','color',col(i,:),'markersize',30,'displayname',num2str(termDate));
            figure(fig_termTime); hold on; plot(termDate,x(xn),'.','color',col(i,:),'markersize',30);
        else
            figure(fig_termMap); hold on; plot(termx,termy,'.','color',col(i,:),'markersize',30,'handlevisibility','off');
            figure(fig_termTime); hold on; plot(termDate,x(xn),'.','color',col(i,:),'markersize',30);
        end
        
        x_term(i,:) = ([termDate x(xn)]);
        
        term(i).decidate = termDate;
        term(i).X = termx;
        term(i).Y = termy;
        term(i).x = x(xn);
        
    end 
 
    % Save term as .mat variable
    save('CraneTerminusPosition_2002-2019.mat','term');
    
    % Save figure as image
    cd([homepath,'figures']);
    saveas(fig_termMap,'CraneTerminusPosition_2002-2019.png');
    saveas(fig_termTime,'CraneTerminusAlongCL_2002-2019.png');
    disp('Figures 4 & 5 saved');
 
%% 3. Plot observation timeseries, saving figures
% Save an image of the figure for every year to make a time lapse

clear all; close all;

    % Add path to folders holding required functions 
    homepath = '/Users/raineyaberle/Desktop/Research/';
    addpath([homepath,'matlabFunctions']);
    addpath([homepath,'matlabFunctions/hugheylab-nestedSortStruct']);
    
    % Load Crane centerline
    cd([homepath,'CraneGlacier_modeling/inputs-outputs']);
    cl.X = load('Crane_centerline.mat').x; cl.Y = load('Crane_centerline.mat').y;    
        % Define x as distance along centerline
        x = []; x(1)=0;
        for i=2:(length(cl.X))
            x(i)=sqrt((cl.X(i)-cl.X(i-1))^2+(cl.Y(i)-cl.Y(i-1))^2)+x(i-1);
        end
        
    % Load bed elevations
    hb.hb0 = load('Crane_OIBPicks_2018').CraneOIBPicks_2018_full.hb_geoid;
    hb.X =  load('Crane_OIBPicks_2018').CraneOIBPicks_2018_full.Easting;
    hb.Y =  load('Crane_OIBPicks_2018').CraneOIBPicks_2018_full.Northing;
    % Interpolate to centerline
    hb.hb = interp1(hb.Y,hb.hb0,cl.Y);
    
    % Load surface 
    h = load('Crane_SurfaceObservations_2009-2018.mat').h;
    
    % Load terminus positions
    cd([homepath,'CraneGlacier_modeling/inputs-outputs']);
    termX = load('LarsenB_centerline.mat').centerline.termx; 
    termY = load('LarsenB_centerline.mat').centerline.termy;
    termx = x(dsearchn([cl.X cl.Y],[termX' termY']));
    termdate = load('LarsenB_centerline.mat').centerline.termdate;
    % reformat surface dates to match terminus dates
    hdates = zeros(1,length(h));
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
    end
    clear DV*  
    
    % Load velocities
    u = load('Crane_CenterlineSpeeds_2007-2017.mat').U; 

    % Plot surfaces, velocities, & bed
    n = [0 0 1 2 3 5 10 13 15 28 35 36]; % indices of surfaces
    m = [2 4 6 8 15 17 19 22 24 26 28 0]; % indices for velocities     
    col = parula(length(n)+3); % color scale for plotting
    
    figure(5); clf
    set(gcf,'Units','centimeters','position',[5 5 20 25]); 
    subplot(2,1,1);
        ylabel('Elevation (m)'); xlabel('Distance Along Centerline (km)');
        set(gca,'fontname','Arial','fontsize',12,'linewidth',2,'fontweight','bold');
        xlim([0 75]); ylim([-1500 1000]); title('a) Observed Glacier Geometry');
        hold on; grid on; legend('Location','northeast');
    subplot(2,1,2);
        ylabel('Speed (km/yr)'); xlabel('Distance Along Centerline (km)');
        set(gca,'fontname','Arial','fontsize',12,'linewidth',2,'fontweight','bold');
        xlim([0 75]); ylim([0 3.5]); title('b) Observed Glacier Speed'); 
        hold on; grid on; legend('Location','northeast');
    
    % cd to folder for image saves
    cd('/Users/raineyaberle/Desktop/Research/write-ups,etc/AGU2020/ObsVid/');
    
    for i=1:length(n)
        if i==1
            figure(5);
            subplot(2,1,1);
                % plot bed
                plot(x./10^3,hb.hb,'-k','linewidth',2,'displayname','Bed (2018)');
        end % plot bed
        if n(i)==0
            figure(5);
            subplot(2,1,2);
                % plot speed
                plot(x./10^3,u(m(i)).speed.*3.1536e7./10^3,'-','color',col(i,:),'linewidth',2,...
                'displayname',num2str(fix(u(m(i)).year)));
        elseif m(i)==0 % if there are no speed observations
            figure(5);
            subplot(2,1,1); 
            % plot surface
                hi = h(n(i)).surface(1:h(n(i)).Iterm); 
                % plot
                plot(x(1:h(n(i)).Iterm)./10^3,hi,'-','Color',col(i,:),...
                    'linewidth',2,'markersize',10,'displayname',h(n(i)).date(1:4));
            % plot calving front
                % if h=NaN at calving front
                if isnan(h(n(i)).surface(h(n(i)).Iterm))
                    plot([x(h(n(i)).Iterm)./10^3 x(h(n(i)).Iterm)./10^3],...
                    [hb.hb(h(n(i)).Iterm)+10 0],...
                    'linewidth',2,'color',col(i,:),'handlevisibility','off');                    
                else 
                     plot([x(h(n(i)).Iterm)./10^3 x(h(n(i)).Iterm)./10^3],...
                    [hb.hb(h(n(i)).Iterm)+10 h(n(i)).surface(h(n(i)).Iterm)],...
                    'linewidth',2,'color',col(i,:),'handlevisibility','off');
                end 
        else
            figure(5);
            subplot(2,1,1);
                % plot
                hi = h(n(i)).surface(1:h(n(i)).Iterm); 
                plot(x(1:h(n(i)).Iterm)./10^3,hi,'-','Color',col(i,:),...
                    'linewidth',2,'markersize',10,'displayname',h(n(i)).date(1:4));
            % plot calving front 
                % if h=NaN at calving front
                if isnan(h(n(i)).surface(h(n(i)).Iterm))
                    plot([x(h(n(i)).Iterm)./10^3 x(h(n(i)).Iterm)./10^3],...
                    [hb.hb(h(n(i)).Iterm)+10 0],...
                    'linewidth',2,'color',col(i,:),'handlevisibility','off');                    
                else 
                     plot([x(h(n(i)).Iterm)./10^3 x(h(n(i)).Iterm)./10^3],...
                    [hb.hb(h(n(i)).Iterm)+10 h(n(i)).surface(h(n(i)).Iterm)],...
                    'linewidth',2,'color',col(i,:),'handlevisibility','off');
                end 
            subplot(2,1,2);
            % plot speed
                plot(x./10^3,u(m(i)).speed.*3.1536e7./10^3,'-','color',col(i,:),'linewidth',2,...
                'displayname',num2str(fix(u(m(i)).year)));
        end 
        
        % save figure on each iteration to make time lapse
        figname = [num2str((i+2007-1)),'.png'];
        saveas(gcf,figname,'png');
        disp([figname,' saved']);
    end
    
%% 4. Plot study area

homepath='/Users/raineyaberle/Desktop/Research/';
addpath([homepath,'matlabFunctions']);
addpath([homepath,'matlabFunctions/line2arrow-kakearney-pkg-8aead6f/axescoord2figurecoord']);
addpath([homepath,'matlabFunctions/line2arrow-kakearney-pkg-8aead6f/parsepv']);
addpath([homepath,'matlabFunctions/line2arrow-kakearney-pkg-8aead6f/line2arrow']);

save_figure=1; % =1 to save figure 

% load Landsat image
    cd([homepath,'Imagery/LC08_L1GT_218106_20191013_20191018_01_T2']);
    landsat = dir('*B8.TIF');
    [LS,R] = readgeoraster(landsat.name); [ny,nx] = size(LS);
 
    % polar stereographic coordinates of image boundaries
    x = linspace(min(R.XWorldLimits),max(R.XWorldLimits),nx); 
    y = linspace(min(R.YWorldLimits),max(R.YWorldLimits),ny);
    
    % display image
    figure(10); clf; hold on; 
    set(gcf,'position',[400 200 600 500]);
    imagesc(x./10^3,y./10^3,flipud(LS)); colormap("gray");
    axis equal;
    % Specify the region of interest (ROI)
    xlim([-2.4525e3 -2.3725e3]); ylim([1.215e3 1.285e3]);
    set(gca,'fontsize',14,'fontname','Arial','linewidth',2);
    grid on; xlabel('Easting (km)'); ylabel('Northing (km)');
    
    % plot centerline
    cd([homepath,'CraneGlacier_modeling/inputs-outputs']);
    cl.x = load('Crane_centerline.mat').x; cl.y = load('Crane_centerline.mat').y;
    figure(10); 
    % plot line with arrows to show flow direction
    h1 = line(cl.x(1:40)./10^3,cl.y(1:40)./10^3,'linewidth',3,'color',[0.3010, 0.7450, 0.9330]);
    h2 = line(cl.x(40:110)./10^3,cl.y(40:110)./10^3,'linewidth',3,'color',[0.3010, 0.7450, 0.9330]);
    h3 = line(cl.x(110:end)./10^3,cl.y(110:end)./10^3,'linewidth',3,'color',[0.3010, 0.7450, 0.9330]);
    line2arrow(h1,'color',[0.3010, 0.7450, 0.9330],'headwidth',20,'headlength',15);
    line2arrow(h2,'color',[0.3010, 0.7450, 0.9330],'headwidth',20,'headlength',15);    
    line2arrow(h3,'color',[0.3010, 0.7450, 0.9330],'headwidth',20,'headlength',15);
    
    % insert text labels
    % Crane
    text(-2.425e3,1.25e3,'Crane','color','w','fontsize',18,'fontweight','bold');  
    % Former Larsen B Ice Shelf
    txt = sprintf('Former Larsen B \n      Ice Shelf');
    text(-2.4e3,1.281e3,txt,'color','w','fontsize',14,'fontweight','bold');
    % Jorum
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Jorum', ...
        'HeadStyle','none','LineStyle', 'none', 'TextRotation',45,'Position',[.48 .85 0 0],...
        'FontSize',18,'FontName','Arial','fontweight','bold','TextColor','w');    
    % Flask
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Flask', ...
        'HeadStyle','none','LineStyle', 'none', 'TextRotation',45,'Position',[.875 .2 0 0],...
        'FontSize',18,'FontName','Arial','fontweight','bold','TextColor','w');    
    % Mapple
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Mapple', ...
        'HeadStyle','none','LineStyle', 'none', 'TextRotation',66,'Position',[.63 .6 0 0],...
        'FontSize',13,'FontName','Arial','fontweight','bold','TextColor','w');    
    % Melville
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Melville', ...
        'HeadStyle','none','LineStyle', 'none', 'TextRotation',70,'Position',[.69 .6 0 0],...
        'FontSize',13,'FontName','Arial','fontweight','bold','TextColor','w');           
    % Pequod
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Pequod', ...
        'HeadStyle','none','LineStyle', 'none', 'TextRotation',66,'Position',[.78 .65 0 0],...
        'FontSize',13,'FontName','Arial','fontweight','bold','TextColor','w');               
    % Starbuck
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Starbuck', ...
        'HeadStyle','none','LineStyle', 'none', 'TextRotation',66,'Position',[.87 .59 0 0],...
        'FontSize',13,'FontName','Arial','fontweight','bold','TextColor','w');                   
    % Stubb
    
    % plot LIMA inset in figure
    ax=axes('pos',[0.15 0.13 0.25 0.25]);
    set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[]);
    cd([homepath,'Imagery']);
    LIMA = imread('LIMA.jpg');
    imshow(LIMA);
    hold on; plot(1100,4270,'o','color','y','markersize',10,'linewidth',3);
    
    if save_figure
        % save image
        cd([homepath,'figures']);
        saveas(gcf,'StudyArea.png','png');
        disp(['Figure saved in ',pwd]);
    end


    