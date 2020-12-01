%Rainey Aberle
%May 2020
%Script to calculate temperature-dependent rate factor A 
%using RACMO long-term annual air temperature

clear all; close all;

homepath = '/Users/raineyaberle/Desktop/Research/';
cd(homepath);

%Calculate mean long-term annual air temp and interpolate
%along Crane centerline (RACMO 2000-2015)

    %Load centerline
    cd([homepath,'CraneGlacier_modeling/inputs-outputs']);
    x_cl = load('Crane_centerline.mat').x; y_cl = load('Crane_centerline.mat').y;
    
    addpath([homepath,'matlabFunctions']);    
    cl_lonlat = ps2wgs(x_cl,y_cl,'StandardParallel',-71,'StandardMeridian',0);
    
    % Define x as distance along centerline
    x = []; x(1)=0; 
    for i=2:length(x_cl)
        x(i) = sqrt((x_cl(i)-x_cl(i-1))^2+(y_cl(i)-y_cl(i-1))^2)+x(i-1);
    end 

    %Load 2011 OIB Crane surface elevation
    h_cl_11 = load('Crane_SurfaceObservations_2009-2018.mat').h(3);
    h_cl_11 = h_cl_11.surface;
    
    %Load RACMO air temperature and height
    cd([homepath,'CraneGlacier_modeling/data/RACMO2.3']);
    T_lat = ncread('RACMO2.3p2_XPEN055_T2m_monthly_1979_2016.nc','lat'); %degrees north
    T_lon = ncread('RACMO2.3p2_XPEN055_T2m_monthly_1979_2016.nc','lon'); %degrees east
    T = squeeze(ncread('RACMO2.3p2_XPEN055_T2m_monthly_1979_2016.nc','t2m')); % degrees K
    h_lat = ncread('Height_latlon_XPEN055.nc','lat'); %degrees north
    h_lon = ncread('Height_latlon_XPEN055.nc','lon'); %degrees east
    h = ncread('Height_latlon_XPEN055.nc','height'); %m above the surface
    h_cl_RACMO = griddata(h_lon,h_lat,h,cl_lonlat(:,1),cl_lonlat(:,2));
    
    % Convert to polar stereographic coordinates
    [h_x,h_y] = wgs2ps(h_lon,h_lat,'StandardParallel',-71,'StandardMeridian',0);    
    
        %Grab RACMO grid
        for i=1:length(cl_lonlat(:,1))
            lat_diff = abs(cl_lonlat(i,2)*ones(size(T_lat)) - T_lat);
            lon_diff = abs(cl_lonlat(i,1)*ones(size(T_lon)) - T_lon);
            diff_map = sqrt(lat_diff.^2+lon_diff.^2);
            RACMO_ref = find(diff_map==min(min(diff_map)));
            [RACMOy(i),RACMOx(i)] = ind2sub(size(squeeze(nanmean(T(:,:,270:450),4))),RACMO_ref);
        end
       
    % Calculate mean annual air temp for full RACMO grid (2000-2015)
    T_m=zeros(length(T(:,1,1)),length(T(1,:,1)));
    for i=1:length(T(:,1,1))  
        for j=1:length(T(1,:,1))
            T_m(i,j) = nanmean(T(i,j,253:444))-273.15;
        end 
    end
    
    figure; hold on; set(gcf,'Units','centimeters','Position',[0 8 18 16]);
    surf(T_lon,T_lat,T_m); view(2);
    plot3(cl_lonlat(:,1),cl_lonlat(:,2),ones(length(cl_lonlat(:,1)))*3000,'-m','LineWidth',4,'DisplayName','Crane Centerline');
    set(gca,'FontName','Calibri','FontSize',12);
    xlabel('Lon'); ylabel('Lat'); title('RACMO 2000-15 Mean Air Temp'); 
    c = colorbar; set(get(c,'title'),'string','^oC');
    ylim([-66 -65]); xlim([-63.5 -61.5]); 
    hold off;
    
    % Interpolate mean 2000-2015 air temperature along Crane centerline
    % (RACMO months 253:444)
    T_cl=[];
        for i=1:length(RACMOy)
            for j=253:444
                T_cl(i,j) = T(RACMOy(i),RACMOx(i),j);
            end    
        end  
        
       %Take the average at each point along centerline
        T_cl(:,1:252)=[]; %Started adding points in column 253
        T_cl(:,1) = nanmean(T_cl,2);
        T_cl(:,2:end) = []; %deg K
        T_cl = T_cl - 273.15; %degrees C
     
    % Adjust temperature for elevation dependence using a dry adiabatic
    % lapse rate
    T_adj = zeros(length(h_cl_11),1);
    for i=1:length(x)
    
        % Grab points within a certain distance from centerline, interpolate
        maxDist = 10e3; % ~3 RACMO grid points
        nearbyPts = []; % Hold points within a certain distance
        numPts = 0; % Total number of points found
    
        xi = x_cl(i); yi = y_cl(i);
        for j=1:length(h(:,1))
            for k=1:length(h(1,:))
                dist = sqrt((xi-h_x(j,k))^2 + (yi-h_y(j,k))^2);
                if dist<=maxDist
                    numPts = numPts+1;
                    nearbyPts(numPts,1:2) = ([h(j,k) T_m(j,k)]); 
                end 
            end 
        end 
        
        % Calculate a linear trendline for nearby points
        P = polyfit(nearbyPts(:,1),nearbyPts(:,2),1);
        int = -10:1800; 
        yfit = P(1)*int+P(2);
    
        % Apply regression slope to current grid cell, solve for eqn
        m = P(1); b = T_cl(i)-m*h_cl_RACMO(i);
        Tfit = [h_cl_RACMO m.*h_cl_RACMO+b];
        T_adj(i) = m*h_cl_11(i)+b; 
        
    end 
        
    %Dry adiabatic lapse rate    
    lr = 9.8e-3; %lapse rate (degrees C m^-1)
    
    %Calculate new air temperature using lapse rate and RACMO reference height
    h_cl_m = movmean(h_cl_11,10,'omitnan');
    T_cl_11 = (h_cl_RACMO-transpose(h_cl_11))*lr+T_cl;        
    
    %Plot
    figure; hold on; set(gcf,'Units','centimeters','Position',[17 8 17 16]);
    plot(x,T_cl_11,'-b','linewidth',2);
    grid on; set(gca,'fontname','Calibri','fontsize',12);
    xlabel('Distance Along Centerline (m)'); ylabel('Temperature (^oC)');
    title('Mean 2000-15 Air Temp Along Crane Centerline');
    hold off;
        
%Calculate rate factor A assuming an Arrhenius relationship
    T_cl_m = movmean(T_cl,10,'omitnan') + 273.15; %K
    Q = 115; %J/mol
    R = 8.314; %J/(K*mol)
    A = (3.5*10^-25).*exp(-Q./(R.*(T_cl_m))); %Pa^-3 s^-1
    
    figure; set(gcf,'Units','centimeters','Position',[34 8 17 16]);
    plot(x,A,'-m','LineWidth',3); grid on;
    set(gca,'FontSize',12,'FontName','Calibri'); 
    xlabel('Distance Along Centerline (m)'); ylabel('A (Pa^-^3 s^-^1)');
    title('Creep Parameter Along Crane Centerline');
    
% Save rate factor
cd([homepath,'CraneGlacier_modeling/inputs-outputs']);
save('Crane_RateFactorA.mat','A');
disp('rate factor variable saved');
