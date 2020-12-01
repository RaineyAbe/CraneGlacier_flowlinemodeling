% Rainey Aberle
% Fall 2019 / Spring 2020
% Script to grab OIB surface, thickness, and bed elevations along Crane Glacier
% and a tributary (TribA) centerlines 
% Run each section in order (lower sections refer to variables previously defined). 
    
    % Note: OIB columns
    % Elevation = 5
    % Surface elevation = surface(5) - elevation(7)
    % Thickness = 4
    % Bed = elevation(5) - bottom(8)
    % Convert elevations to geoid
    % No data values = -9999

close all; clear all;

% Set up figures

    fig1 = figure;  % OIB flight paths
        set(gca,'FontName','Calibri','FontSize',14);
        set(gcf,'Units','centimeters','Position',[15 8 24 21]);
        xlabel('Easting (m)'); ylabel('Northing (m)');
        title('Crane OIB Flight Paths'); legend; 
        grid on; grid minor; axis equal; 
    fig2 = figure;  % Surface 
        set(gca,'FontName','Calibri','FontSize',14);
        set(gcf,'Units','centimeters','Position',[0 10 17 16]);
        xlabel('Distance Along Centerline (m)'); ylabel('Elevation (m)');
        title('Crane OIB Surface'); legend; grid on;    
    fig3 = figure;  % Thickness
        set(gca,'FontName','Calibri','FontSize',14);
        set(gcf,'Units','centimeters','Position',[17 10 17 16]);        
        xlabel('Distance Along Centerline (m)'); ylabel('Elevation (m)');
        title('Crane OIB Thickness'); legend('Location','northwest'); grid on;
    fig4 = figure;  % Bed
        set(gca,'FontName','Calibri','FontSize',14);
        set(gcf,'Units','centimeters','Position',[36 10 17 16]);        
        xlabel('Distance Along Centerline (m)'); ylabel('Elevation (m)');
        title('Crane OIB Bed'); legend; grid on;
        
    fig5 = figure; % Bed observations on map
        set(gca,'FontName','Calibri','FontSize',14);
        set(gcf,'Units','centimeters','Position',[15 8 24 21]);        
        xlabel('Easting (m)'); ylabel('Northing (m)');
        title('Crane OIB Bed Observations'); legend; grid on;        

% Load OIB CSVs and Landsat Image, display flight paths

    % Load CSVs (See ConvertCoordinates.m if coords are wonky),
    % Convert to polar stereographic coordinates, plot onto Landsat image

    % Display Landsat image
        cd /Users/raineyaberle/Desktop/Research/Imagery/LC08_L1GT_218106_20191013_20191018_01_T2
        disp('Loading Landsat image...');
        landsat = dir('*B8.TIF');
        [LS,R] = readgeoraster(landsat.name); [ny,nx] = size(LS);
        % Increase intensity to increase brightness
        LS = LS.*1.3;
        
        % Polar stereographic coordinates of image boundaries
        x = linspace(min(R.XWorldLimits),max(R.XWorldLimits),nx); 
        y = linspace(min(R.YWorldLimits),max(R.YWorldLimits),ny);
        
        figure(1)
        hold on; imagesc(x,y,flipud(LS)); colormap("gray");
        % Specify the region of interest (ROI)
        x1 = -2.44e6;x2 = -2.385e6; y1 = 1.22e6;y2 = 1.28e6;
        xlim([x1 x2]); ylim([y1 y2]); 
        
        figure(5)
        grid on; hold on; imagesc(x,y,flipud(LS)); colormap("gray");
        % Specify the region of interest (ROI)
        x1 = -2.44e6;x2 = -2.385e6; y1 = 1.22e6;y2 = 1.28e6;
        xlim([x1 x2]); ylim([y1 y2]);
        
     % Load centerline
     cd /Users/raineyaberle/Desktop/Research/CraneGlacier_modeling/inputs-outputs
     x_cl = load('Crane_centerline.mat').x; y_cl = load('Crane_centerline.mat').y;
     
        % Define x as distance along centerline
        x=[]; x(1)=0;
        for i=2:length(x_cl)
            x(i) = sqrt((x_cl(i)-x_cl(i-1)).^2+(y_cl(i)-y_cl(i-1)).^2)+x(i-1);
        end
        
        % Convert to lat lon to grab geoid height
        cd ..; cd general_matlabcodes
        [lon_cl,lat_cl] = ps2wgs(x_cl,y_cl);
        h_geoid = geoidheight(lat_cl,lon_cl);

    % Loop through all flight paths
    disp('Plotting flight paths...')
    cd /Users/raineyaberle/Desktop/Research/Crane_modeling/BedandSurface/OIB/ConvertedOIB
    
    years = [2009 2010 2011 2016 2017 2018];
    files = dir('Converted*.csv');
    
    col = parula(length(files)+1); %Color scheme for plotting
    
    % Use data points within a certain distance of centerline
    maxDist = 1200;

    for i=1:length(years)
        
        %Load .csv file (See ConvertCoordinates.m if coordinates are wonky)
        cd /Users/raineyaberle/Desktop/Research/Crane_modeling/BedandSurface/OIB/ConvertedOIB
        ibcsv = files(i);
        ib = readmatrix(ibcsv.name);
        ib_lat = ib(:,1); ib_lon = ib(:,2);    

        %Create new matrix w polar stereographic (ps) coordinates
        cd /Users/raineyaberle/Desktop/Research/general_matlabcodes/

        P = wgs2ps(ib_lon,ib_lat,'StandardParallel',-71,'StandardMeridian',0);
        ps_ib = ib; 
        ps_ib(:,1) = P(:,1); 
        ps_ib(:,2) = P(:,2);
       
        %Plot OIB flight path
        figure(1)
        plot(ps_ib(:,1), ps_ib(:,2),'Color',col(i,:),'LineWidth',3,'DisplayName',num2str(years(i)));
        hold on; drawnow
        set(gca,'FontName','Calibri','FontSize',14);
        xlabel('Easting (m)'); ylabel('Northing (m)');
        title('OIB Flight Paths at Crane Glacier');
        
        % Extract geometry observations 
        % (by finding OIB point nearest to each point along centerline
        ps_ib(ps_ib==-9999)=NaN; %Replace no data values with NaN
        surf = []; thick = []; bed = [];
        maxDist = 1e3;
        for j=1:length(x_cl)
            n = dsearchn(ps_ib(:,1:2),[x_cl(j) y_cl(j)]);
            dist = sqrt((ps_ib(n,1)-x_cl(j))^2+(ps_ib(n,2)-y_cl(j))^2);
            if dist>=maxDist
                surf = ([surf; ps_ib(n,1:2) NaN]);
                thick = ([thick; ps_ib(n,1:2) NaN]);
                bed = ([bed; ps_ib(n,1:2) NaN]);                
            else
                surf = ([surf; ps_ib(n,1:2) ps_ib(n,5)-ps_ib(n,7)-h_geoid(j)]);
                thick = ([thick; ps_ib(n,1:2) ps_ib(n,4)]);
                bed = ([bed; ps_ib(n,1:2) ps_ib(n,5)-ps_ib(n,8)-h_geoid(j)]);
            end 
        end   
        
        % Save ps_ib variables
        Crane_OIB(i).year = years(i);
        Crane_OIB(i).bed = bed(:,3);
        Crane_OIB(i).surf = surf(:,3);
        Crane_OIB(i).thick = thick(:,3);
        Crane_OIB(i).X = x_cl; 
        Crane_OIB(i).Y = y_cl;
        if years(i)==2010
            ps_ib_2010 = ps_ib; 
        elseif years(i)==2011
            ps_ib_2011 = ps_ib;
        end 
        
        % Plot results
        figure(2)
        hold on; plot(x,surf(:,3),'.','Color',col(i,:),'DisplayName',num2str(years(i)),'MarkerSize',15);
        drawnow
        
        figure(3)
        hold on; plot(x,thick(:,3),'.','Color',col(i,:),'DisplayName',num2str(years(i)),'MarkerSize',15);
        drawnow
        
        figure(4)
        hold on; plot(x,bed(:,3),'.','Color',col(i,:),'DisplayName',num2str(years(i)),'MarkerSize',15);
        drawnow 
        
        % Option: Plot real bed observations on Landsat image
        n = find(~isnan(bed(:,3)));
        figure(5)
        hold on; plot(bed(n,1),bed(n,2),'.y','MarkerSize',20,'DisplayName',num2str(years(i)));
        
        disp(num2str(years(i)));        
        
    end 
    
     % Save figure as image
%      cd /Users/raineyaberle/Desktop/Research/figures
%      saveas(fig1,'Crane_OIBFlightPaths.png');
%      saveas(fig2,'Crane_OIBSurface.png');
%      saveas(fig3,'Crane_OIBThickness.png');
%      saveas(fig4,'Crane_OIBBed.png');
%      saveas(fig5,'Crane_OIBRealBedObservations.png');
%      disp('Figures 1-5 saved as .png files');
     
    % Save observations as .mat file
%    cd ..; cd Crane_modeling/BedandSurface
%    save('Crane_OIBObservations_2009-2018.mat','Crane_OIB');
%    disp('Observations saved as .mat file');
     
%%
% Load OIB surface and bed elevations for Crane Tributary (TribA)
    %Note: Only 2010 and 2011 missions pass over

close all;

    figure; imagesc(x,y,flipud(LS)); caxis([0 250]); colormap("gray");
    set(gca, 'YDir','normal','FontSize',16,'FontName','Times New Roman');
    xlabel('Easting (m)'); ylabel('Northing (m)');
    grid on; grid minor; axis equal; hold on;
    
    %Zoom to Crane 
    xlim([x1 x2]); ylim([y1 y2]); 
    
    %Load Tributary A centerline and flux gate, and plot onto image 
    cd /Users/raineyaberle/Desktop/Research/BedandSurface
    load('Crane_TributaryA.mat','xa_cl','ya_cl','xa_w','ya_w');
    plot(xa_cl,ya_cl,'-.m','LineWidth',3,'DisplayName','TribA Centerline'); 
    plot(xa_w,ya_w,'-.g','LineWidth',3,'DisplayName','TribA Flux Gate');
    plot(x_cl,y_cl,'-.r','LineWidth',3,'DisplayName','Crane Centerline');
    legend('Location','east'); hold off;
    
    %Find OIB points within a certain distance of TribA centerline
    maxDist = 300; 
    
        %2010
        TribA_surf_2010 = []; TribA_thick_2010 = []; TribA_bed_2010 = [];
        for i=1:length(ps_ib_2010)
            test=[];
            xi=ps_ib_2010(i,1); yi=ps_ib_2010(i,2);
            for j=1:length(xa_cl)
                xj = xa_cl(j); yj = ya_cl(j);
                dist = sqrt((xi-xj)^2+(yi-yj)^2);
                if dist<maxDist
                    test(j)=1;
                else 
                    test(j)=NaN;
                end 
            end
            if(any(test==1))
                TribA_surf_2010 = ([TribA_surf_2010; ps_ib_2010(i,1:2) ps_ib_2010(i,5)-ps_ib_2010(i,7)]);
                TribA_thick_2010 = ([TribA_thick_2010; ps_ib_2010(i,1:2) ps_ib_2010(i,4)]);
                TribA_bed_2010 = ([TribA_bed_2010; ps_ib_2010(i,1:2) ps_ib_2010(i,5)-ps_ib_2010(i,8)]);            end
        end 
        disp('2010');
        
        %2011
        TribA_surf_2011 = []; TribA_thick_2011 = []; TribA_bed_2011 = [];
        for i=1:length(ps_ib_2011)
            test=[];
            xi=ps_ib_2011(i,1); yi=ps_ib_2011(i,2);
            for j=1:length(xa_cl)
                xj = xa_cl(j); yj = ya_cl(j);
                dist = sqrt((xi-xj)^2+(yi-yj)^2);
                if dist<maxDist
                    test(j)=1;
                else 
                    test(j)=NaN;
                end 
            end
            if(any(test==1))
                TribA_surf_2011 = ([TribA_surf_2011; ps_ib_2011(i,1:2) ps_ib_2011(i,5)-ps_ib_2011(i,7)]);
                TribA_thick_2011 = ([TribA_thick_2011; ps_ib_2011(i,1:2) ps_ib_2011(i,4)]);
                TribA_bed_2011 = ([TribA_bed_2011; ps_ib_2011(i,1:2) ps_ib_2011(i,5)-ps_ib_2011(i,8)]);
            end
        end 
        disp('2011');  
        
        %Plot TribA OIB surface, thickness, and bed elevations
        disp('Plotting h, H, & hb...');
        figure; hold on;
        plot(TribA_surf_2010(:,1),TribA_surf_2010(:,3),'.b','DisplayName','2010','MarkerSize',8);
        plot(TribA_surf_2011(:,1),TribA_surf_2011(:,3),'.c','DisplayName','2011','MarkerSize',8);
        grid on; xlabel('Easting (m)'); ylabel('Elevation (m)'); legend;
        set(gca,'FontName','Calibri','FontSize',14);
        title('TribA OIB Surface Elevations'); hold off;

        figure; hold on;
        plot(TribA_thick_2010(:,1),TribA_thick_2010(:,3),'.b','DisplayName','2010','MarkerSize',8);
        plot(TribA_thick_2011(:,1),TribA_thick_2011(:,3),'.c','DisplayName','2011','MarkerSize',8);
        grid on; xlabel('Easting (m)'); ylabel('Elevation (m)'); legend;
        set(gca,'FontName','Calibri','FontSize',14);
        title('TribA OIB Thickness'); hold off;
        
        figure; hold on;
        plot(TribA_bed_2010(:,1),TribA_bed_2010(:,3),'.b','DisplayName','2010','MarkerSize',8);
        plot(TribA_bed_2011(:,1),TribA_bed_2011(:,3),'.c','DisplayName','2011','MarkerSize',8);
        grid on; xlabel('Easting (m)'); ylabel('Elevation (m)'); legend;
        set(gca,'FontName','Calibri','FontSize',14);
        title('TribA OIB Bed Elevations'); hold off;