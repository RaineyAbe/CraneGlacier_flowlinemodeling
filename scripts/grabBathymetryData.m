% Rainey Aberle
% May 2020
% Script to grab NBP0603 bathymetry data near Crane terminus

close all;

% Load centerline
    cd /Users/raineyaberle/Desktop/Research/Crane_modeling
    x_cl = load('Crane_centerline.mat').x; y_cl = load('Crane_centerline.mat').y;
    cd /Users/raineyaberle/Desktop/Research/general_matlabcodes
    cl_lonlat = ps2wgs(x_cl,y_cl,'StandardParallel',-71,'StandardMeridian',0);
    
% 2011 terminus position
    term_11 = [cl_lonlat(138,1) cl_lonlat(138,2)];
    
% Load and plot bathumetry data

    % Crane Entrance 
    figure; hold on; set(gcf,'Units','centimeters','Position',[0 10 18 15]);
    set(gca,'FontName','Calibri','FontSize',16);
    
    cd /Users/raineyaberle/Desktop/Research/BedandSurface/MGDS_Download/NBP0603
    [Lon_entr,Lat_entr,Z_entr] = grdread2('CraneEntrance.grd');
    surf(Lon_entr,Lat_entr,Z_entr,'edgecolor','none'); view(2);
    plot3(cl_lonlat(:,1),cl_lonlat(:,2),ones(length(x_cl),1).*2000,'-.m','LineWidth',2);
    plot3(term_11(1),term_11(2),2000,'*m','MarkerSize',15);
    xlabel('Lon'); ylabel('Lat'); title('Crane Entrance'); grid on; 
    colorbar;

    % Crane Trough
    figure; hold on; set(gcf,'Units','centimeters','Position',[18 10 18 15]);
    set(gca,'FontName','Calibri','FontSize',16);
    [Lon_trough,Lat_trough,Z_trough] = grdread2('CraneTrough.grd');
    surf(Lon_trough,Lat_trough,Z_trough,'edgecolor','none'); view(2);
    plot3(cl_lonlat(:,1),cl_lonlat(:,2),ones(length(x_cl),1).*2000,'-.m','LineWidth',2);
    plot3(term_11(1),term_11(2),2000,'*m','MarkerSize',15);
    xlabel('Lon'); ylabel('Lat'); title('Crane Trough'); grid on; 
    colorbar;
    
    % CT and CS
    figure; hold on; set(gcf,'Units','centimeters','Position',[36 10 18 15]);
    set(gca,'FontName','Calibri','FontSize',16);
    [Lon_ct,Lat_ct,Z_ct] = grdread2('CTandCS.grd');
    surf(Lon_ct,Lat_ct,Z_ct,'edgecolor','none'); view(2);
    plot3(cl_lonlat(:,1),cl_lonlat(:,2),ones(length(x_cl),1).*2000,'-.m','LineWidth',2);
    plot3(term_11(1),term_11(2),2000,'*m','MarkerSize',15);
    xlabel('Lon'); ylabel('Lat'); title('CT and CS'); grid on; 
    colorbar;

%%
%Interpolate depths that intersect centerline

% Define x as distance along centerline
    x=[]; x(1)=0;
    for i=2:length(x_cl)
        x(i) = sqrt((x_cl(i)-x_cl(i-1)).^2+(y_cl(i)-y_cl(i-1)).^2)+x(i-1);
    end 

% Find points along centerline with data coverage, interpolate data
n1 = find(cl_lonlat(:,1)>=Lon_trough,1);
    cl_trough = interp2(Lon_trough,Lat_trough,Z_trough,cl_lonlat(n1:end,1),cl_lonlat(n1:end,2));
    cl_trough = [ones(length(x)-length(cl_trough),1); cl_trough];
    cl_trough(cl_trough==1) = NaN;
n2 = find(cl_lonlat(:,1)>=Lon_entr,1);
    cl_entr = interp2(Lon_entr,Lat_entr,Z_entr,cl_lonlat(n2:end,1),cl_lonlat(n2:end,2));
    cl_entr = [ones(length(x)-length(cl_entr),1); cl_entr];
    cl_entr(cl_entr==1) = NaN;
n3 = find(cl_lonlat(:,1)>=Lon_ct,1);
    cl_ct = interp2(Lon_ct,Lat_ct,Z_ct,cl_lonlat(n3:end,1),cl_lonlat(n3:end,2));
    cl_ct = [ones(length(x)-length(cl_ct),1); cl_ct];
    cl_ct(cl_ct==1) = NaN;

% Plot
figure; hold on; 
plot(x,cl_trough,'-.b','linewidth',2,'DisplayName','Trough');
plot(x,cl_entr,'-.m','linewidth',2,'displayname','Entrance');
plot(x,cl_ct,'-.r','linewidth',2,'displayname','CT');
xlabel('Distance Along Centerline (m)'); ylabel('Depth (m)');
set(gca,'FontSize',16,'FontName','Calibri','YDir','Reverse'); grid on; 
legend('Location','north');

% Save figure as image
cd /Users/raineyaberle/Desktop/Research/figures
saveas(gcf,'Crane_BathymetryObservations.png');
disp('Figure saved as .png');

% Save variables in .mat file
cd .. ; cd BedandSurface
save('Crane_BathymetryData.mat','x','cl_trough','cl_entr','cl_ct');
hold off;
disp('Variables saved as .mat file');
    
    
