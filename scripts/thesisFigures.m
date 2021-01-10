%% Make Figures!
% RKA 2020-2021
% Script to plot and save figures for thesis

% Figure 1. Map of the study area
% Figure 2. Observed Conditions Time Series
% Figure 3. Sensitivity Tests

close all; clear all;
homepath = '/Users/raineyaberle/Desktop/Research/CraneGlacier_modeling/';

% Load Crane centerline
    cd([homepath,'inputs-outputs']);
    cl.X = load('Crane_centerline.mat').x; cl.Y = load('Crane_centerline.mat').y;
    
        % Define x as distance along centerline
        cl.x = zeros(1,length(cl.X));
        for i=2:(length(cl.X))
            cl.x(i)=sqrt((cl.X(i)-cl.X(i-1))^2+(cl.Y(i)-cl.Y(i-1))^2)+cl.x(i-1);
        end
        
%% Figure 1. Map of the study area

% Add paths to functions
    addpath([homepath,'../matlabFunctions']);
    addpath([homepath,'../matlabFunctions/line2arrow-kakearney-pkg-8aead6f/axescoord2figurecoord']);
    addpath([homepath,'../matlabFunctions/line2arrow-kakearney-pkg-8aead6f/parsepv']);
    addpath([homepath,'../matlabFunctions/line2arrow-kakearney-pkg-8aead6f/line2arrow']);

% load Landsat image
    cd([homepath,'../Imagery/LC08_L1GT_218106_20191013_20191018_01_T2']);
    landsat = dir('*B8.TIF');
    [LS.im,LS.R] = readgeoraster(landsat.name); [LS.ny,LS.nx] = size(LS.im);
 
    % polar stereographic coordinates of image boundaries
    LS.x = linspace(min(LS.R.XWorldLimits),max(LS.R.XWorldLimits),LS.nx); 
    LS.y = linspace(min(LS.R.YWorldLimits),max(LS.R.YWorldLimits),LS.ny);
    
    % display image
    figure(1); clf; hold on; 
    set(gcf,'position',[400 200 600 500]);
    imagesc(LS.x./10^3,LS.y./10^3,flipud(LS.im)); colormap("gray");
    % Specify the region of interest (ROI)
    xlim([-2.4525e3 -2.3725e3]); ylim([1.215e3 1.285e3]);
    set(gca,'fontsize',14,'fontname','Arial','linewidth',2);
    grid on; xlabel('Easting (km)'); ylabel('Northing (km)');
    
    % plot line with arrows to show flow direction
    h1 = line(cl.X(1:40)./10^3,cl.Y(1:40)./10^3,'linewidth',3,'color',[0.3010, 0.7450, 0.9330]);
    h2 = line(cl.X(40:110)./10^3,cl.Y(40:110)./10^3,'linewidth',3,'color',[0.3010, 0.7450, 0.9330]);
    h3 = line(cl.X(110:end)./10^3,cl.Y(110:end)./10^3,'linewidth',3,'color',[0.3010, 0.7450, 0.9330]);
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
    
    % plot LIMA inset in figure
    ax=axes('pos',[0.15 0.13 0.25 0.25]);
    set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[]);
    cd([homepath,'../Imagery']);
    LIMA = imread('LIMA.jpg');
    imshow(LIMA);
    hold on; plot(1100,4270,'o','color','y','markersize',10,'linewidth',3);

%% Figure 2. Observed Conditions Time Series

   cd([homepath,'inputs-outputs']);

    col = parula(length(2002:2019)); % color scheme for plotting
    
    % 1. Ice surface elevation
    h = load('Crane_SurfaceObservations_2009-2018.mat').h;
    col_h = [3 4 5 5 6 6 6 6 6 7 7 7 8 8 9 9 9 10.*ones(1,11) 11*ones(1,7) 12]+5;

    % 2. Glacier bed (OIB)
    hb = load('Crane_ObservedBed_Tate.mat').hb.hb0; % OIB
    bathym = load('Crane_BathymetryData.mat').cl_trough; % Rebesco et al. (2014)
    
    % 3. Glacier terminus position
    termX = load('LarsenB_centerline.mat').centerline.termx;
    termY = load('LarsenB_centerline.mat').centerline.termy;
    termx = cl.x(dsearchn([cl.X cl.Y],[termX' termY']));
    termdate = load('LarsenB_centerline.mat').centerline.termdate;
    col_term = [1 3 4 6 8 13.*ones(1,6) 14.*ones(1,7) 15.*ones(1,6) 16*ones(1,14) 17*ones(1,22) 18]; 
    
    % 4. Ice surface speed
    U = load('Crane_CenterlineSpeeds_2007-2017.mat').U; 
    col_U = [1 1 2 2 3 3 4 4 5 5 5 5 5 6:11]+5;
    
    % 5. Glacier width
    W = load('Crane_CalculatedWidth.mat').width.W;
    
    % 6. SMB
    SMB = load('Crane_downscaledSMB_2009-2016.mat').SMB;
    col_SMB = [8:15];
    
    % Plot
    figure(2); clf
    set(gcf,'Position',[100 100 800 500]);
    subplot(2,2,1); hold on; % geometry
        set(gca,'linewidth',2,'fontsize',12); grid on;
        ylabel('Elevation (m)'); 
        xlim([0 55]); ylim([-1200 1200]);
        for i=1:length(h)
            if i==16
            else
                plot(cl.x./10^3,h(i).surface,'linewidth',1.5,...
                    'color',col(col_h(i),:),'HandleVisibility','off');
            end
        end
        plot(cl.x./10^3,hb,'-k','linewidth',1.5,'displayname','Bed');
        %plot(cl.x./10^3,-bathym,'--k','linewidth',1.5,'displayname','Bathymetry');
        % Add text label
        ta = text((min(get(gca,'XLim')-0.25*(max(get(gca,'XLim')))-min(get(gca,'XLim')))),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
            ' a ','edgecolor','k','fontsize',13,'fontweight','bold','linewidth',1.5);          
    subplot(2,2,2); hold on; % speed
        set(gca,'linewidth',2,'fontsize',12); grid on;
        ylabel('Speed (m/yr)');
        xlim([0 55]);
        for i=1:length(U)
            plot(cl.x./10^3,U(i).speed.*3.1536e7,'linewidth',1.5,'color',col(col_h(i),:));
        end
        % Add text label
        tb = text((min(get(gca,'XLim')-0.25*(max(get(gca,'XLim')))-min(get(gca,'XLim')))),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
            ' b ','edgecolor','k','fontsize',13,'fontweight','bold','linewidth',1.5);      
    subplot(2,2,3); hold on; % terminus position
        set(gca,'linewidth',2,'fontsize',12); grid on;
        xlabel('Distance Along Centerline (km)'); ylabel('Date');
        xlim([0 55]);
        for i=1:length(termx)
            plot(termx(i)./10^3,termdate(i),'ok','MarkerFaceColor',col(col_term(i),:),'markersize',8);
        end
        % Add text label
        tc = text((min(get(gca,'XLim')-0.25*(max(get(gca,'XLim')))-min(get(gca,'XLim')))),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
            ' c ','edgecolor','k','fontsize',13,'fontweight','bold','linewidth',1.5);     
    subplot(2,2,4); hold on; % SMB
        set(gca,'linewidth',2,'fontsize',12); grid on;
        xlabel('Distance Along Centerline (km)'); ylabel('Mean Annual SMB (m/yr)');
        xlim([0 55]); set(gca,'clim',[2002 2019]); 
        c = colorbar('Limits',[2002 2019],'Ticks',[2002 2010 2019],'Position',[.92 .32 .03 .3410]);
        for i=1:length(SMB)
            plot(cl.x(1:137)./10^3,SMB(i).smb_interp,'color',col(col_SMB(i),:),'linewidth',1.5);
        end
        % Add text label
        td = text((min(get(gca,'XLim')-0.25*(max(get(gca,'XLim')))-min(get(gca,'XLim')))),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
            ' d ','edgecolor','k','fontsize',13,'fontweight','bold','linewidth',1.5);          
    % Save figure
    cd([homepath,'../figures']);
    saveas(gcf,'observedTimeSeries.png','png');
    disp(['observedTimeSeries.png saved in: ',pwd]);
    
    cd([homepath,'../write-ups/Thesis/figures']);
    saveas(gcf,'observedTimeSeries.png','png');
    disp(['observedTimeSeries.png saved in: ',pwd]);
    
%% Figure 3. Sensitivity Tests

    close all;
    cd([homepath,'scripts/100yrScenario']);
    
    save_figure = 1; % = 1 to save figure
    
    % Load saved figures
    im_a = hgload('SMR0_SMB-1e-07.fig');
    im_b = hgload('SMR0_SMB-3e-07.fig');
    im_c = hgload('SMR0_SMB1e-07.fig');
    im_d = hgload('SMR-1e-07_SMB0.fig');
    im_e = hgload('SMR-3e-07_SMB0.fig');
    im_f = hgload('SMR1e-07_SMB0.fig');

    % Create figure
    figure(10); clf; set(gcf,'Position',[0 0 1100 1000]);
    
    % Copy figures onto subplots
    figure(10); 
    %SMB
    SP(1)=subplot(2,3,1); hold on; grid on;
        copyobj(allchild(get(im_a,'CurrentAxes')),SP(1));
        xlim([20 50]); ylim([-1000 600]);
        set(gca,'fontsize',12,'linewidth',2);
        ylabel('Elevation (m)');
        title('-3.2 m/a');
        % Add text label
        ta = text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
            ' a ','edgecolor','k','fontsize',13,'fontweight','bold','linewidth',1.5);    
    SP(2)=subplot(2,3,2); hold on; grid on;
        copyobj(allchild(get(im_b,'CurrentAxes')),SP(2));
        xlim([20 50]); ylim([-1000 600]);
        set(gca,'fontsize',12,'linewidth',2);
        title('-9.5 m/a');
        % Add text label            
        tb = text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
            ' b ','edgecolor','k','fontsize',13,'fontweight','bold','linewidth',1.5);     
    SP(3)=subplot(2,3,3); hold on; grid on;
        copyobj(allchild(get(im_c,'CurrentAxes')),SP(3));
        xlim([20 50]); ylim([-1000 600]);
        set(gca,'fontsize',12,'linewidth',2);
        title('+3.2 m/a');
        % Add text label            
         tc = text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.86+min(get(gca,'XLim')),...
             (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
             ' c ','edgecolor','k','fontsize',13,'fontweight','bold','linewidth',1.5);     
    %SMR
    SP(4)=subplot(2,3,4); hold on; grid on;
        copyobj(allchild(get(im_d,'CurrentAxes')),SP(4));
        xlim([20 50]); ylim([-1000 600]);
        set(gca,'fontsize',12,'linewidth',2);
        xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)'); 
        % Add text label            
        td = text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
            ' d ','edgecolor','k','fontsize',13,'fontweight','bold','linewidth',1.5);                
    SP(5)=subplot(2,3,5); hold on; grid on;
        copyobj(allchild(get(im_e,'CurrentAxes')),SP(5));
        xlim([20 50]); ylim([-1000 600]);
        set(gca,'fontsize',12,'linewidth',2);
        xlabel('Distance Along Centerline (km)');
        % Add text label            
        te = text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
            ' e ','edgecolor','k','fontsize',13,'fontweight','bold','linewidth',1.5);     
    SP(6)=subplot(2,3,6); hold on; grid on;
        copyobj(allchild(get(im_f,'CurrentAxes')),SP(6));
        xlim([20 50]); ylim([-1000 600]);
        set(gca,'fontsize',12,'linewidth',2);
        xlabel('Distance Along Centerline (km)');
        legend('no change','change','Position',[0.9 0.48 0.0836 0.0490],...
            'fontsize',12);
        % Add text label            
        tf = text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.88+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
            ' f ','edgecolor','k','fontsize',13,'fontweight','bold','linewidth',1.5);                            
   % Insert text annotations
   tSMB = annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','SMB', ...
        'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.05 .79 0 0],...
        'FontSize',18,'FontName','Arial','fontweight','bold','TextColor',[0.75 0 0]);
   tSMR = annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','SMR', ...
        'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[.05 .31 0 0],...
        'FontSize',18,'FontName','Arial','fontweight','bold','TextColor',[0.75 0 0]);

% Save figure 
if save_figure
    % Save in figures folder
    cd([homepath,'../figures/']);
    saveas(gcf,'Crane_sensitivityTests.png','png');
    disp(['figure saved in: ',pwd]);
    
    % Save in thesis figures folder
    cd([homepath,'../write-ups/Thesis/figures/']);
    saveas(gcf,'Crane_sensitivityTests.png','png');
    disp(['figure saved in: ',pwd]);  
    
end


