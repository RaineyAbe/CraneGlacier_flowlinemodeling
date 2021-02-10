%% Make Figures!
% RKA 2020-2021
% Script to plot and save figures for thesis

% Figure 1. Map of the study area - THIS VERSION NOT USED
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

% load sensitivity test no change variables
cd([homepath,'inputs-outputs']);
load('Crane_sensitivityTests_noChange.mat');

cd([homepath,'scripts/100yrScenario/results']);

save_figure = 1; % = 1 to save figure

% Set up figures
figure(10); clf % SMR
set(gcf,'Position',[100 300 1300 400]);
subplot(1,3,1); hold on; grid on; % a) geometry
    set(gca,'fontsize',14,'linewidth',2);
    xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m.a.s.l)'); 
    xlim([35 60]); ylim([-1200 600]);
    plot(x1./10^3,hb1,'-k','linewidth',2);
    % Add text label            
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
        (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' a ','edgecolor','k','fontsize',14,'fontweight','bold','linewidth',1.5);  
subplot(1,3,2); hold on; grid on; % b) speed
    set(gca,'fontsize',14,'linewidth',2);
    xlabel('Distance Along Centerline (km)'); ylabel('Speed (m a^{-1})');  
    xlim([0 60]); ylim([0 1500]);
    % Add text label            
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
        (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' b ','edgecolor','k','fontsize',14,'fontweight','bold','linewidth',1.5); 
subplot(1,3,3); hold on; grid on; % c) \Delta H
    set(gca,'fontsize',14,'linewidth',2);
    xlabel('\Delta SMR (m a^{-1})'); ylabel('\DeltaH (m)');
    xlim([-4 2]); ylim([-1.5 1.5]);
    % Add text label            
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
        (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' c ','edgecolor','k','fontsize',14,'fontweight','bold','linewidth',1.5);         
figure(11); clf % SMB
set(gcf,'Position',[200 100 1300 400]);
subplot(1,3,1); hold on; grid on; % c) geometry
    set(gca,'fontsize',14,'linewidth',2);
    xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m.a.s.l)'); 
    xlim([35 60]); ylim([-1200 600]);        
    plot(x1./10^3,hb1,'-k','linewidth',2);
    % Add text label            
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
        (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' a ','edgecolor','k','fontsize',14,'fontweight','bold','linewidth',1.5);         
subplot(1,3,2); hold on; grid on; % d) speed
    set(gca,'fontsize',14,'linewidth',2);
    xlabel('Distance Along Centerline (km)'); ylabel('Speed (m a^{-1})');     
    xlim([0 60]); ylim([0 1500]); grid on;
    % Add text label            
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
        (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' b ','edgecolor','k','fontsize',14,'fontweight','bold','linewidth',1.5); 
subplot(1,3,3); hold on; grid on; % c) \Delta H
    set(gca,'fontsize',14,'linewidth',2);
    xlabel('\Delta SMB (m a^{-1})'); ylabel('\DeltaH (m)');
    xlim([-4 2]); ylim([-10 25]);
    % Add text label            
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
        (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' c ','edgecolor','k','fontsize',14,'fontweight','bold','linewidth',1.5); 
    
% Loop through results and plot
files = dir('*geom.mat');
Ismr=[]; Ismb=[];
for i=1:length(files)
    if contains(files(i).name,'SMB0')
        Ismr=[Ismr; i];
    elseif contains(files(i).name,'SMR0')
        Ismb=[Ismb; i];
    end
end

col=parula(length(Ismr)+1); % color scheme for plotting
% SMR
for i=1:length(Ismr)
    load(files(Ismr(i)).name);
    figure(10); hold on;
    subplot(1,3,1); hold on; 
        % ice surface
        plot(x1(1:cf)./10^3,hf(1:cf),'color',col(i,:),'linewidth',2);
        % calving front
        plot([x1(cf) x1(cf)]./10^3,[hf(cf) hf(cf)-Hf(cf)],'color',col(i,:),'linewidth',2);
        % floating bed
        plot(x1(glf:cf)./10^3,hf(glf:cf)-Hf(glf:cf),'color',col(i,:),'linewidth',2);
    subplot(1,3,2); hold on; 
        % ice surface speed
        plot(x1(1:cf)./10^3,Uf(1:cf).*3.1536e7,'color',col(i,:),'linewidth',2);
    subplot(1,3,3); hold on; 
        % dH
        dH = mean(H1(1:c1)-Hf(1:c1),'omitnan'); 
        n1 = regexp(files(Ismr(i)).name,'R')+1; % indices to start of dSMR value
        n2 = regexp(files(Ismr(i)).name,'_S')-1; % indices to end of dSMR value
        plot(str2double(files(Ismr(i)).name(n1:n2))*3.1536e7,dH,'o',...
            'markersize',10,'linewidth',2,'color',col(i,:));
        % colormap
        colormap(flipud(col)); c=colorbar; caxis([-4 1]);   
        set(get(c,'Title'),'String','(m a^{-1})');
end

% SMB
for i=1:length(Ismb)
    load(files(Ismb(i)).name);
    figure(11); hold on;
    subplot(1,3,1); hold on; 
        % ice surface
        plot(x1(1:cf)./10^3,hf(1:cf),'color',col(i,:),'linewidth',2);
        % calving front
        plot([x1(cf) x1(cf)]./10^3,[hf(cf)-Hf(cf) hf(cf)],'color',col(i,:),'linewidth',2);
        % floating bed
        plot(x1(glf:cf)./10^3,hf(glf:cf)-Hf(glf:cf),'color',col(i,:),'linewidth',2);
    subplot(1,3,2); hold on; 
        % ice surface speed
        plot(x1(1:cf)./10^3,Uf(1:cf).*3.1536e7,'color',col(i,:),'linewidth',2);
    subplot(1,3,3); hold on;
        % dH
        dH = mean(H1(1:c1)-Hf(1:c1),'omitnan'); 
        n1 = regexp(files(Ismb(i)).name,'B')+1; % indices to start of dSMB value
        n2 = regexp(files(Ismb(i)).name,'_g')-1; % indices to end of dSMB value
        plot(str2double(files(Ismb(i)).name(n1:n2))*3.1536e7,dH,'o',...
            'markersize',10,'linewidth',2,'color',col(i,:));
        % colormap
        colormap(flipud(col)); colorbar; caxis([-4 1]); 
        set(get(c,'Title'),'String','(m a^{-1})');            
end
    
% Save figure 
if save_figure
    % Save in figures folder
    cd([homepath,'../figures/']);
    figure(10);
    saveas(gcf,'Crane_sensitivityTests_SMR.png','png');
    figure(11);
    saveas(gcf,'Crane_sensitivityTests_SMB.png','png');    
    disp(['figure saved in: ',pwd]);
    
    % Save in thesis figures folder
    cd([homepath,'../write-ups/Thesis/figures/']);
    figure(10);
    saveas(gcf,'Crane_sensitivityTests_SMR.png','png');
    figure(11);
    saveas(gcf,'Crane_sensitivityTests_SMB.png','png');    
    disp(['figure saved in: ',pwd]);  
    
end


