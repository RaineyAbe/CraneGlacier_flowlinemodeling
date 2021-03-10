%% Make Figures!
% RKA 2021
% Script to plot and save figures for thesis

% Figure 1. Map of the study area - THIS VERSION NOT USED
% Figure 2. Observed Conditions Time Series
% Figure 3. Sensitivity Tests

close all; clear all;
homepath = '/Users/raineyaberle/Desktop/Research/CraneGlacier_modeling/';

% Load Crane centerline
    cd([homepath,'inputs-outputs']);
    cl.Xi = load('Crane_centerline.mat').x; cl.Yi = load('Crane_centerline.mat').y;
    
    % Define x as distance along centerline
    cl.xi = zeros(1,length(cl.Xi));
    for i=2:(length(cl.Xi))
        cl.xi(i)=sqrt((cl.Xi(i)-cl.Xi(i-1))^2+(cl.Yi(i)-cl.Yi(i-1))^2)+cl.xi(i-1);
    end
    
    % Regrid the centerline to 300m equal spacing between points
    cl.x = 0:300:cl.xi(end);
    cl.X = interp1(cl.xi,cl.Xi,cl.x); cl.Y = interp1(cl.xi,cl.Yi,cl.x);
        
%% Figure 1. Map of the study area

save_figure = 0; % = 1 to save figure

% Add paths to functions
    addpath([homepath,'../matlabFunctions']);
    addpath([homepath,'../matlabFunctions/cmocean_v2.0/cmocean']);
    addpath([homepath,'../matlabFunctions/AntarcticMappingTools_v5.17/AntarcticMappingTools']);
    addpath([homepath,'../matlabFunctions/textborder']);

% Specify the ROI
    xlimits=[-2.4525e6 -2.3725e6]; ylimits = [1.215e6 1.285e6];
    
% load & display Landsat image
    cd([homepath,'../Imagery/LC08_L1GT_218106_20191013_20191018_01_T2']);
    landsat = dir('*B8.TIF');
    [LS.im,LS.R] = readgeoraster(landsat.name); [LS.ny,LS.nx] = size(LS.im);
    % polar stereographic coordinates of image boundaries
    LS.x = linspace(min(LS.R.XWorldLimits),max(LS.R.XWorldLimits),LS.nx); 
    LS.y = linspace(min(LS.R.YWorldLimits),max(LS.R.YWorldLimits),LS.ny);
    % display image
    figure(1); clf; hold on; 
    set(gcf,'position',[452 75 773 622],'color','w');
    set(gca,'YTick',[],'XTick',[]);
    ax1=axes;
    im1 = imagesc(ax1,LS.x,LS.y,flipud(LS.im*1.1)); colormap('gray');
    set(ax1,'YDir','normal','XLim',xlimits','YLim',ylimits,'FontSize',14,...
        'linewidth',2,'fontsize',14);
    xlabel('Easting (m)'); ylabel('Northing (m)');

% load & display REMA as contours
    cd([homepath,'../Imagery']);
    [REMA.im,REMA.R] = readgeoraster('REMA_clipped.tif');
    REMA.im(REMA.im==-9999)=NaN;
    [REMA.ny,REMA.nx] = size(REMA.im);
    % polar stereographic coordinates of image boundaries
    REMA.x = linspace(min(REMA.R.XWorldLimits),max(REMA.R.XWorldLimits),REMA.nx); 
    REMA.y = linspace(min(REMA.R.YWorldLimits),max(REMA.R.YWorldLimits),REMA.ny);
    hold on; contour(ax1,REMA.x,REMA.y,flipud(REMA.im),0:500:2000,'-k','linewidth',2,'ShowText','on');
    
% load & display ITS_LIVE velocity map
    cd([homepath,'data/velocities']);
    v.v = ncread('ANT_G0240_2017.nc','v'); v.v(v.v==-3267)=NaN;
    v.x = ncread('ANT_G0240_2017.nc','x'); v.y = ncread('ANT_G0240_2017.nc','y');
    ax2 = axes('XTick',[],'YTick',[]);
    im2=imagesc(ax2,v.x,v.y,v.v'); colormap(cmocean('haline'));
    im2.AlphaData=0.5; caxis([0 1100]);
    set(ax2,'YDir','normal','Visible','off','fontsize',14); 
    l=legend('Position',[0.15 0.72 0.12 0.07],'Color',[175 175 175]/255,...
        'fontname','arial');
    hold on; set(gca,'XLim',xlimits,'YLim',ylimits); 
    
% plot OIB flights
    cd([homepath,'../figures/1_StudyArea']);
    OIBfiles = dir('IRMCR2*.csv');
    for i=1:length(OIBfiles)
        OIB.file = readmatrix(OIBfiles(i).name);
        [OIB.x,OIB.y] = wgs2ps(OIB.file(:,2),OIB.file(:,1),'StandardParallel',-71,'StandardMeridian',0);
        if i==1
            hold on; plot(ax2,OIB.x,OIB.y,'color','w',...
                'displayname','NASA OIB','linewidth',1);
        else
            hold on; plot(ax2,OIB.x,OIB.y,'color','w',...
                'HandleVisibility','off','linewidth',1);            
        end
    end
    
% plot centerline points
    plot(cl.X,cl.Y,'o','color','m','markerfacecolor','m',...
        'markeredgecolor','k','markersize',5,'displayname','Centerline');
  
% Link axes
    linkaxes([ax1 ax2]);
    
% Add lat lon coordinates to top and right axes
    ax1; graticuleps(-66:0.5:-60,-68:0.5:-60,'color','k','linewidth',1,...
        'HandleVisibility','off');

% Add colorbar
    c = colorbar('Position',[0.2 0.53 0.02 0.15],'fontsize',12,'fontweight',...
        'bold','color','w','fontname','arial');
    set(get(c,'Title'),'String','Speed (m a^{-1})','color','w','fontname',...
        'arial','Position',[5 95 0]);

% insert text labels
    % Crane
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Crane', ...
       'HeadStyle','none','LineStyle', 'none', 'TextRotation',65,'Position',[.53 0.65 0 0],...
       'FontSize',16,'FontName','Arial','fontweight','bold','TextColor','w');
    % Former Larsen B Ice Shelf
    txt = sprintf('Former Larsen \n    B Ice Shelf');
    text(-2.388e6,1.28e6,txt,'color','w','fontsize',16,'fontweight','bold');
    % Jorum
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Jorum', ...
        'HeadStyle','none','LineStyle', 'none', 'TextRotation',35,'Position',[.52 .9 0 0],...
        'FontSize',16,'FontName','Arial','fontweight','bold','TextColor','w');    
    % Flask
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Flask', ...
        'HeadStyle','none','LineStyle', 'none', 'TextRotation',45,'Position',[.9 .28 0 0],...
        'FontSize',16,'FontName','Arial','fontweight','bold','TextColor','w');    
    % Mapple
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Mapple', ...
        'HeadStyle','none','LineStyle', 'none', 'TextRotation',66,'Position',[.625 .6 0 0],...
        'FontSize',14,'FontName','Arial','fontweight','bold','TextColor','w');    
    % Melville
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Melville', ...
        'HeadStyle','none','LineStyle', 'none', 'TextRotation',70,'Position',[.7 .6 0 0],...
        'FontSize',14,'FontName','Arial','fontweight','bold','TextColor','w');           
    % Pequod
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Pequod', ...
        'HeadStyle','none','LineStyle', 'none', 'TextRotation',66,'Position',[.75 .59 0 0],...
        'FontSize',14,'FontName','Arial','fontweight','bold','TextColor','w');               
    % Starbuck
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Starbuck', ...
        'HeadStyle','none','LineStyle', 'none', 'TextRotation',38,'Position',[.765 .445 0 0],...
        'FontSize',14,'FontName','Arial','fontweight','bold','TextColor','w');   
    % Tributary A
    txt = sprintf('A');
        text(-2.413e6,1.232e6,txt,'color','w','fontsize',16,'fontweight','bold');       
    % Tributary B
    txt = sprintf('B');
        text(-2.415e6,1.237e6,txt,'color','w','fontsize',16,'fontweight','bold');    
    % Tributary C
    txt = sprintf('C');
        text(-2.417e6,1.251e6,txt,'color','w','fontsize',16,'fontweight','bold');       
    
    % plot LIMA inset in figure
    ax4=axes('pos',[0.15 0.13 0.2 0.2]);
    set(ax4,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[]);
    cd([homepath,'../Imagery']);
    LIMA = imread('LIMA.jpg');
    imshow(LIMA);
    hold on; plot(1100,4270,'o','color','y','markersize',10,'linewidth',3);
    
if save_figure
    set(gcf,'InvertHardCopy','off'); % save colors as is
    cd([homepath,'../figures']);
    saveas(gcf,'StudyArea.png','png');
    cd([homepath,'../write-ups/Thesis/figures']);
    saveas(gcf,'1_StudyArea.png','png');
    disp('figure saved.');
end

%% Figure 2. Observed Conditions Time Series

close all;

cd([homepath,'inputs-outputs']);

save_figure = 0; % = 1 to save figure

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
termx = cl.xi(dsearchn([cl.Xi cl.Yi],[termX' termY']));
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
            plot(cl.xi./10^3,h(i).surface,'linewidth',1.5,...
                'color',col(col_h(i),:),'HandleVisibility','off');
        end
    end
    plot(cl.xi./10^3,hb,'-k','linewidth',1.5,'displayname','Bed');
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
        plot(cl.xi./10^3,U(i).speed.*3.1536e7,'linewidth',1.5,'color',col(col_h(i),:));
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
    colormap(col); c = colorbar('Limits',[2002 2019],'Ticks',[2002 2010 2019],'Position',[.92 .32 .03 .3410]);
    for i=1:length(SMB)
        plot(cl.xi(1:137)./10^3,SMB(i).smb_interp,'color',col(col_SMB(i),:),'linewidth',1.5);
    end
    % Add text label
    td = text((min(get(gca,'XLim')-0.25*(max(get(gca,'XLim')))-min(get(gca,'XLim')))),...
        (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' d ','edgecolor','k','fontsize',13,'fontweight','bold','linewidth',1.5);          

% Save figure
if save_figure
    cd([homepath,'../figures']);
    saveas(gcf,'observedTimeSeries.png','png');
    disp(['observedTimeSeries.png saved in: ',pwd]);

    cd([homepath,'../write-ups/Thesis/figures']);
    saveas(gcf,'observedTimeSeries.png','png');
    disp(['2_observedTimeSeries.png saved in: ',pwd]);
end
    
%% Figure 3. Sensitivity Tests

close all;

% load sensitivity test no change variables
cd([homepath,'inputs-outputs']);
load('Crane_sensitivityTests_noChange.mat');

cd([homepath,'scripts/modelingWorkflow/results']);
addpath([homepath,'../matlabFunctions/cmocean_v2.0/cmocean']);

save_figure = 0; % = 1 to save figure

% Set up figures
figure(10); clf % SMR
set(gcf,'Position',[100 300 1300 400]);
subplot(1,3,1); hold on; grid on; % a) geometry
    set(gca,'fontsize',16,'linewidth',2);
    xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m.a.s.l)'); 
    xlim([35 60]); ylim([-1200 600]);
    plot(x1./10^3,hb1,'-k','linewidth',2);
    % Add text label            
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
        (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' a ','edgecolor','k','fontsize',16,'fontweight','bold','linewidth',1.5);  
subplot(1,3,2); hold on; grid on; % b) speed
    set(gca,'fontsize',16,'linewidth',2);
    xlabel('Distance Along Centerline (km)'); ylabel('Speed (m a^{-1})');  
    xlim([0 60]); ylim([0 2200]);
    % Add text label            
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
        (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' b ','edgecolor','k','fontsize',16,'fontweight','bold','linewidth',1.5); 
subplot(1,3,3); hold on; grid on; % c) \Delta H
    set(gca,'fontsize',16,'linewidth',2);
    xlabel('\Delta SMR (m a^{-1})'); ylabel('\DeltaH (m)');
    xlim([-1.5 6.5]); ylim([-5 5]);
    % Add text label            
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
        (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' c ','edgecolor','k','fontsize',16,'fontweight','bold','linewidth',1.5);         
figure(11); clf % SMB
set(gcf,'Position',[200 100 1300 400]);
subplot(1,3,1); hold on; grid on; % c) geometry
    set(gca,'fontsize',16,'linewidth',2);
    xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m.a.s.l)'); 
    xlim([35 60]); ylim([-1200 600]);        
    plot(x1./10^3,hb1,'-k','linewidth',2);
    % Add text label            
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
        (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' d ','edgecolor','k','fontsize',16,'fontweight','bold','linewidth',1.5);         
subplot(1,3,2); hold on; grid on; % d) speed
    set(gca,'fontsize',16,'linewidth',2);
    xlabel('Distance Along Centerline (km)'); ylabel('Speed (m a^{-1})');     
    xlim([0 60]); ylim([0 2200]); grid on;
    % Add text label            
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
        (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' e ','edgecolor','k','fontsize',16,'fontweight','bold','linewidth',1.5); 
subplot(1,3,3); hold on; grid on; % c) \Delta H
    set(gca,'fontsize',16,'linewidth',2);
    xlabel('\Delta SMB (m a^{-1})'); ylabel('\DeltaH (m)');
    xlim([-6.5 1.5]); ylim([-40 15]);
    % Add text label            
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
        (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' f ','edgecolor','k','fontsize',16,'fontweight','bold','linewidth',1.5); 
    
% Loop through results, split into SMR & SMB change scenarios
files = dir('*geom.mat');
for i=1:length(files)
    if contains(files(i).name,'SMB0') && contains(files(i).name,'_geom.mat')
        files(i).change = -str2double(files(i).name(regexp(files(i).name,'R')+1:...
            regexp(files(i).name,'_S')-1)).*3.1536e7; 
        files(i).changeIn = 'SMR';
    elseif contains(files(i).name,'SMR0') && contains(files(i).name,'_geom.mat')
        files(i).change = str2double(files(i).name(regexp(files(i).name,'B')+1:...
            regexp(files(i).name,'_g')-1)).*3.1536e7;  
        files(i).changeIn = 'SMB';
    end
end
    
% Sort by changeIn and change
files = struct2table(files);
files = sortrows(files,[8,7]);
% save indices for SMB & SMR
Ismb = find(strcmp('SMB',table2array(files(:,8))));
Ismr = find(strcmp('SMR',table2array(files(:,8))));
files = table2struct(files);

col=cmocean('thermal',length(Ismb)+2); % color scheme for plotting

% SMR
for i=1:length(Ismr)
    load(files(Ismr(i)).name);
    figure(10); hold on;
    subplot(1,3,1); hold on; 
        % ice surface
        plot(x1(1:c2)./10^3,h2(1:c2),'color',col(i,:),'linewidth',2);
        % calving front
        plot([x1(c2) x1(c2)]./10^3,[h2(c2) h2(c2)-H2(c2)],'color',col(i,:),'linewidth',2);
        % floating bed
        plot(x1(gl2:c2)./10^3,h2(gl2:c2)-H2(gl2:c2),'color',col(i,:),'linewidth',2);
    subplot(1,3,2); hold on; 
        % ice surface speed
        plot(x1(1:c2)./10^3,U2(1:c2).*3.1536e7,'color',col(i,:),'linewidth',2);
    subplot(1,3,3); hold on; 
        % dH
        dH = nanmean(H2(1:c2))-nanmean(H1(1:c1)); 
        plot(files(Ismr(i)).change,dH,'o',...
            'markersize',10,'linewidth',2,'color',col(i,:));
        % colormap
        colormap(col); c=colorbar; caxis([-1 8]);   
        set(get(c,'label'),'String','\DeltaSMR (m a^{-1})');
end

% SMB
for i=1:length(Ismb)
    load(files(Ismb(i)).name);
    figure(11); hold on;
    subplot(1,3,1); hold on; 
        % ice surface
        plot(x2(1:c2)./10^3,h2(1:c2),'color',col(i,:),'linewidth',2);
        % calving front
        plot([x2(c2) x2(c2)]./10^3,[h2(c2)-H2(c2) h2(c2)],'color',col(i,:),'linewidth',2);
        % floating bed
        plot(x2(gl2:c2)./10^3,h2(gl2:c2)-H2(gl2:c2),'color',col(i,:),'linewidth',2);
    subplot(1,3,2); hold on; 
        % ice surface speed
        plot(x2(1:c2)./10^3,U2(1:c2).*3.1536e7,'color',col(i,:),'linewidth',2);
    subplot(1,3,3); hold on;
        % dH
        dH = nanmean(H2(1:c2))-nanmean(H1(1:c1)); 
        plot(files(Ismb(i)).change,dH,'o',...
            'markersize',10,'linewidth',2,'color',col(i,:));
        % colormap
        colormap(col); c=colorbar; caxis([-6 3]); 
        set(get(c,'label'),'String','\DeltaSMB (m a^{-1})');            
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
    saveas(gcf,'3_Crane_sensitivityTests_SMR.png','png');
    figure(11);
    saveas(gcf,'3_Crane_sensitivityTests_SMB.png','png');    
    disp(['figure saved in: ',pwd]);  
    
end


