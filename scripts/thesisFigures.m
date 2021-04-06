%% Create Figures for MS Thesis
% RKA 2021

% Figure 1. Map of the study area
% Figure 2. Observed Conditions Time Series
% Figure 3. Beta Solution
% Figure 4. Sensitivity Tests
% Figure 5. Regional glacier terminus time series

close all; clear all;

% Define home path
homepath = '/Users/raineyaberle/Desktop/Research/CraneGlacier_flowlinemodeling/';

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

close all;

save_figure = 0;    % = 1 to save figure
fontsize = 16;      % font size
markersize = 8;    % marker size
linewidth = 2;      % line width

% Add paths to functions
    addpath([homepath,'../matlabFunctions']);
    addpath([homepath,'../matlabFunctions/cmocean_v2.0/cmocean']);
    addpath([homepath,'../matlabFunctions/AntarcticMappingTools_v5.17/AntarcticMappingTools']);

% Specify xticks/yticks coordinates and axes limits
    xticks = linspace(-2.455e6,-2.375e6,5); yticks=linspace(1.21e6,1.29e6,5);
    xlimits=[xticks(1) xticks(end)]; ylimits = [yticks(1) yticks(end)];
    
% load & display Landsat image
    cd([homepath,'../Imagery/LC08_L1GT_218106_20191013_20191018_01_T2']);
    landsat = dir('*B8.TIF');
    [LS.im,LS.R] = readgeoraster(landsat.name); [LS.ny,LS.nx] = size(LS.im);
    % polar stereographic coordinates of image boundaries
    LS.x = linspace(min(LS.R.XWorldLimits),max(LS.R.XWorldLimits),LS.nx); 
    LS.y = linspace(min(LS.R.YWorldLimits),max(LS.R.YWorldLimits),LS.ny);
    % display image on ax1
    figure(1); clf; hold on; 
    set(gcf,'position',[452 75 1000 800],'color','w');
    ax1=gca; 
    im1 = imagesc(ax1,LS.x,LS.y,flipud(LS.im*1.1)); colormap('gray');
    % set axes properties
    ax1.XTick=xticks; ax1.XTickLabel=string(xticks./10^3);
    ax1.YTick=yticks; ax1.YTickLabel=string(yticks./10^3);
    set(ax1,'YDir','normal','XLim',xlimits','YLim',ylimits,'FontSize',14,...
        'linewidth',2,'fontsize',fontsize,'Position',[0.165 0.08 0.67 0.84]); 
    xlabel('Easting (km)'); ylabel('Northing (km)');
    
% load & display REMA as contours
    cd([homepath,'../Imagery']);
    [REMA.im,REMA.R] = readgeoraster('REMA_clipped.tif');
    REMA.im(REMA.im==-9999)=NaN;
    [REMA.ny,REMA.nx] = size(REMA.im);
    % polar stereographic coordinates of image boundaries
    REMA.x = linspace(min(REMA.R.XWorldLimits),max(REMA.R.XWorldLimits),REMA.nx); 
    REMA.y = linspace(min(REMA.R.YWorldLimits),max(REMA.R.YWorldLimits),REMA.ny);
    % display contours on ax1
    hold on; contour(ax1,REMA.x,REMA.y,flipud(REMA.im),0:500:2000,'-k','linewidth',2,'ShowText','on');
    
% load & display ITS_LIVE velocity map
    cd([homepath,'data/velocities']);
    v.v = ncread('ANT_G0240_2017.nc','v'); v.v(v.v==-3267)=NaN;
    v.x = ncread('ANT_G0240_2017.nc','x'); v.y = ncread('ANT_G0240_2017.nc','y');
    % display velocity map on ax2
    ax2 = axes; 
    im2=imagesc(ax2,v.x,v.y,v.v'); colormap(cmocean('haline'));
    im2.AlphaData=0.5; caxis([0 1100]);
    % set axes properties
    ax2.XLim=xlimits; ax2.YLim=ylimits; 
    set(ax2,'YDir','normal','Visible','off','fontsize',12,'XTick',[],'YTick',[]); hold on; 
    % add legend
    l=legend('Position',[0.23 0.8 0.12 0.07],'Color',[175 175 175]/255,...
        'fontname','arial','fontsize',fontsize);
    
% Add lat lon coordinates to top and right axes
    latlon = graticuleps(-66:0.25:-64,-64:0.5:-62,'color',[175 175 175]/255,'linewidth',1,...
        'HandleVisibility','off');
    ax2.Position=[0.08 0.08 0.84 0.84];     
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','62.00^oS','textcolor',[175 175 175]/255,...
        'HeadStyle','none','LineStyle', 'none','Position',[.9 .64 0 0],'FontSize',fontsize,'FontName','Arial');     
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','62.50^oS','textcolor',[175 175 175]/255, ...
        'HeadStyle','none','LineStyle', 'none','Position',[.9 .36 0 0],'FontSize',fontsize,'FontName','Arial');     
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','65.25^oW','textcolor',[175 175 175]/255, ...
        'HeadStyle','none','LineStyle', 'none','Position',[.64 .94 0 0],'FontSize',fontsize,'FontName','Arial');     
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','65.00^oW','textcolor',[175 175 175]/255, ...
        'HeadStyle','none','LineStyle', 'none','Position',[.36 .94 0 0],'FontSize',fontsize,'FontName','Arial'); 
    
% plot OIB flights
    cd([homepath,'../figures/1_StudyArea']);
    OIBfiles = dir('IRMCR2*.csv');
    for i=1:length(OIBfiles)
        OIB.file = readmatrix(OIBfiles(i).name);
        [OIB.x,OIB.y] = wgs2ps(OIB.file(:,2),OIB.file(:,1),'StandardParallel',-71,'StandardMeridian',0);
        if i==1
            hold on; plot(ax2,OIB.x,OIB.y,'color','w',...
                'displayname','NASA OIB','linewidth',linewidth-1);
        else
            hold on; plot(ax2,OIB.x,OIB.y,'color','w',...
                'HandleVisibility','off','linewidth',linewidth-1);            
        end
    end
    
% plot centerline points
    plot(cl.X,cl.Y,'o','color','m','markerfacecolor','m',...
        'markeredgecolor','k','markersize',markersize,'displayname','Centerline');
    
% plot centerline distance markers
    dist = 0:10e3:cl.x(end); Idist = dsearchn(cl.x',dist'); % ticks every 10 km
    plot(cl.X(Idist),cl.Y(Idist),'o','markerfacecolor',[253,224,221]/255,...
        'markeredgecolor','k','markersize',markersize,'handlevisibility','off');
    labels = string(dist./10^3);
    text(cl.X(Idist)+1500,cl.Y(Idist),labels,'color',[253,224,221]/255,...
        'fontweight','bold','fontname','arial','fontsize',fontsize);

% add colorbar
    c = colorbar('Position',[0.29 0.55 0.02 0.15],'fontsize',fontsize,'fontweight',...
        'bold','color','w','fontname','arial');
    set(get(c,'Title'),'String','Speed (m a^{-1})','color','w','fontname',...
        'arial','Position',[5 125 0]);

% insert text labels
    % Crane
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Crane', ...
       'HeadStyle','none','LineStyle', 'none', 'TextRotation',70,'Position',[.54 0.65 0 0],...
       'FontSize',fontsize-1,'FontName','Arial','fontweight','bold','TextColor',[200 200 200]/255);
    % Former Larsen B Ice Shelf
    %txt = sprintf('Former Larsen \n    B Ice Shelf');
    %text(-2.405e6,1.282e6,txt,'color',[200 200 200]/255,'fontsize',12,'fontweight','bold');
    % Jorum
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Jorum', ...
        'HeadStyle','none','LineStyle', 'none', 'TextRotation',35,'Position',[.52 .85 0 0],...
        'FontSize',fontsize-1,'FontName','Arial','fontweight','bold','TextColor',[200 200 200]/255);    
    % Flask
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Flask', ...
       'HeadStyle','none','LineStyle', 'none', 'TextRotation',55,'Position',[.77 .14 0 0],...
       'FontSize',fontsize-1,'FontName','Arial','fontweight','bold','TextColor',[200 200 200]/255);    
    % Mapple
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Mapple', ...
        'HeadStyle','none','LineStyle', 'none', 'TextRotation',66,'Position',[.625 .62 0 0],...
        'FontSize',fontsize-1,'FontName','Arial','fontweight','bold','TextColor',[200 200 200]/255);    
    % Melville
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Melville', ...
        'HeadStyle','none','LineStyle', 'none', 'TextRotation',70,'Position',[.68 .6 0 0],...
        'FontSize',fontsize-1,'FontName','Arial','fontweight','bold','TextColor',[200 200 200]/255);           
    % Pequod
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Pequod', ...
        'HeadStyle','none','LineStyle', 'none', 'TextRotation',66,'Position',[.73 .59 0 0],...
        'FontSize',fontsize-1,'FontName','Arial','fontweight','bold','TextColor',[200 200 200]/255);               
    % Starbuck
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Starbuck', ...
        'HeadStyle','none','LineStyle', 'none', 'TextRotation',38,'Position',[.78 .465 0 0],...
        'FontSize',fontsize-1,'FontName','Arial','fontweight','bold','TextColor',[200 200 200]/255);   
    % Tributary A
    txt = sprintf('A');
        text(-2.413e6,1.232e6,txt,'color',[200 200 200]/255,'fontsize',fontsize-1,'fontweight','bold','fontname','arial');       
    % Tributary B
    txt = sprintf('B');
        text(-2.415e6,1.237e6,txt,'color',[200 200 200]/255,'fontsize',fontsize-1,'fontweight','bold','fontname','arial');    
    % Tributary C
    txt = sprintf('C');
        text(-2.417e6,1.251e6,txt,'color',[200 200 200]/255,'fontsize',fontsize-1,'fontweight','bold','fontname','arial');       
    
    % plot LIMA inset in figure
    ax4=axes('pos',[0.2 0.12 0.2 0.2]);
    set(ax4,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[]);
    cd([homepath,'../Imagery']);
    LIMA = imread('LIMA.jpg');
    imshow(LIMA);
    hold on; plot(1100,4270,'o','color','y','markersize',10,'linewidth',3);
    
if save_figure
    set(gcf,'InvertHardCopy','off'); % save colors as is
    cd([homepath,'../figures']);
    saveas(gcf,'fig1.1_studyArea.png','png');
    cd([homepath,'../write-ups/Thesis/figures']);
    saveas(gcf,'fig1.1_studyArea.png','png');
    cd([homepath,'../write-ups/Thesis/MSThesis_RaineyAberle_LateX']);
    saveas(gcf,'fig1.1_studyArea.png','png');
    disp('figure 1 saved.');    
end

%% Figure 2. Observed Conditions Time Series

close all;

cd([homepath,'inputs-outputs']);

save_figure = 1; % = 1 to save figure
fontsize = 16; % fontsize for plots
linewidth = 2.5; % line width for plots
markersize = 10; % marker size for plots
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
set(gcf,'Position',[100 100 1000 600]);
ax1=axes('Position',[0.1 0.6 0.35 0.35],'linewidth',2,'fontsize',fontsize); % geometry
    hold on; grid on; ylabel('Elevation (m)'); 
    xlim([0 60]); ylim([-1200 1200]);
    for i=1:length(h)
        if i==16
        else
            plot(ax1,cl.xi./10^3,h(i).surface,'linewidth',linewidth,...
                'color',col(col_h(i),:),'HandleVisibility','off');
        end
    end
    plot(ax1,cl.xi./10^3,hb,'-k','linewidth',linewidth,'displayname','Bed');
    % Add text label
    text(55.5,(max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' a ','edgecolor','k','fontsize',fontsize,'fontweight','bold','linewidth',linewidth-1,'backgroundcolor','w');    
ax2=axes('Position',[0.55 0.6 0.35 0.35],'linewidth',2,'fontsize',fontsize); % speed
    hold on; grid on; ylabel('Speed (m a^{-1})');
    xlim([0 60]);
    for i=1:length(U)
        plot(ax2,cl.xi./10^3,U(i).speed.*3.1536e7,'linewidth',linewidth,'color',col(col_h(i),:));
    end
    % Add text label
    text(55.5,(max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' b ','edgecolor','k','fontsize',fontsize,'fontweight','bold','linewidth',linewidth-1,'backgroundcolor','w');      
ax3=axes('Position',[0.1 0.1 0.35 0.35],'linewidth',2,'fontsize',fontsize); % terminus position
    hold on; grid on; xlabel('Distance Along Centerline (km)'); ylabel('Date');
    xlim([0 60]);
    for i=1:length(termx)
        plot(ax3,termx(i)./10^3,termdate(i),'ok','MarkerFaceColor',col(col_term(i),:),'markersize',markersize);
    end
    % Add text label
    text(55.5,(max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' c ','edgecolor','k','fontsize',fontsize,'fontweight','bold','linewidth',linewidth-1,'backgroundcolor','w');     
ax4=axes('Position',[0.55 0.1 0.35 0.35],'linewidth',2,'fontsize',fontsize); % speed
    hold on; grid on; xlabel('Distance Along Centerline (km)'); ylabel('Mean Annual SMB (m a^{-1})');
    xlim([0 60]); set(gca,'clim',[2002 2019]); 
    % Add colorbar
    colormap(col); 
    c = colorbar('Limits',[2002 2019],'Ticks',[2002 2010 2019],'Position',[.92 .32 .03 .3410],...
        'fontsize',fontsize);
    for i=1:length(SMB)
        plot(ax4,cl.xi(1:137)./10^3,SMB(i).smb_interp,'color',col(col_SMB(i),:),'linewidth',linewidth);
    end
    % Add text label
    text(55.5,(max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' d ','edgecolor','k','fontsize',fontsize,'fontweight','bold','linewidth',linewidth-1,'backgroundcolor','w');          

% Save figure
if save_figure
    cd([homepath,'../figures']);
    saveas(gcf,'fig3.1_centerlineObs.png','png');
    cd([homepath,'../write-ups/Thesis/figures']);
    saveas(gcf,'fig3.1_centerlineObs.png','png');
    cd([homepath,'../write-ups/Thesis/MSThesis_RaineyAberle_LateX']);
    saveas(gcf,'fig3.1_centerlineObs.png','png');    
    disp('figure 2 saved.');
end

%% Figure 3. Beta Solution

close all;

save_figure = 1;                % = 1 to save figure
fonsize = 18;                   % font size
fontname = 'Arial';             % font name
linewidth = 2.5;                % line width
markersize = 30;                % marker size

% load betaSolutions
cd([homepath,'inputs-outputs/']);
beta = load('betaSolutions.mat').betaSolutions;
RMSE = zeros(1,length(beta)); dx0 = zeros(1,length(beta));
for i=1:length(beta)
    RMSE(i) = beta(i).RMSE;
    dx0(i) = beta(i).dx;
end

% load stress coupling length results
load('Crane_SCL_results.mat');

% load observed speed
x0 = load('Crane_flowline_initialization.mat').x0;
U0 = load('Crane_flowline_initialization.mat').U0;
c0 = load('Crane_flowline_initialization.mat').c0;

% determine best beta spatial resolution
Ibest = find(abs(gradient(gradient(RMSE)))<0.1e-5,1,'first');
if length(Ibest)>1 % if lowest RMSE exists for multiple dx, use lowest dx
    Ibest(2:end)=[];
end

% plot results
figure(3); clf
set(gcf,'Position',[272 143 1063 654]);
ax1=axes; set(gca,'position',[0.08 0.58 0.4 0.4]);
    set(gca,'fontsize',fontsize,'fontname',fontname,'linewidth',linewidth);
    xlabel('Distance Along Centerline (km)'); ylabel('\beta (s^{1/m} m^{-1/m})'); 
    hold on; grid on; xlim([0 45]);
    plot(beta(Ibest).x./10^3,beta(Ibest).beta,'linewidth',2); 
    text(42,(max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' a ','edgecolor','k','fontsize',fontsize,'fontweight','bold','linewidth',1.5,'backgroundcolor','w');     
ax2=axes; set(gca,'position',[0.58 0.58 0.4 0.4]);
    set(gca,'fontsize',fontsize,'fontname',fontname,'linewidth',2);
    xlabel('dx (km)'); ylabel('RMSE (m a^{-1})'); 
    hold on; plot(dx0./10^3,RMSE.*3.1536e7,'.','markersize',markersize); grid on;
    plot(dx0(Ibest)/10^3,RMSE(Ibest)*3.1536e7,'*','markersize',markersize-10,'linewidth',linewidth);  
    text(1.9,(max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' b ','edgecolor','k','fontsize',fontsize,'fontweight','bold','linewidth',1.5,'backgroundcolor','w');     
ax3=axes; set(gca,'position',[0.08 0.08 0.4 0.4]);
    set(gca,'fontsize',fontsize,'fontname',fontname,'linewidth',2);
    xlabel('Distance Along Centerline (km)'); ylabel('U (m a^{-1})');
    hold on; grid on; legend('Location','best'); xlim([0 45]);
    plot(x0(1:c0)./10^3,U0(1:c0).*3.1536e7,'k','displayname','observed','linewidth',linewidth);
    plot(beta(Ibest).x./10^3,beta(Ibest).U.*3.1536e7,'--k','displayname','solution','linewidth',2); 
    text(42,(max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' c ','edgecolor','k','fontsize',fontsize,'fontweight','bold','linewidth',1.5,'backgroundcolor','w');     
ax4=axes; set(gca,'position',[0.58 0.08 0.4 0.4]);
    set(gca,'fontsize',fontsize,'fontname',fontname,'linewidth',linewidth);
    xlabel('Distance Along Centerline (km)'); ylabel('R_{xx} (kPa)');
    hold on; grid on; 
    plot(x./10^3,Rxx./10^3,'linewidth',linewidth); 
    ax=get(gca); % get current axes
    for i=1:length(ITRxx)
        plot([x(ITRxx(i)) x(ITRxx(i))]./10^3,[ax.YLim(1) ax.YLim(2)],'--k','linewidth',linewidth-1.5);
    end  
    xlim([0 45]);
    text(5,-1000,['\Deltax_{mean} = ',num2str(round(mean(diff(x(ITRxx))))/10^3),'km'],...
        'edgecolor','k','fontsize',fontsize,'fontweight','bold','linewidth',linewidth-0.5,...
        'backgroundcolor','w','color','k','fontname',fontname);    
    text(42,(max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' d ','edgecolor','k','fontsize',fontsize-1,'fontweight','bold','linewidth',linewidth-0.5,...
        'backgroundcolor','w','fontname',fontname);    
    text(28,-200,'\Deltax','fontsize',13,'fontweight','bold','color','k');
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','}', ...
       'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[0.835 0.25 0 0],...
       'FontSize',28,'FontName',fontname,'TextColor','k');
        
if save_figure
    figure(3);
    % save in figures folder
    cd([homepath,'../figures/']);
    saveas(gcf,'fig3.2_betaSolution.png','png');  
    % save in thesis figures folder
    cd([homepath,'../write-ups/Thesis/figures/']);
    saveas(gcf,'fig3.2_betaSolution.png','png');
    % save in thesis latex folder
    cd([homepath,'../write-ups/Thesis/MSThesis_RaineyAberle_LateX/']);
    saveas(gcf,'fig3.2_betaSolution.png','png');    
    disp('figure 3 saved.');    
end
    
%% Figure 4. Sensitivity Tests

close all;

save_figure = 0;    % = 1 to save figure
fontsize = 16;      % font size
fontname = 'Arial'; % font name
linewidth = 2;      % line width
markersize = 8;    % marker size

% load sensitivity test no change variables
cd([homepath,'inputs-outputs']);
load('Crane_2100_noChange.mat');
% load glacier width
W = load('Crane_flowline_initialization.mat').W0;
% load initial fwd
fwd0 = load('fwdbest.mat').fwdbest;

cd([homepath,'scripts/modelingWorkflow/results']);
addpath([homepath,'../matlabFunctions/cmocean_v2.0/cmocean']);

% set up axes
figure(4); clf 
set(gcf,'Position',[100 200 1000 1000]);
% SMR
ax1=axes('position',[0.08 0.73 0.25 0.25]); hold on; grid on; % a) geometry
    set(gca,'fontsize',fontsize,'linewidth',2,'YTick',-1500:500:1500);
    ylabel('Elevation (m)'); 
    xlim([35 60]); ylim([-1200 600]);
    plot(x1./10^3,hb1,'-k','linewidth',2);
    % add text label            
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
        (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' a ','edgecolor','k','backgroundcolor','w','fontsize',fontsize,...
        'fontweight','bold','linewidth',linewidth-1);  
ax2=axes('position',[0.4 0.73 0.25 0.25]); hold on; grid on; % b) speed
    set(gca,'fontsize',fontsize,'linewidth',2,'YTick',0:500:2500);        
    ylabel('U (m a^{-1})');  
    xlim([0 60]); ylim([0 3000]);
    % add text label            
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
        (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' b ','edgecolor','k','backgroundcolor','w','fontsize',fontsize,...
        'fontweight','bold','linewidth',linewidth-1); 
ax3=axes('position',[0.72 0.73 0.21 0.25]); hold on; grid on; % c) \Delta H
    set(gca,'fontsize',fontsize,'linewidth',2);
    xlabel('\DeltaSMR (m a^{-1})'); ylabel('\DeltaH_{\mu} (m)');
    xlim([-1.5 6.5]); ylim([-50 20]);               
    % add text label            
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
        (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' c ','edgecolor','k','backgroundcolor','w','fontsize',fontsize,...
        'fontweight','bold','linewidth',linewidth-1);     
% SMB
ax4=axes('position',[0.08 0.4 0.25 0.25]); hold on; grid on; % d) geometry
    set(gca,'fontsize',fontsize,'linewidth',2,'YTick',-1500:500:1500);
    ylabel('Elevation (m)'); 
    xlim([35 60]); ylim([-1200 600]);
    plot(x1./10^3,hb1,'-k','linewidth',2);
    % add text label            
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
        (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' d ','edgecolor','k','backgroundcolor','w','fontsize',fontsize,...
        'fontweight','bold','linewidth',linewidth-1);  
ax5=axes('position',[0.4 0.4 0.25 0.25]); hold on; grid on; % e) speed
    set(gca,'fontsize',fontsize,'linewidth',2,'YTick',0:500:2500);        
    ylabel('U (m a^{-1})');  
    xlim([0 60]); ylim([0 3000]);
    % add text label            
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
        (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' e ','edgecolor','k','backgroundcolor','w','fontsize',fontsize,...
        'fontweight','bold','linewidth',linewidth-1); 
ax6=axes('position',[0.72 0.4 0.21 0.25]); hold on; grid on; % f) \Delta H
    set(gca,'fontsize',fontsize,'linewidth',2);
    xlabel('\DeltaSMB (m a^{-1})'); ylabel('\DeltaH_{\mu} (m)');
    xlim([-6.5 1.5]); ylim([-50 20]);               
    % add text label            
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
        (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' f ','edgecolor','k','backgroundcolor','w','fontsize',fontsize,...
        'fontweight','bold','linewidth',linewidth-1);     
% fwd    
ax7=axes('position',[0.08 0.07 0.25 0.25]); hold on; grid on; % g) geometry
    set(gca,'fontsize',fontsize,'linewidth',2,'YTick',-1500:500:1500);
    xlabel('Distance Along Centerline (m)'); ylabel('Elevation (m)'); 
    xlim([35 60]); ylim([-1200 600]);
    plot(x1./10^3,hb1,'-k','linewidth',2);
    % add text label            
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
        (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' g ','edgecolor','k','backgroundcolor','w','fontsize',fontsize,...
        'fontweight','bold','linewidth',linewidth-1);  
ax8=axes('position',[0.4 0.07 0.25 0.25]); hold on; grid on; % h) speed
    set(gca,'fontsize',fontsize,'linewidth',2,'YTick',0:500:2500);        
    xlabel('Distance Along Centerline (m)'); ylabel('U (m a^{-1})');  
    xlim([0 60]); ylim([0 3000]);
    % add text label            
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
        (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' h ','edgecolor','k','backgroundcolor','w','fontsize',fontsize,...
        'fontweight','bold','linewidth',linewidth-1); 
ax9=axes('position',[0.72 0.07 0.21 0.25]); hold on; grid on; % i) \Delta H
    set(gca,'fontsize',fontsize,'linewidth',2);
    xlabel('fwd (m)'); ylabel('\DeltaH_{\mu} (m)');
    xlim([20 28]); ylim([-50 20]);               
    % add text label            
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
        (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' i ','edgecolor','k','backgroundcolor','w','fontsize',fontsize,...
        'fontweight','bold','linewidth',linewidth-1);         
    
% loop through results, split into SMR, SMB, and fwd change scenarios
files = dir('*geom.mat');
for i=1:length(files)
    if contains(files(i).name,'SMB0') && contains(files(i).name,'_geom.mat')
        if contains(files(i).name,'SMR0') && contains(files(i).name,'SMB0')
            files(i).change = 0; 
            files(i).changeIn = 'SMR'; 
            files(length(files)+1).change = 0;
            files(length(files)).changeIn = 'SMB';
            files(length(files)).name = files(i).name; 
        else
            files(i).change = -str2double(files(i).name(regexp(files(i).name,'R')+1:...
                regexp(files(i).name,'_S')-1)).*3.1536e7; 
        end
        files(i).changeIn = 'SMR';                    
    elseif contains(files(i).name,'SMR0') && contains(files(i).name,'_geom.mat')
        files(i).change = str2double(files(i).name(regexp(files(i).name,'B')+1:...
            regexp(files(i).name,'_geom')-1)).*3.1536e7; 
        files(i).changeIn = 'SMB';  
    elseif contains(files(i).name,'fwd') && contains (files(i).name,'_geom.mat')
        files(i).change = str2double(files(i).name(regexp(files(i).name,'d_')+2:...
            regexp(files(i).name,'m')-1))-fwd0; % m 
        files(i).changeIn = 'fwd';         
    end
end
    
% sort by changeIn and change
files = struct2table(files);
files = sortrows(files,[8,7]);
% save indices for SMB & SMR
Ismb = find(strcmp('SMB',table2array(files(:,8))));
Ismr = find(strcmp('SMR',table2array(files(:,8))));
Ifwd = find(strcmp('fwd',table2array(files(:,8))));
files = table2struct(files);

col=cmocean('thermal',length(Ismb)+3); % color scheme for plotting
col(1,:) = []; col(end-1:end,:)=[]; % don't use end member colors

% SMR
dL = zeros(length(Ismb),1);
dH = zeros(length(Ismb),1);
dU = zeros(length(Ismb),1);
dV = zeros(length(Ismb),1);
dM = zeros(length(Ismb),1);
for i=1:length(Ismr)
    load(files(Ismr(i)).name);
    figure(4); hold on;
        % ice surface
        plot(ax1,x2(1:c2)./10^3,h2(1:c2),'color',col(i,:),'linewidth',linewidth);
        % calving front
        plot(ax1,[x2(c2) x2(c2)]./10^3,[h2(c2) h2(c2)-H2(c2)],'color',col(i,:),'linewidth',linewidth);
        % floating bed
        plot(ax1,x2(gl2:c2)./10^3,h2(gl2:c2)-H2(gl2:c2),'color',col(i,:),'linewidth',linewidth);
        % ice surface speed
        plot(ax2,x2(1:c2)./10^3,U2(1:c2).*3.1536e7,'color',col(i,:),'linewidth',linewidth);
        % dH
        dH(i) = nanmean(H2(1:c2))-nanmean(H1(1:c1)); 
        plot(ax3,files(Ismr(i)).change,dH(i),'o',...
            'markersize',10,'linewidth',linewidth,'color',col(i,:));
        % colormap
        colormap(col);
        cb1=colorbar('position',[0.94 0.73 0.015 0.25],'fontname',fontname,'Ticks',0:1/7:1,'TickLabels',string(-1:6));
        set(get(cb1,'label'),'String','\DeltaSMR (m a^{-1})','fontsize',fontsize-3);
    % save results 
    dU(i) = (nanmean(U2(1:c2))-nanmean(U1(1:c1))).*3.1536e7; % m/a
    dL(i) = x2(c2)-x1(c1); % m   
    dV(i) = dH(i)*nanmean(W(1:135))*x1(c1); 
    dM(i) = dV(i)*917*10e-12; % metric Gt
end
% create table to store result quantities
T_smr = table(dL,dU,dH,dV,dM);

% SMB
dL = zeros(length(Ismb),1);
dH = zeros(length(Ismb),1);
dU = zeros(length(Ismb),1);
dV = zeros(length(Ismb),1);
dM = zeros(length(Ismb),1);
for i=1:length(Ismb)
    load(files(Ismb(i)).name);
    figure(4); hold on;
        % ice surface
        plot(ax4,x2(1:c2)./10^3,h2(1:c2),'color',col(i,:),'linewidth',linewidth);
        % calving front
        plot(ax4,[x2(c2) x2(c2)]./10^3,[h2(c2)-H2(c2) h2(c2)],'color',col(i,:),'linewidth',linewidth);
        % floating bed
        plot(ax4,x2(gl2:c2)./10^3,h2(gl2:c2)-H2(gl2:c2),'color',col(i,:),'linewidth',linewidth);
        % ice surface speed
        plot(ax5,x2(1:c2)./10^3,U2(1:c2).*3.1536e7,'color',col(i,:),'linewidth',linewidth);
        % dH
        dH(i) = nanmean(H2(1:c2))-nanmean(H1(1:c1)); % m
        plot(ax6,files(Ismb(i)).change,dH(i),'o',...
            'markersize',10,'linewidth',linewidth,'color',col(i,:));
        % colormap
        colormap(col);
        cb2=colorbar('position',[0.94 0.4 0.015 0.25],'fontname',fontname,'Ticks',0:1/7:1,'TickLabels',string(-6:1)); 
        set(get(cb2,'label'),'String','\DeltaSMB (m a^{-1})','fontsize',fontsize-2);    
    % save results 
    dU(i) = (nanmean(U2(1:c2))-nanmean(U1(1:c1))).*3.1536e7; % m/a
    dL(i) = x2(c2)-x1(c1); % m
    dV(i) = dH(i)*nanmean(W(1:135))*x1(c1); 
    dM(i) = dV(i)*917*10e-12; % metric Gt
end
% create table to store result quantities
T_smb = table(dL,dU,dH,dV,dM);

% fwd
dL = zeros(length(Ismb),1);
dH = zeros(length(Ismb),1);
dU = zeros(length(Ismb),1);
dV = zeros(length(Ismb),1);
dM = zeros(length(Ismb),1);
col = cmocean('thermal',length(Ifwd)+3); col(1,:)=[];col(end-1:end,:)=[];
for i=1:length(Ifwd)
    load(files(Ifwd(i)).name);
    figure(4); hold on;
        % ice surface
        plot(ax7,x2(1:c2)./10^3,h2(1:c2),'color',col(i,:),'linewidth',linewidth);
        % calving front
        plot(ax7,[x2(c2) x2(c2)]./10^3,[h2(c2)-H2(c2) h2(c2)],'color',col(i,:),'linewidth',linewidth);
        % floating bed
        plot(ax7,x2(gl2:c2)./10^3,h2(gl2:c2)-H2(gl2:c2),'color',col(i,:),'linewidth',linewidth);
        % ice surface speed
        plot(ax8,x2(1:c2)./10^3,U2(1:c2).*3.1536e7,'color',col(i,:),'linewidth',linewidth);
        % dH
        dH(i) = nanmean(H2(1:c2))-nanmean(H1(1:c1)); % m
        plot(ax9,files(Ifwd(i)).change+fwd0,dH(i),'o',...
            'markersize',10,'linewidth',linewidth,'color',col(i,:));
        % colormap
        colormap(col);
        cb3=colorbar('position',[0.94 0.08 0.015 0.25],'fontname',fontname,'Ticks',0:1/6:1,'TickLabels',string(-3:3)); 
        set(get(cb3,'label'),'String','\Deltafwd (m)','fontsize',fontsize-3);    
    % save results 
    dU(i) = (nanmean(U2(1:c2))-nanmean(U1(1:c1))).*3.1536e7; % m/a
    dL(i) = x2(c2)-x1(c1); % m
    dV(i) = dH(i)*nanmean(W(1:135))*x1(c1); 
    dM(i) = dV(i)*917*10e-12; % metric Gt
end
% create table to store result quantities
T_fwd = table(dL,dU,dH,dV,dM);

% save figures 
if save_figure
    figure(4);
    % save in figures folder
    cd([homepath,'../figures/']);
    saveas(gcf,'fig4.1_sensitivityTests.png','png');  
    % save in thesis figures folders
    cd([homepath,'../write-ups/Thesis/figures/']);
    saveas(gcf,'sensitivityTests.png','png');
    cd([homepath,'../write-ups/Thesis/MSThesis_RaineyAberle_LateX']);
    saveas(gcf,'fig4.1_sensitivityTests.png','png');    
    disp('figure 4 saved.');
end

%% Figure 5. Regional Glacier Terminus Time Series

close all;

cd([homepath,'data/terminus/regional/']);

save_figure = 0; % = 1 to save figure
fontsize = 16; 
fontname = 'Arial';
linewidth = 3; 

% add path to matlab functions
addpath([homepath,'matlabFunctions/']);

% load Landsat images
cd([homepath,'data/Imagery/']);
% LSA 
landsatA = dir('*216105*B8.TIF');
    [LSA.im,LSA.R] = readgeoraster(landsatA.name); [LSA.ny,LSA.nx] = size(LSA.im);
    % polar stereographic coordinates of image boundaries
    LSA.x = linspace(min(LSA.R.XWorldLimits),max(LSA.R.XWorldLimits),LSA.nx); 
    LSA.y = linspace(min(LSA.R.YWorldLimits),max(LSA.R.YWorldLimits),LSA.ny);
% LSB
landsatB = dir('*218106*B8.TIF');
    [LSB.im,LSB.R] = readgeoraster(landsatB.name); [LSB.ny,LSB.nx] = size(LSB.im);
    % polar stereographic coordinates of image boundaries
    LSB.x = linspace(min(LSB.R.XWorldLimits),max(LSB.R.XWorldLimits),LSB.nx); 
    LSB.y = linspace(min(LSB.R.YWorldLimits),max(LSB.R.YWorldLimits),LSB.ny);

% load LIMA
LIMA = imread('LIMA.jpg');

% load terminus coordinate shapefiles
cd([homepath,'data/terminus/regional/']);

% set up figure, subplots, and plot
figure(6); clf; 
set(gcf,'Position',[100 200 1000 1000]); 
% a) Drygalski
cd('Drygalski'); % enter folder 
    files = dir('*.shp'); % load file
    % loop through files to grab centerline
    for j=1:length(files)
        if contains(files(j).name,'CL')
            file = shaperead(files(j).name);
            cl.lon = file.X; % Lon
            cl.lat = file.Y; % Lat
            [cl.X,cl.Y] = wgs2ps(cl.lon,cl.lat,'StandardParallel',-71,'StandardMeridian',0);
            % define x as distance along centerline
            x = zeros(1,length(cl.X));
            for k=2:length(x)
                x(k) = sqrt((cl.X(k)-cl.X(k-1))^2+(cl.Y(k)-cl.Y(k-1))^2)+x(k-1); %(m)
            end
        end
    end
    % initialize f
    clear f; f(1).lon = NaN; f(1).lat = NaN; f(1).date = NaN;
    % loop through files to grab terminus positions
    for j=1:length(files)
        if ~contains(files(j).name,'CL')
            file = shaperead(files(j).name);
            % loop through file
            for k=1:length(file)
                f(length(f)+1).lon = file(k).X; % Lon
                f(length(f)).lat = file(k).Y; % Lat
                % convert to polar stereographic coordinates
                [f(length(f)).X,f(length(f)).Y] = wgs2ps(f(length(f)).lon,f(length(f)).lat,'StandardParallel',-71,'StandardMeridian',0);
                f(length(f)).x = x(dsearchn([cl.X' cl.Y'],[nanmean(f(length(f)).X) nanmean(f(length(f)).Y)])); % x
            end
        end 
    end
    % plot
    col = parula(length(f)); % color scheme for plotting
    ax1=axes('position',[0.08 0.58 0.23 0.4]); hold on; grid on;      
        set(gca,'fontname',fontname,'fontsize',fontsize);
        xlabel('Easting (km)'); ylabel('Northing (km)');
        % plot landsat image
        colormap(ax1,'gray');
        imagesc(LSA.x/10^3,LSA.y/10^3,flipud(LSA.im*1.1)); 
        % plot centerline
        plot(cl.X(4:7)/10^3,cl.Y(4:7)/10^3,'k','linewidth',linewidth);
        annotation('textarrow',[0.236 0.237],[0.85 0.8516],'headwidth',18,'headlength',16);
        % plot terminus positions
        for j=1:length(f)
            plot(f(j).X/10^3,f(j).Y/10^3,'color',col(j,:),'linewidth',linewidth);
        end
        axis tight; text((ax1.XLim(2)-ax1.XLim(1))*0.88+ax1.XLim(1),(max(ax1.YLim)-min(ax1.YLim))*0.925+min(ax1.YLim),...
            ' a ','edgecolor','k','fontsize',fontsize,'fontweight','bold','linewidth',linewidth-1,'backgroundcolor','w'); 
        cd .. % exit folder
% b) Hektoria & Green
cd('HekGreen'); % enter folder 
    files = dir('*.shp'); % load file
    % loop through files to grab centerline
    for j=1:length(files)
        if contains(files(j).name,'CL')
            file = shaperead(files(j).name);
            cl.lon = file.X; % Lon
            cl.lat = file.Y; % Lat
            [cl.X,cl.Y] = wgs2ps(cl.lon,cl.lat,'StandardParallel',-71,'StandardMeridian',0);
            % define x as distance along centerline
            x = zeros(1,length(cl.X));
            for k=2:length(x)
                x(k) = sqrt((cl.X(k)-cl.X(k-1))^2+(cl.Y(k)-cl.Y(k-1))^2)+x(k-1); %(m)
            end
        end
    end
    % initialize f
    clear f; f(1).lon = NaN; f(1).lat = NaN; f(1).date = NaN;
    % loop through files to grab terminus positions
    for j=1:length(files)
        if ~contains(files(j).name,'CL')
            file = shaperead(files(j).name);
            % loop through file
            for k=1:length(file)
                f(length(f)+1).lon = file(k).X; % Lon
                f(length(f)).lat = file(k).Y; % Lat
                % convert to polar stereographic coordinates
                [f(length(f)).X,f(length(f)).Y] = wgs2ps(f(length(f)).lon,f(length(f)).lat,'StandardParallel',-71,'StandardMeridian',0);
                f(length(f)).x = x(dsearchn([cl.X' cl.Y'],[nanmean(f(length(f)).X) nanmean(f(length(f)).Y)])); % x
            end
        end 
    end
    % plot
    col = parula(length(f)); % color scheme for plotting
    ax2=axes('position',[0.39 0.58 0.24 0.4]); hold on; grid on;      
        set(gca,'fontname',fontname,'fontsize',fontsize);
        xlabel('Easting (km)'); ylabel('Northing (km)');
        % plot centerline
        plot(cl.X(2:8)/10^3,cl.Y(2:8)/10^3,'k','linewidth',linewidth);
        annotation('textarrow',[0.633 0.635],[0.807 0.808],'headwidth',18,'headlength',16);
        % plot terminus positions
        for j=1:length(f)
            plot(f(j).X/10^3,f(j).Y/10^3,'color',col(j,:),'linewidth',linewidth);
        end
        axis tight; text((ax2.XLim(2)-ax2.XLim(1))*0.88+ax2.XLim(1),(max(ax2.YLim)-min(ax2.YLim))*0.925+min(ax2.YLim),...
            ' b ','edgecolor','k','fontsize',fontsize,'fontweight','bold','linewidth',linewidth-1,'backgroundcolor','w'); 
        cd .. % exit folder
% c) Jorum
cd('Jorum'); % enter folder 
    files = dir('*.shp'); % load file
    % loop through files to grab centerline
    for j=1:length(files)
        if contains(files(j).name,'CL')
            file = shaperead(files(j).name);
            cl.lon = file.X; % Lon
            cl.lat = file.Y; % Lat
            [cl.X,cl.Y] = wgs2ps(cl.lon,cl.lat,'StandardParallel',-71,'StandardMeridian',0);
            % define x as distance along centerline
            x = zeros(1,length(cl.X));
            for k=2:length(x)
                x(k) = sqrt((cl.X(k)-cl.X(k-1))^2+(cl.Y(k)-cl.Y(k-1))^2)+x(k-1); %(m)
            end
        end
    end
    % initialize f
    clear f; f(1).lon = NaN; f(1).lat = NaN; f(1).date = NaN;
    % loop through files to grab terminus positions
    for j=1:length(files)
        if ~contains(files(j).name,'CL')
            file = shaperead(files(j).name);
            % loop through file
            for k=1:length(file)
                f(length(f)+1).lon = file(k).X; % Lon
                f(length(f)).lat = file(k).Y; % Lat
                % convert to polar stereographic coordinates
                [f(length(f)).X,f(length(f)).Y] = wgs2ps(f(length(f)).lon,f(length(f)).lat,'StandardParallel',-71,'StandardMeridian',0);
                f(length(f)).x = x(dsearchn([cl.X' cl.Y'],[nanmean(f(length(f)).X) nanmean(f(length(f)).Y)])); % x
            end
        end 
    end
    % plot
    col = parula(length(f)); % color scheme for plotting
    ax3=axes('position',[0.72 0.58 0.23 0.4]); hold on; grid on;      
        set(gca,'fontname',fontname,'fontsize',fontsize);
        xlabel('Easting (km)'); ylabel('Northing (km)'); ylim([1277.5 1281.5])
        % plot centerline
        plot(cl.X(5:10)/10^3,cl.Y(5:10)/10^3,'k','linewidth',linewidth);
        annotation('arrow',[0.888 0.889],[0.983 0.985],'headwidth',18,'headlength',16);
        % plot terminus positions
        for j=1:length(f)
            plot(f(j).X/10^3,f(j).Y/10^3,'color',col(j,:),'linewidth',linewidth);
        end
        axis tight; text((ax3.XLim(2)-ax3.XLim(1))*0.88+ax3.XLim(1),(max(ax3.YLim)-min(ax3.YLim))*0.925+min(ax3.YLim),...
            ' c ','edgecolor','k','fontsize',fontsize,'fontweight','bold','linewidth',linewidth-1,'backgroundcolor','w'); 
        cd .. % exit folder
% regional map
ax4 = axes; hold on; imshow(LIMA);
    xlim([0.8056e3 1.4356e3]); ylim([4.146e3 4.846e3]);
    % plot borders around image
    plot([ax4.XLim(1) ax4.XLim(1)],ax4.YLim,'k','linewidth',linewidth-1);
    plot([ax4.XLim(2) ax4.XLim(2)],ax4.YLim,'k','linewidth',linewidth-1);    
    plot(ax4.XLim,[ax4.YLim(1) ax4.YLim(1)],'k','linewidth',linewidth-1);   
    plot(ax4.XLim,[ax4.YLim(2) ax4.YLim(2)],'k','linewidth',linewidth-1);   
    set(gca,'position',[0.11 0.08 0.25 0.4]);
    % plot location text labels
    text(869,4247,' a ','edgecolor','k','fontsize',fontsize-3,'fontweight','bold','linewidth',1,'backgroundcolor','w');
    text(890,4356,' b ','edgecolor','k','fontsize',fontsize-3,'fontweight','bold','linewidth',1,'backgroundcolor','w');
    text(960,4560,' c ','edgecolor','k','fontsize',fontsize-3,'fontweight','bold','linewidth',1,'backgroundcolor','w');    
    text(1030,4712,' d ','edgecolor','k','fontsize',fontsize-3,'fontweight','bold','linewidth',1,'backgroundcolor','w');    
    text(1170,4767,' e ','edgecolor','k','fontsize',fontsize-3,'fontweight','bold','linewidth',1,'backgroundcolor','w');    
    ax4.Position=[0.0 0.05 0.4 0.4];
    % add colorbar
    colormap = parula; c = colorbar('southoutside'); 
    c.FontName=fontname; c.FontSize = fontsize-1;
    c.Ticks = [0 0.5  1]; c.TickLabels = [{'2014'},{'2017'},{'2021'}];
% d) Crane
cd('Crane'); % enter folder 
    files = dir('*.shp'); % load file
    % loop through files to grab centerline
    for j=1:length(files)
        if contains(files(j).name,'CL')
            file = shaperead(files(j).name);
            cl.lon = file.X; % Lon
            cl.lat = file.Y; % Lat
            [cl.X,cl.Y] = wgs2ps(cl.lon,cl.lat,'StandardParallel',-71,'StandardMeridian',0);
            % define x as distance along centerline
            x = zeros(1,length(cl.X));
            for k=2:length(x)
                x(k) = sqrt((cl.X(k)-cl.X(k-1))^2+(cl.Y(k)-cl.Y(k-1))^2)+x(k-1); %(m)
            end
        end
    end
    % initialize f
    clear f; f(1).lon = NaN; f(1).lat = NaN; f(1).date = NaN;
    % loop through files to grab terminus positions
    for j=1:length(files)
        if ~contains(files(j).name,'CL')
            file = shaperead(files(j).name);
            % loop through file
            for k=1:length(file)
                f(length(f)+1).lon = file(k).X; % Lon
                f(length(f)).lat = file(k).Y; % Lat
                % convert to polar stereographic coordinates
                [f(length(f)).X,f(length(f)).Y] = wgs2ps(f(length(f)).lon,f(length(f)).lat,'StandardParallel',-71,'StandardMeridian',0);
                f(length(f)).x = x(dsearchn([cl.X' cl.Y'],[nanmean(f(length(f)).X) nanmean(f(length(f)).Y)])); % x
            end
        end 
    end
    % plot
    col = parula(length(f)); % color scheme for plotting
    ax5=axes('position',[0.39 0.08 0.23 0.4]); hold on; grid on;      
        set(gca,'fontname',fontname,'fontsize',fontsize);
        xlabel('Easting (km)'); ylabel('Northing (km)');
        % plot landsat image
        imagesc(LS.x/10^3,LS.y/10^3,flipud(LS.im*1.2)); colormap('gray'); 
        % plot centerline
        plot(cl.X(8:12)/10^3,cl.Y(8:12)/10^3,'k','linewidth',linewidth);
        %annotation('arrow',[0.515 0.52],[0.45 0.465],'headwidth',18,'headlength',16);
        %xlim([-2.4114e3  -2.4019e3]); ylim([1.264e3 1.2636e3]);
        % plot terminus positions
        for j=1:length(f)
            plot(f(j).X/10^3,f(j).Y/10^3,'color',col(j,:),'linewidth',linewidth);
        end
        axis tight; text((ax5.XLim(2)-ax5.XLim(1))*0.88+ax5.XLim(1),(max(ax5.YLim)-min(ax5.YLim))*0.925+min(ax5.YLim),...
            ' d ','edgecolor','k','fontsize',fontsize,'fontweight','bold','linewidth',linewidth-1,'backgroundcolor','w'); 
        cd .. % exit folder
% e) Flask
cd('Flask'); % enter folder 
    files = dir('*.shp'); % load file
    % loop through files to grab centerline
    for j=1:length(files)
        if contains(files(j).name,'CL')
            file = shaperead(files(j).name);
            cl.lon = file.X; % Lon
            cl.lat = file.Y; % Lat
            [cl.X,cl.Y] = wgs2ps(cl.lon,cl.lat,'StandardParallel',-71,'StandardMeridian',0);
            % define x as distance along centerline
            x = zeros(1,length(cl.X));
            for k=2:length(x)
                x(k) = sqrt((cl.X(k)-cl.X(k-1))^2+(cl.Y(k)-cl.Y(k-1))^2)+x(k-1); %(m)
            end
        end
    end
    % initialize f
    clear f; f(1).lon = NaN; f(1).lat = NaN; f(1).date = NaN;
    % loop through files to grab terminus positions
    for j=1:length(files)
        if ~contains(files(j).name,'CL')
            file = shaperead(files(j).name);
            % loop through file
            for k=1:length(file)
                f(length(f)+1).lon = file(k).X; % Lon
                f(length(f)).lat = file(k).Y; % Lat
                % convert to polar stereographic coordinates
                [f(length(f)).X,f(length(f)).Y] = wgs2ps(f(length(f)).lon,f(length(f)).lat,'StandardParallel',-71,'StandardMeridian',0);
                f(length(f)).x = x(dsearchn([cl.X' cl.Y'],[nanmean(f(length(f)).X) nanmean(f(length(f)).Y)])); % x
            end
        end 
    end
    % plot
    col = parula(length(f)); % color scheme for plotting
    ax6=axes('position',[0.72 0.08 0.23 0.4]); hold on; grid on;      
        set(gca,'fontname',fontname,'fontsize',fontsize);
        xlabel('Easting (km)'); ylabel('Northing (km)');
        % plot centerline
        plot(cl.X(9:12)/10^3,cl.Y(9:12)/10^3,'k','linewidth',linewidth);
        annotation('arrow',[0.804 0.803],[0.472 0.482],'headwidth',18,'headlength',16);
        % plot terminus positions
        for j=1:length(f)
            plot(f(j).X/10^3,f(j).Y/10^3,'color',col(j,:),'linewidth',linewidth);
        end
        axis tight; 
        text((ax6.XLim(2)-ax6.XLim(1))*0.88+ax6.XLim(1),(max(ax6.YLim)-min(ax6.YLim))*0.925+min(ax6.YLim),...
            ' e ','edgecolor','k','fontsize',fontsize,'fontweight','bold','linewidth',linewidth-1,'backgroundcolor','w'); 
        cd .. % exit folder

% save figure
if save_figure
    cd([homepath,'figures']);
    saveas(gcf,'regionalTermini.png','png');
    cd([homepath,'../write-ups/Thesis/figures/']);
    saveas(gcf,'regionalTermini.png','png');
    disp('figure 6 saved');
end



