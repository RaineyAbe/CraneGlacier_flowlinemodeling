%% Create Figures for MS Thesis
% RKA 2021

% Figures:
%   - Map of the study area
%   - Observed Conditions Time Series
%   - Cumulative Strain Rates / A adjustment
%   - Beta Solution
%   - Sensitivity Tests
%   - Regional glacier terminus time series

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
        
%% Map of the study area

close all;

save_figure = 1;    % = 1 to save figure
fontsize = 20;      % font size
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
    set(gcf,'position',[557    75   724   622],'color','w');
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
    ax1.Position = [0.135 0.1 0.73 0.85];
    ax2.Position=[0.001 0.1 1 0.85];     
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','62.00^oS','textcolor',[175 175 175]/255,...
        'HeadStyle','none','LineStyle', 'none','Position',[.97 .67 0 0],'FontSize',fontsize,'FontName','Arial');     
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','62.50^oS','textcolor',[175 175 175]/255, ...
        'HeadStyle','none','LineStyle', 'none','Position',[.97 .39 0 0],'FontSize',fontsize,'FontName','Arial');     
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','65.25^oW','textcolor',[175 175 175]/255, ...
        'HeadStyle','none','LineStyle', 'none','Position',[.64 .97 0 0],'FontSize',fontsize,'FontName','Arial');     
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','65.00^oW','textcolor',[175 175 175]/255, ...
        'HeadStyle','none','LineStyle', 'none','Position',[.36 .97 0 0],'FontSize',fontsize,'FontName','Arial'); 
    
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
    for i=1:length(labels)
        labels(i) = strcat(labels(i)," km");
    end
    text(cl.X(Idist)+1500,cl.Y(Idist),labels,'color',[253,224,221]/255,...
        'fontweight','bold','fontname','arial','fontsize',fontsize);

% add colorbar
    c = colorbar('Position',[0.29 0.58 0.02 0.15],'fontsize',fontsize,'fontweight',...
        'bold','color','w','fontname','arial');
    set(get(c,'Title'),'String','Speed (m a^{-1})','color','w','fontname',...
        'arial','Position',[5 100 0]);

% insert text labels
    % Crane
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Crane', ...
       'HeadStyle','none','LineStyle', 'none', 'TextRotation',70,'Position',[.52 0.65 0 0],...
       'FontSize',fontsize-1,'FontName','Arial','fontweight','bold','TextColor',[200 200 200]/255);
    % Former Larsen B Ice Shelf
    %txt = sprintf('Former Larsen \n    B Ice Shelf');
    %text(-2.405e6,1.282e6,txt,'color',[200 200 200]/255,'fontsize',12,'fontweight','bold');
    % Jorum
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Jorum', ...
        'HeadStyle','none','LineStyle', 'none', 'TextRotation',35,'Position',[.52 .855 0 0],...
        'FontSize',fontsize-1,'FontName','Arial','fontweight','bold','TextColor',[200 200 200]/255);    
    % Flask
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Flask', ...
       'HeadStyle','none','LineStyle', 'none', 'TextRotation',55,'Position',[.8 .19 0 0],...
       'FontSize',fontsize-1,'FontName','Arial','fontweight','bold','TextColor',[200 200 200]/255);    
    % Mapple
    %annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Mapple', ...
    %    'HeadStyle','none','LineStyle', 'none', 'TextRotation',66,'Position',[.65 .7 0 0],...
    %    'FontSize',fontsize-1,'FontName','Arial','fontweight','bold','TextColor',[200 200 200]/255);    
    % Melville
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Melville', ...
        'HeadStyle','none','LineStyle', 'none', 'TextRotation',70,'Position',[.68 .6 0 0],...
        'FontSize',fontsize-1,'FontName','Arial','fontweight','bold','TextColor',[200 200 200]/255);           
    % Pequod
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Pequod', ...
        'HeadStyle','none','LineStyle', 'none', 'TextRotation',66,'Position',[.74 .59 0 0],...
        'FontSize',fontsize-1,'FontName','Arial','fontweight','bold','TextColor',[200 200 200]/255);               
    % Starbuck
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Starbuck', ...
        'HeadStyle','none','LineStyle', 'none', 'TextRotation',38,'Position',[.78 .47 0 0],...
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
    L = imread('LIMA.jpg');
    imshow(L);
    hold on; plot(1100,4270,'o','color','y','markersize',10,'linewidth',3);
    
if save_figure
    set(gcf,'InvertHardCopy','off'); % save colors as is
    cd([homepath,'../figures']);
    saveas(gcf,'fig1.1_studyArea.png','png');
    cd([homepath,'../write-ups/Thesis/figures']);
    saveas(gcf,'fig1.1_studyArea.png','png');
    disp('figure 1 saved.');    
end

%% Observed Conditions Time Series

close all;

cd([homepath,'inputs-outputs']);

save_figure = 1; % = 1 to save figure
fontsize = 18; % fontsize for plots
font = 'Arial';
linewidth = 2; % line width for plots
markersize = 10; % marker size for plots
col1 = parula(length(2002:2019)); % color scheme for plotting

% 1. Ice surface elevation
h = load('Crane_surfaceElevationObs.mat').h;
col_h = [3 4 5 5 6 6 6 6 6 7 7 7 8 8 9 9 9 10.*ones(1,11) 11*ones(1,7) 12]+5;

% 2. Glacier bed (OIB)
hb = load('Crane_ObservedBed_Tate.mat').hb.hb0; % OIB
bathym = load('Crane_bathymetryData.mat').cl_trough; % Rebesco et al. (2014)

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
SMB = load('Crane_downscaledSMB_2009-2019.mat').SMB;
col_SMB = [8:18];

% Plot
figure(2); clf
set(gcf,'Position',[100 100 1000 600]);
ax2=axes('Position',[0.1 0.6 0.35 0.35],'linewidth',2,'fontsize',fontsize,'fontname',font); % geometry
    hold on; grid on; ylabel('Elevation (m)'); 
    xlim([0 60]); ylim([-1200 1200]);
    for i=1:length(h)
        if i==16
        else
            plot(ax2,cl.xi./10^3,h(i).surface,'linewidth',linewidth,...
                'color',col1(col_h(i),:),'HandleVisibility','off');
        end
    end
    plot(ax2,cl.xi./10^3,hb,'-k','linewidth',linewidth,'displayname','Bed');
    % Add text label
    text(55,(max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        '(a)','fontsize',fontsize,'linewidth',linewidth-1,'backgroundcolor','w','fontname',font);    
ax2=axes('Position',[0.55 0.6 0.35 0.35],'linewidth',2,'fontsize',fontsize,'fontname',font); % speed
    hold on; grid on; ylabel('Speed (m a^{-1})');
    xlim([0 60]);
    for i=1:length(U)
        plot(ax2,cl.xi./10^3,U(i).speed.*3.1536e7,'linewidth',linewidth,'color',col1(col_h(i),:));
    end
    % Add text label
    text(55,(max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        '(b)','fontsize',fontsize,'linewidth',linewidth-1,'backgroundcolor','w','fontname',font);      
ax3=axes('Position',[0.1 0.1 0.35 0.35],'linewidth',2,'fontsize',fontsize,'fontname',font); % terminus position
    hold on; grid on; xlabel('Distance Along Centerline (km)'); ylabel('Date');
    xlim([0 60]);
    for i=1:length(termx)
        plot(ax3,termx(i)./10^3,termdate(i),'ok','MarkerFaceColor',col1(col_term(i),:),'markersize',markersize);
    end
    % Add text label
    text(55,(max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        '(c)','fontsize',fontsize,'linewidth',linewidth-1,'backgroundcolor','w','fontname',font);     
ax4=axes('Position',[0.55 0.1 0.35 0.35],'linewidth',2,'fontsize',fontsize,'fontname',font); % speed
    hold on; grid on; xlabel('Distance Along Centerline (km)'); ylabel('Mean Annual SMB (m a^{-1})');
    xlim([0 60]); set(gca,'clim',[2002 2019]); 
    % Add colorbar
    colormap(col1); 
    c = colorbar('Limits',[2002 2019],'Ticks',[2002 2010 2019],'Position',[.92 .32 .03 .3410],...
        'fontsize',fontsize);
    for i=1:length(SMB)
        plot(ax4,cl.xi./10^3,SMB(i).smb_interp,'color',col1(col_SMB(i),:),'linewidth',linewidth);
    end
    % Add text label
    text(55,(max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        '(d)','fontsize',fontsize,'linewidth',linewidth-1,'backgroundcolor','w','fontname',font);          

% Save figure
if save_figure
    cd([homepath,'../write-ups/Thesis/figures']);
    saveas(gcf,'centerlineObservations.png','png');   
    disp('figure 2 saved.');
end

%% Cumulative Strain Rates / A adjustment

close all; 

save_figure = 1; 
fontsize = 18;
fontname = 'Arial';
linewidth = 2;

% Load A, A_adj, and eta_dot_cum
cd([homepath,'inputs-outputs']);
A = load('Crane_rateFactorA.mat').A; % rate factor
A_adj = load('Crane_adjustedRateFactor.mat').A_adj; % adjusted rate factor
eta_dot_cum = load('Crane_adjustedRateFactor.mat').eta_dot_cum; % annual cumulative strain rates

% plot
col1 = parula(length(eta_dot_cum(:,1))); % color scheme for plotting

figure(3); clf; 
set(gcf,'Position',[100 100 1000 400]);
ax1 = axes('position',[0.06 0.15 0.39 0.78]); 
    set(ax1,'linewidth',2,'fontsize',fontsize,'fontname',fontname);
    xlabel('Distance Along Centerline (km)'); xlim([0 45]);  
    %ylabel('$$ \Sigma ( \dot{\eta} ) (s^{-1})$$','Interpreter','latex','fontname','Arial','fontsize',18);
    ylabel('Cumulative Strain Rate');
    hold on; grid on;
    for i=2:length(eta_dot_cum(:,1))
        plot(ax1,cl.xi(1:135)/10^3,eta_dot_cum(i,1:135),'color',col1(i,:),'linewidth',linewidth-1); drawnow
    end
    plot(ax1,cl.xi(1:135)/10^3,nanmean(eta_dot_cum(:,1:135),1),'-k','linewidth',linewidth)
    text(55,(max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.08+min(get(gca,'YLim')),...
        '(a)','fontsize',fontsize,'linewidth',linewidth-1,'backgroundcolor','w','fontname',fontname);     
ax2 = axes('position',[0.53 0.15 0.39 0.78]);
    set(ax2,'linewidth',2,'fontsize',fontsize,'fontname',fontname);
    xlabel('Distance Along Centerline (km)'); ylabel('A (Pa^{-3} a^{-1})');
    hold on; grid on; xlim([0 45]); 
    for i=2:length(A_adj)
        plot(ax2,cl.xi(1:135)/10^3,A_adj(1:135)*3.1536e7,'-k','linewidth',linewidth,'displayname','A_{adj}');
    end
    plot(ax2,cl.xi(1:135)/10^3,A(1:135)*3.1536e7,'--k','linewidth',linewidth,'displayname','A');
    text(55,(max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.08+min(get(gca,'YLim')),...
        '(b)','fontsize',fontsize,'linewidth',linewidth-1,'backgroundcolor','w','fontname',fontname); 
cb = colorbar('position',[0.93 0.25 0.02 0.5],'fontname',fontname,...
    'fontsize',fontsize-2,'Ticks',0:0.5:1,'TickLabels',[{'2009'},{'2013'},{'2017'}]);

% save figure
if save_figure
    figure(3); 
    % save in thesis figures folder
    cd([homepath,'../write-ups/Thesis/figures/']);
    saveas(gcf,'cumulativeStrains.png','png');  
    disp('figure 3 saved.'); 
end
        
%% Beta Solution

close all;

save_figure = 1;                % = 1 to save figure
fontsize = 18;                   % font size
fontname = 'Arial';             % font name
linewidth = 2.5;                % line width
markersize = 25;                % marker size

% load variables
cd([homepath,'inputs-outputs/']);
beta = load('Crane_betaSolution.mat').beta;
xcf = load('Crane_betaSolution.mat').xcf;
U = load('Crane_betaSolution.mat').Un;
x0 = load('Crane_flowlineModelInitialization.mat').x0;
c0 = load('Crane_flowlineModelInitialization.mat').c0;
load('Crane_SCL_results.mat'); % stress-coupling length results
U_2018 = load('Crane_centerlineSpeedsWidthAveraged_2007-2018.mat').U_widthavg(20).speed;

% plot results
figure(4); clf
set(gcf,'Position',[272 143 1000 500]);
ax1=axes; set(gca,'position',[0.08 0.15 0.36 0.8]);
    set(gca,'fontsize',fontsize,'fontname',fontname,'linewidth',2);
    hold on; grid on; legend('Location','northwest'); xlim([0 50]);    
    xlabel('Distance Along Centerline (km)'); ylabel('\beta (s^{1/m} m^{-1/m})'); 
    plot(x0(1:c0)./10^3,beta(1:c0),'linewidth',2,'displayname','\beta _{sol}');     
    yyaxis right; ylabel('U (m a^{-1})'); ylim([0 750]);
    plot(cl.xi(1:135)./10^3,U_2018(1:135).*3.1536e7,'k','displayname','U_{obs}','linewidth',linewidth);
    plot(betax(1:dsearchn(betax',x0(c0)))./10^3,U(1:dsearchn(betax',x0(c0))).*3.1536e7,'--k','displayname','U_{mod}','linewidth',2); 
    text(42,(max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.05+min(get(gca,'YLim')),...
        '(a)','fontsize',fontsize,'linewidth',1,'backgroundcolor','w');   
ax2=axes; set(gca,'position',[0.62 0.15 0.36 0.8]);
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
        'fontsize',fontsize,'linewidth',1,'backgroundcolor','w','color','k','fontname',fontname);    
    text(42,(max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.05+min(get(gca,'YLim')),...
        '(b)','fontsize',fontsize-1,'linewidth',1,'backgroundcolor','w','fontname',fontname);    
    text(27.8,-200,'\Deltax','fontsize',13,'fontweight','bold','color','k');
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','}', ...
       'HeadStyle','none','LineStyle', 'none', 'TextRotation',90,'Position',[0.85 0.565 0 0],...
       'FontSize',24,'FontName',fontname,'TextColor','k','fontweight','bold');
   
if save_figure
    figure(4);
    % save in figures folder
    cd([homepath,'../figures/']);
    saveas(gcf,'betaSolution.png','png');  
    % save in thesis figures folder
    cd([homepath,'../write-ups/Thesis/figures/']);
    saveas(gcf,'betaSolution.png','png');  
    disp('figure 4 saved.');    
end
    
%% Sensitivity Tests

%close all;

save_figures = 1;    % = 1 to save figure
fontsize = 18;      % font size
fontname = 'Arial'; % font name
linewidth = 2;      % line width
markersize = 10;    % marker size

% load sensitivity test no change variables
cd([homepath,'inputs-outputs']);
load('Crane_2100_noChange.mat');
% load observations/parameters
x_cl = load('Crane_centerline.mat').x; y_cl = load('Crane_centerline.mat').y; 
    % define x as distance along centerline
    x0 = zeros(1,length(y_cl));
    for i=2:length(x0)
        x0(i) = sqrt((x_cl(i)-x_cl(i-1))^2+(y_cl(i)-y_cl(i-1))^2)+x0(i-1);
    end
load('Crane_flowlineModelInitialization.mat')
% define time stepping (s)
dt = 0.01*3.1536e7;
t_start = 0*3.1536e7;
t_end = 91*3.1536e7;
t = (t_start:dt:t_end);
clear dt t_start t_end

cd([homepath,'scripts/3_sensitivityTests/results/']);
addpath([homepath,'../matlabFunctions/cmocean_v2.0/cmocean']);

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
            files(length(files)+1).change = 0;
            files(length(files)).changeIn = 'fwd';
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
Ismb = flipud(find(strcmp('SMB',table2array(files(:,8)))));
Ismr = find(strcmp('SMR',table2array(files(:,8))));
Ifwd = find(strcmp('fwd',table2array(files(:,8))));
files = table2struct(files);

% define color schemes
col1=cmocean('thermal',length(Ismb)+3); % color scheme for figures 5 & 6 
col1(1,:) = []; col1(end-1:end,:)=[]; % don't use end member colors
col2=[216/255,179/255,101/255; 0.605513566501496,0.603942186152517,0.600896223533862; 90/255,180/255,172/255]; % color scheme for figure 7

% set up axes
loop=1; % loop to minimize plotting commands
while loop==1
    figure(5); clf 
    set(gcf,'Position',[100 200 1000 1000]);
    % SMR
    ax1=axes('position',[0.08 0.73 0.225 0.25]); hold on; grid on; % a) geometry
        set(gca,'fontsize',fontsize,'linewidth',2,'YTick',-1500:500:1500);
        ylabel('Elevation (m)'); 
        xlim([40 70]); ylim([-800 200]);
        plot(x1./10^3,hb1,'-k','linewidth',2);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.89+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.91+min(get(gca,'YLim')),...
            '(a)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1);  
    ax2=axes('position',[0.385 0.73 0.225 0.25]); hold on; grid on; % b) speed
        set(gca,'fontsize',fontsize,'linewidth',2,'YTick',0:200:1000);        
        ylabel('U (m a^{-1})');  
        xlim([40 70]); ylim([300 700]);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.89+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.91+min(get(gca,'YLim')),...
            '(b)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1); 
    ax3=axes('position',[0.69 0.73 0.225 0.25]); hold on; grid on; % c) gl/cf
        set(gca,'fontsize',fontsize,'linewidth',2);
        ylabel('SMR_{max} (m a^{-1})');
        xlim([40 70]); ylim([-smr0*3.1536e7-1 -smr0*3.1536e7+11]);               
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.89+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.91+min(get(gca,'YLim')),...
            '(c)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1);     
        % add colorbar
        colormap(col1);
        cb1=colorbar('position',[0.935 0.76 0.012 0.2],'fontname',fontname,'fontsize',fontsize-3,'Ticks',0:1/5:1,'TickLabels',string(0:2:10));
        set(get(cb1,'label'),'String','\DeltaSMR (m a^{-1})','fontsize',fontsize-3);
    % SMB
    ax4=axes('position',[0.08 0.41 0.225 0.25]); hold on; grid on; % d) geometry
        set(gca,'fontsize',fontsize,'linewidth',2,'YTick',-1500:500:1500);
        ylabel('Elevation (m)'); 
        xlim([40 70]); ylim([-800 200]);
        plot(x1./10^3,hb1,'-k','linewidth',2);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.89+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.91+min(get(gca,'YLim')),...
            '(d)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1);  
    ax5=axes('position',[0.385 0.41 0.225 0.25]); hold on; grid on; % e) speed
        set(gca,'fontsize',fontsize,'linewidth',2,'YTick',0:200:1000);        
        ylabel('U (m a^{-1})');  
        xlim([40 70]); ylim([300 700]);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.89+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.91+min(get(gca,'YLim')),...
            '(e)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1); 
    ax6=axes('position',[0.69 0.41 0.225 0.25]); hold on; grid on; % f) gl/cf
        set(gca,'fontsize',fontsize,'linewidth',2);
        ylabel('SMB_{mean} (m a^{-1})');
        xlim([40 70]); ylim([nanmean(smb0)*3.1536e7-11 nanmean(smb0)*3.1536e7+1]);              
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.89+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.91+min(get(gca,'YLim')),...
            '(f)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1);  
        % add colorbar
        colormap(col1);
        cb2=colorbar('position',[0.935 0.44 0.012 0.2],'fontname',fontname,'fontsize',fontsize-3,'Ticks',0:1/5:1,'TickLabels',string(0:-2:-10));
        set(get(cb2,'label'),'String','\DeltaSMB (m a^{-1})','fontsize',fontsize-3);
    % fwd    
    ax7=axes('position',[0.08 0.09 0.225 0.25]); hold on; grid on; % g) geometry
        set(gca,'fontsize',fontsize,'linewidth',2,'YTick',-1500:500:1500);
        xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)'); 
        xlim([40 70]); ylim([-800 200]);
        plot(x1./10^3,hb1,'-k','linewidth',2);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.89+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.91+min(get(gca,'YLim')),...
            '(g)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1);  
    ax8=axes('position',[0.385 0.09 0.225 0.25]); hold on; grid on; % h) speed
        set(gca,'fontsize',fontsize,'linewidth',2,'YTick',0:200:1000);        
        xlabel('Distance Along Centerline (km)'); ylabel('U (m a^{-1})');  
        xlim([40 70]); ylim([300 700]);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.89+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.91+min(get(gca,'YLim')),...
            '(h)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1); 
    ax9=axes('position',[0.69 0.09 0.225 0.25]); hold on; grid on; % i) gl/cf
        set(gca,'fontsize',fontsize,'linewidth',2);
        xlabel('Distance Along Centerline (km)'); ylabel('FWD (m)');
        xlim([40 70]); ylim([fwd0-0.5 fwd0+5.5]);         
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.89+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.91+min(get(gca,'YLim')),...
            '(i)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1); 
        % add colorbar
        colormap(col1);
        cb3=colorbar('position',[0.935 0.12 0.012 0.2],'fontname',fontname,'fontsize',fontsize-3,'Ticks',0:1/2:1,'TickLabels',string(0:0.5:5)); 
        set(get(cb3,'label'),'String','\DeltaFWD (m)','fontsize',fontsize-3);    
    % ice mass discharge across grounding line
    figure(6); clf;
    set(gcf,'position',[50 300 1000 400]); 
    ax10 = axes('position',[0.08 0.18 0.27 0.75]); hold on; 
        set(gca,'fontsize',fontsize,'linewidth',2); grid on; 
        xlabel('SMR_{max} (m a^{-1})'); ylabel('F_{gl} (Gt a^{-1})'); 
        xlim([-smr0*3.1536e7-1 -smr0*3.1536e7+11]); ylim([0.5 1.1]); 
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.08+min(get(gca,'YLim')),...
            '(a)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1); 
        % add colorbar
        colormap(col1); 
        cb4=colorbar('northoutside','fontname',fontname,'fontsize',fontsize-2,'Ticks',0:1/5:1,'TickLabels',string(0:2:10));
        set(get(cb4,'label'),'String','\DeltaSMR (m a^{-1})','fontsize',fontsize-2);
    ax11 = axes('position',[0.39 0.18 0.27 0.75]); hold on;
        set(gca,'fontsize',fontsize,'linewidth',2); grid on;
        xlabel('SMB_{mean} (m a^{-1})');
        xlim([nanmean(smb0(1:135)*3.1536e7-11) nanmean(smb0(1:135)*3.1536e7+1)]); ylim([0.5 1.1]);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.08+min(get(gca,'YLim')),...
            '(b)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1); 
        % add colorbar
        colormap(col1);
        cb5=colorbar('northoutside','fontname',fontname,'fontsize',fontsize-2,'Ticks',0:1/5:1,'TickLabels',string(0:-2:-10));
        set(get(cb5,'label'),'String','\DeltaSMB (m a^{-1})','fontsize',fontsize-2);
    ax12 = axes('position',[0.7 0.18 0.27 0.75]); hold on;
        set(gca,'fontsize',fontsize,'linewidth',2); grid on;
        xlabel('FWD (m)'); 
        xlim([fwd0-0.5 fwd0+5.5]); ylim([0.5 1.1]);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.08+min(get(gca,'YLim')),...
            '(c)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1); 
         % add colorbar
        colormap(col1);
        cb6=colorbar('northoutside','fontname',fontname,'fontsize',fontsize-2,'Ticks',0:1/2:1,'TickLabels',string(0:0.5:5)); 
        set(get(cb6,'label'),'String','\DeltaFWD (m)','fontsize',fontsize-2);
    figure(7); clf;  
    set(gcf,'Position',[200 200 1000 1000]);   
    ax13 = axes('position',[0.08 0.09 0.85 0.85]);
        set(gca,'fontsize',fontsize,'linewidth',2); grid on;
        xlabel('Year'); ylabel('Grounding Line Discharge (Gt a^{-1})');
        ylim([0.4 2.6]);
        hold on; grid on; legend('Position',[0.13 0.75 0.18 0.125]);
        xlim([2008 2101]);
        yyaxis right; ylabel('Calving Front Position (km)');
    loop=loop+1; % exit loop        
end

% SMR
dsmr = zeros(length(Ismr),1); % change in SMR
dHgl = zeros(length(Ismr),1); % change in thickness at the grounding line
dUgl = zeros(length(Ismr),1); % change in ice speed at the grounding line
dL = zeros(length(Ismr),1); % change in length
dgl = zeros(length(Ismr),1); % change in grounding line position
Ugl = zeros(length(Ismr),1); % speed at grounding line 
Fsmr = zeros(length(Ismr),1); % grounding line discharge
for i=1:length(Ismr)
    load(files(Ismr(i)).name);
    % calculate grounding line discharge
    W = interp1(x0,W0,x2); % interpolate width on spatial grid
    % F (Gt/a) = (U m/s)*(A m^2)*(917 kg/m^3)*(3.1536e7 s/a)*(1e-12 Gt/kg)
    % use the mean of a rectangular and an ellipsoidal bed
    Fsmr(i) = (H2(gl2)*U2(gl2)*W(gl2)*pi*1/4)*917*1e-12*3.1536e7; % Gt/a
    figure(5); hold on;
        % ice surface
        plot(ax1,x2(1:c2)./10^3,h2(1:c2),'color',col1(i,:),'linewidth',linewidth);
        % calving front
        plot(ax1,[x2(c2) x2(c2)]./10^3,[h2(c2) h2(c2)-H2(c2)],'color',col1(i,:),'linewidth',linewidth);
        % floating bed
        plot(ax1,x2(gl2:c2)./10^3,h2(gl2:c2)-H2(gl2:c2),'color',col1(i,:),'linewidth',linewidth);
        % ice surface speed
        plot(ax2,x2(1:c2)./10^3,U2(1:c2).*3.1536e7,'color',col1(i,:),'linewidth',linewidth);
        % dH
        dHgl(i) = H2(gl2)-H1(gl1); 
        %figure(10); hold on; plot(H2,'color',col(i,:),'linewidth',2);
        % cf and gl positions
        plot(ax3,x2(gl2)/1e3,files(Ismr(i)).change-smr0*3.1536e7,'x',...
            'markersize',markersize,'linewidth',linewidth,'color',col1(i,:));
        plot(ax3,x2(c2)/1e3,files(Ismr(i)).change-smr0*3.1536e7,'o',...
            'markersize',markersize,'linewidth',linewidth,'color',col1(i,:));        
    % ice mass discharge
    figure(6);
    plot(ax10,files(Ismr(i)).change-smr0*3.1536e7,Fsmr(i),'o','markersize',markersize,'linewidth',linewidth,'color',col1(i,:));
    figure(7);
    if i==1
        yyaxis left; plot(ax13,t/3.1536e7+2009,movmean(Fgl2,100),'-k','linewidth',linewidth,'displayname','no change');  
        yyaxis right; plot(ax13,t/3.1536e7+2009,movmean(XCF2/10^3,100),'--k','linewidth',linewidth,'HandleVisibility','off');
    elseif i==length(Ismr)
        yyaxis left; plot(ax13,t/3.1536e7+2009,movmean(Fgl2,100),'-','color',col2(1,:),'linewidth',linewidth,'displayname','+1^oC ocean temp.');
        yyaxis right; plot(ax13,t/3.1536e7+2009,movmean(XCF2/10^3,100),'--','color',col2(1,:),'linewidth',linewidth,'HandleVisibility','off');        
    end
    % save results 
    dsmr(i) = files(Ismr(i)).change; % m/a
    dUgl(i) = (U2(gl2)-U1(gl1)).*3.1536e7; % m/a     
    dL(i) = x2(c2)-x1(c1); % m   
    dgl(i) = x2(gl2)-x1(gl1); % m
    Ugl(i) = U2(gl2)*3.1536e7; % m/a
end
% create table to store result quantities
varNames = {'dsmr','dL','dxgl','dHgl','dUgl','Fsmr'};
T_smr = table(round(dsmr),round(dL),round(dgl),round(dHgl),round(dUgl),Fsmr,'VariableNames',varNames);

% SMB
dsmb = zeros(length(Ismb),1); % change in SMR
dHgl = zeros(length(Ismb),1); % change in thickness at the grounding line
dUgl = zeros(length(Ismb),1); % change in ice speed at the grounding line
dL = zeros(length(Ismb),1); % change in length
dgl = zeros(length(Ismb),1); % change in grounding line position
Ugl = zeros(length(Ismb),1); % speed at grounding line 
Fsmb = zeros(length(Ismb),1); % grounding line discharge
for i=1:length(Ismb)
    load(files(Ismb(i)).name);
    % calculate grounding line discharge
    W = interp1(x0,W0,x2); % interpolate width on spatial grid
    % F (Gt/a) = (U m/s)*(A m^2)*(917 kg/m^3)*(3.1536e7 s/a)*(1e-12 Gt/kg)
    % use the mean of a rectangular and an ellipsoidal bed
    Fsmb(i) = (H2(gl2)*U2(gl2)*W(gl2)*pi*1/4)*917*1e-12*3.1536e7; % Gt/a
    figure(5); hold on;
        % ice surface
        plot(ax4,x2(1:c2)./10^3,h2(1:c2),'color',col1(i,:),'linewidth',linewidth);
        % calving front
        plot(ax4,[x2(c2) x2(c2)]./10^3,[h2(c2) h2(c2)-H2(c2)],'color',col1(i,:),'linewidth',linewidth);
        % floating bed
        plot(ax4,x2(gl2:c2)./10^3,h2(gl2:c2)-H2(gl2:c2),'color',col1(i,:),'linewidth',linewidth);
        % ice surface speed
        plot(ax5,x2(1:c2)./10^3,U2(1:c2).*3.1536e7,'color',col1(i,:),'linewidth',linewidth);
        % dH
        dHgl(i) = H2(gl2)-H1(gl1); 
        % cf and gl positions
        plot(ax6,x2(gl2)/1e3,nanmean(files(Ismb(i)).change+smb0(1:135)*3.1536e7),'x',...
            'markersize',10,'linewidth',linewidth,'color',col1(i,:));
        plot(ax6,x2(c2)/1e3,nanmean(files(Ismb(i)).change+smb0(1:135)*3.1536e7),'o',...
            'markersize',10,'linewidth',linewidth,'color',col1(i,:));        
    % ice mass discharge
    figure(6);
    plot(ax11,nanmean(files(Ismb(i)).change+smb0(1:135)*3.1536e7),Fsmb(i),'o','markersize',markersize,'linewidth',linewidth,'color',col1(i,:));
    figure(7);
    if i==length(Ismb)
        yyaxis left; plot(ax13,t/3.1536e7+2009,movmean(Fgl2,100),'-','color',col2(2,:),'linewidth',linewidth,'displayname','+1^oC air temp.');
        yyaxis right; plot(ax13,t/3.1536e7+2009,movmean(XCF2/10^3,100),'--','color',col2(2,:),'linewidth',linewidth,'HandleVisibility','off');
    end
    % save results 
    dsmb(i) = files(Ismb(i)).change; % m/a
    dUgl(i) = (U2(gl2)-U1(gl1)).*3.1536e7; % m/a     
    dL(i) = x2(c2)-x1(c1); % m   
    dgl(i) = x2(gl2)-x1(gl1); % m
    Ugl(i) = U2(gl2)*3.1536e7; % m/a
end
% create table to store result quantities
varNames = {'dsmb','dL','dxgl','dHgl','dUgl','Fsmb'};
T_smb = table(round(dsmb),round(dL),round(dgl),round(dHgl),round(dUgl),Fsmb,'VariableNames',varNames);

% fwd
dfwd = zeros(length(Ifwd),1); % change in fwd
dHgl = zeros(length(Ifwd),1); % change in mean thickness
dUgl = zeros(length(Ifwd),1); % change in mean ice speed
dL = zeros(length(Ifwd),1); % change in length
dgl = zeros(length(Ifwd),1); % change in grounding line position
Ugl = zeros(length(Ifwd),1); % speed at grounding line 
Ffwd = zeros(length(Ifwd),1); % grounding line discharge
col1 = cmocean('thermal',length(Ifwd)+3); col1(1,:)=[];col1(end-1:end,:)=[];
for i=1:length(Ifwd)
    % load file
    load(files(Ifwd(i)).name);
    % calculate grounding line discharge
    W = interp1(x0,W0,x2); % interpolate width on spatial grid
    % F (Gt/a) = (U m/s)*(A m^2)*(917 kg/m^3)*(3.1536e7 s/a)*(1e-12 Gt/kg)
    % use the mean of a rectangular and an ellipsoidal bed
    Ffwd(i) = (H2(gl2)*U2(gl2)*W(gl2))*pi*1/4*917*1e-12*3.1536e7; % Gt/a    
    % plot
    figure(5); hold on;
        % ice surface
        plot(ax7,x2(1:c2)./10^3,h2(1:c2),'color',col1(i,:),'linewidth',linewidth);
        % calving front
        plot(ax7,[x2(c2) x2(c2)]./10^3,[h2(c2)-H2(c2) h2(c2)],'color',col1(i,:),'linewidth',linewidth);
        % floating bed
        plot(ax7,x2(gl2:c2)./10^3,h2(gl2:c2)-H2(gl2:c2),'color',col1(i,:),'linewidth',linewidth);
        % ice surface speed
        plot(ax8,x2(1:c2)./10^3,U2(1:c2).*3.1536e7,'color',col1(i,:),'linewidth',linewidth);
        % dH
        dHgl(i) = H2(gl2)-H1(gl1); % m
        % cf and gl positions
        plot(ax9,x2(gl2)/10^3,files(Ifwd(i)).change+fwd0,'x',...
            'markersize',10,'linewidth',linewidth,'color',col1(i,:));
        plot(ax9,x2(c2)/10^3,files(Ifwd(i)).change+fwd0,'o',...
            'markersize',10,'linewidth',linewidth,'color',col1(i,:));        
    % ice mass discharge
    figure(6); 
    plot(ax12,files(Ifwd(i)).change+fwd0,Ffwd(i),'o','markersize',markersize,'linewidth',linewidth,'color',col1(i,:));    
    figure(7);
    if i==length(Ifwd)
        yyaxis left; plot(ax13,t/3.1536e7+2009,movmean(Fgl2,100),'-','color',col2(3,:),'linewidth',linewidth,'displayname','+5 m FWD');
        yyaxis right; plot(ax13,t/3.1536e7+2009,movmean(XCF2/10^3,100),'--','color',col2(3,:),'linewidth',linewidth,'HandleVisibility','off');
    end
    % save results 
    dfwd(i) = files(Ifwd(i)).change; % m/a
    dUgl(i) = (U2(gl2)-U1(gl1)).*3.1536e7; % m/a     
    dL(i) = x2(c2)-x1(c1); % m   
    dgl(i) = x2(gl2)-x1(gl1); % m
    Ugl(i) = U2(gl2)*3.1536e7; % m/a
end
% create table to store result quantities
varNames = {'dfwd','dL','dxgl','dHgl','dUgl','Ffwd'};
T_fwd = table(dfwd,round(dL),round(dgl),round(dHgl),round(dUgl),Ffwd,'VariableNames',varNames);

% save figures 
if save_figures
    % save in thesis figures folders
    cd([homepath,'../write-ups/Thesis/figures/']);
    figure(5);    
    saveas(gcf,'sensitivityTests_geom+speed.png','png');
    figure(6);
    saveas(gcf,'sensitivityTests_discharge.png','png');
    figure(7);
    saveas(gcf,'sensitivityTests_FglXcf.png','png');
    disp('figures 5-7 saved.');
end

%% Regional Glacier Terminus Time Series

close all;

cd([homepath,'data/terminus/regional/']);

save_figure = 0; % = 1 to save figure
fontsize = 18; 
fontname = 'Arial';
linewidth = 3; 

% add path to required functions
addpath([homepath,'matlabFunctions/']);
addpath([homepath,'../matlabFunctions/line2arrow-kakearney-pkg-8aead6f/']);
addpath([homepath,'../matlabFunctions/line2arrow-kakearney-pkg-8aead6f/line2arrow/']);
addpath([homepath,'..//matlabFunctions/line2arrow-kakearney-pkg-8aead6f/axescoord2figurecoord/']);
addpath([homepath,'../matlabFunctions/line2arrow-kakearney-pkg-8aead6f/parsepv']);

% load Landsat images
cd([homepath,'data/Imagery/']);
% LSA 
landsatA = dir('*217105*B8.TIF');
    [LSA.im,LSA.R] = readgeoraster(landsatA.name); [LSA.ny,LSA.nx] = size(LSA.im);
    % polar stereographic coordinates of image boundaries
    LSA.x = linspace(min(LSA.R.XWorldLimits),max(LSA.R.XWorldLimits),LSA.nx); 
    LSA.y = linspace(min(LSA.R.YWorldLimits),max(LSA.R.YWorldLimits),LSA.ny);
% LSB
landsatB = dir('*217106*B8.TIF');
    [LSB.im,LSB.R] = readgeoraster(landsatB.name); [LSB.ny,LSB.nx] = size(LSB.im);
    % polar stereographic coordinates of image boundaries
    LSB.x = linspace(min(LSB.R.XWorldLimits),max(LSB.R.XWorldLimits),LSB.nx); 
    LSB.y = linspace(min(LSB.R.YWorldLimits),max(LSB.R.YWorldLimits),LSB.ny);

% load LIMA .jpg
cd([homepath,'data/Imagery/']);
L = imread('LIMA.jpg');

% cd to terminus coordinate shapefiles
cd([homepath,'data/terminus/regional/']);

% set up figure, subplots, and plot
figure(7); clf; 
set(gcf,'Position',[100 200 1000 1000]); 
% regional map
ax1 = axes('position',[0.08 0.56 0.23 0.4]); hold on; imshow(L);
    ax1.XLim=[0.7139e3 1.4667e3]; ax1.YLim=[3.8765e3 4.6369e3];
    % plot borders around image
    plot([ax1.XLim(1) ax1.XLim(1)],ax1.YLim,'k','linewidth',linewidth-1);
    plot([ax1.XLim(2) ax1.XLim(2)],ax1.YLim,'k','linewidth',linewidth-1);    
    plot(ax1.XLim,[ax1.YLim(1) ax1.YLim(1)],'k','linewidth',linewidth-1);   
    plot(ax1.XLim,[ax1.YLim(2) ax1.YLim(2)],'k','linewidth',linewidth-1);   
    % plot location text labels
    text(940,3950,'(a)','edgecolor','k','fontsize',fontsize-3,'linewidth',1,'backgroundcolor','w');
    text(1000,4120,'(b)','edgecolor','k','fontsize',fontsize-3,'linewidth',1,'backgroundcolor','w');
    text(1013,4371,'(c)','edgecolor','k','fontsize',fontsize-3,'linewidth',1,'backgroundcolor','w');    
    text(1080,4460,'(d)','edgecolor','k','fontsize',fontsize-3,'linewidth',1,'backgroundcolor','w');    
    text(1085,4550,'(e)','edgecolor','k','fontsize',fontsize-3,'linewidth',1,'backgroundcolor','w');    
    % add colorbar
    colormap(ax1,'parula'); c = colorbar('southoutside'); 
    c.FontName=fontname; c.FontSize = fontsize-1;
    c.Ticks = [0 0.5  1]; c.TickLabels = [{'2014'},{'2017'},{'2021'}];
% a) Edgeworth
cd([homepath,'data/terminus/regional/Edgeworth/']); % enter folder 
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
    col1 = parula(length(f)); % color scheme for plotting
    ax2=axes('position',[0.39 0.58 0.23 0.4]); hold on; grid on;      
        set(gca,'fontname',fontname,'fontsize',fontsize);
        xlabel('Easting (km)'); ylabel('Northing (km)');
        % plot LIMA
        colormap(ax2,'gray');
        imagesc(LSA.x/10^3,LSA.y/10^3,flipud(LSA.im)); 
        ax2.XLim=[-2.4508e3 -2.443e3]; ax2.YLim=[1.4138e3 1.4204e3]; 
        % plot centerline
        l1=line(cl.X(3:8)/10^3,cl.Y(3:8)/10^3,'color','k','linewidth',linewidth);
        % arrow
        line2arrow(l1,'color','k','linewidth',linewidth,'headwidth',20,'headlength',20);
        % plot terminus positions
        for j=1:length(f)
            plot(f(j).X/10^3,f(j).Y/10^3,'color',col1(j,:),'linewidth',linewidth);
        end
        % plot label
        text((ax2.XLim(2)-ax2.XLim(1))*0.88+ax2.XLim(1),(max(ax2.YLim)-min(ax2.YLim))*0.925+min(ax2.YLim),...
            '(a)','edgecolor','k','fontsize',fontsize,'linewidth',1,'backgroundcolor','w'); 
% b) Drygalski
cd([homepath,'data/terminus/regional/Drygalski/']); % enter folder 
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
    col1 = parula(length(f)); % color scheme for plotting
    ax3=axes('position',[0.72 0.58 0.24 0.4]); hold on; grid on;      
        set(gca,'fontname',fontname,'fontsize',fontsize);
        xlabel('Easting (km)'); ylabel('Northing (km)');
        % plot LIMA
        colormap(ax3,'gray');
        imagesc(LSA.x/10^3,LSA.y/10^3,flipud(LSA.im)); 
        ax3.XLim = [-2.4416e3 -2.4247e3]; ax3.YLim = [1.3552e3 1.3721e3];        
        % plot centerline
        l3=line(cl.X(4:7)/10^3,cl.Y(4:7)/10^3,'color','k','linewidth',linewidth);
        % arrow
        line2arrow(l3,'color','k','linewidth',linewidth,'headwidth',20,'headlength',20);        
        % plot terminus positions
        for j=1:length(f)
            plot(f(j).X/10^3,f(j).Y/10^3,'color',col1(j,:),'linewidth',linewidth);
        end
        % plot label
        text((ax3.XLim(2)-ax3.XLim(1))*0.88+ax3.XLim(1),(max(ax3.YLim)-min(ax3.YLim))*0.925+min(ax3.YLim),...
            '(b)','edgecolor','k','fontsize',fontsize,'linewidth',1,'backgroundcolor','w'); 
% c) Hektoria & Green
cd([homepath,'data/terminus/regional/HekGreen/']); % enter folder 
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
    col1 = parula(length(f)); % color scheme for plotting
    ax4=axes('position',[0.08 0.08 0.23 0.4]); hold on; grid on;      
        set(gca,'fontname',fontname,'fontsize',fontsize);
        xlabel('Easting (km)'); ylabel('Northing (km)');
        % plot landsat image
        colormap(ax4,'gray');
        imagesc(LSB.x/10^3,LSB.y/10^3,flipud(LSB.im)); 
        ax4.XLim = [-2.4364e3 -2.4164e3]; ax4.YLim = [1.3002e3 1.3201e3];        
        % plot centerline
        l4=line(cl.X(2:8)/10^3,cl.Y(2:8)/10^3,'color','k','linewidth',linewidth);
        % arrow
        line2arrow(l4,'color','k','linewidth',linewidth,'headwidth',20,'headlength',20);        
        % plot terminus positions
        for j=1:length(f)
            plot(f(j).X/10^3,f(j).Y/10^3,'color',col1(j,:),'linewidth',linewidth);
        end
        text((ax4.XLim(2)-ax4.XLim(1))*0.88+ax4.XLim(1),(max(ax4.YLim)-min(ax4.YLim))*0.925+min(ax4.YLim),...
            '(c)','edgecolor','k','fontsize',fontsize,'linewidth',1,'backgroundcolor','w'); 
% d) Jorum
cd([homepath,'data/terminus/regional/Jorum/']); % enter folder 
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
    col1 = parula(length(f)); % color scheme for plotting
    ax5=axes('position',[0.39 0.08 0.23 0.4]); hold on; grid on;      
        set(gca,'fontname',fontname,'fontsize',fontsize);
        xlabel('Easting (km)'); ylabel('Northing (km)');
        % plot landsat image
        colormap(ax5,'gray');
        imagesc(LSB.x/10^3,LSB.y/10^3,flipud(LSB.im)); 
        ax5.XLim = [-2.4193e3 -2.4134e3]; ax5.YLim = [1.2762e3 1.2821e3];                
        % plot centerline
        l5=line(cl.X(5:10)/10^3,cl.Y(5:10)/10^3,'color','k','linewidth',linewidth);
        % arrow
        line2arrow(l5,'color','k','linewidth',linewidth,'headwidth',20,'headlength',20);        
        % plot terminus positions
        for j=1:length(f)
            plot(f(j).X/10^3,f(j).Y/10^3,'color',col1(j,:),'linewidth',linewidth);
        end
        text((ax5.XLim(2)-ax5.XLim(1))*0.88+ax5.XLim(1),(max(ax5.YLim)-min(ax5.YLim))*0.925+min(ax5.YLim),...
            '(d)','edgecolor','k','fontsize',fontsize,'linewidth',1,'backgroundcolor','w'); 
% e) Crane
cd([homepath,'data/terminus/regional/Crane/']); % enter folder 
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
    col1 = parula(length(f)); % color scheme for plotting
    ax6=axes('position',[0.72 0.08 0.23 0.4]); hold on; grid on;      
        set(gca,'fontname',fontname,'fontsize',fontsize);
        xlabel('Easting (km)'); ylabel('Northing (km)');
        % plot landsat image
        colormap(ax6,'gray');
        imagesc(LSB.x/10^3,LSB.y/10^3,flipud(LSB.im)); 
        ax6.XLim=[-2.4132e3 -2.4005e3]; ax6.YLim=[1.2635e3 1.2762e3];                
        % plot centerline
        l6=line(cl.X(6:13)/10^3,cl.Y(6:13)/10^3,'color','k','linewidth',linewidth);
        % arrow
        line2arrow(l6,'color','k','linewidth',linewidth,'headwidth',20,'headlength',20);        
        % plot terminus positions
        for j=1:length(f)
            plot(f(j).X/10^3,f(j).Y/10^3,'color',col1(j,:),'linewidth',linewidth);
        end
        text((ax6.XLim(2)-ax6.XLim(1))*0.88+ax6.XLim(1),(max(ax6.YLim)-min(ax6.YLim))*0.925+min(ax6.YLim),...
            '(e)','edgecolor','k','fontsize',fontsize,'linewidth',1,'backgroundcolor','w'); 

% save figure
if save_figure
    cd([homepath,'../write-ups/Thesis/figures/']);
    saveas(gcf,'regionalTermini.png','png');
    disp('figure 6 saved');
end



