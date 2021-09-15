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
homepath = '/Users/raineyaberle/Desktop/Research/CraneModeling/CraneGlacier_flowlinemodeling/';
addpath([homepath,'matlabFunctions/']);
addpath([homepath,'matlabFunctions/gridLegend_v1.4/']);
addpath([homepath,'matlabFunctions/cmocean_v2.0/cmocean/']);

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
fontsize = 18;      % font size
fontname = 'Arial'; % font name
markersize = 8;    % marker size
linewidth = 2;      % line width

% Add paths to functions
    addpath([homepath,'../matlabFunctions']);
    addpath([homepath,'../matlabFunctions/cmocean_v2.0/cmocean']);
    addpath([homepath,'../matlabFunctions/AntarcticMappingTools_v5.17/AntarcticMappingTools']);

% Specify xticks/yticks coordinates and axes limits
    xticks = linspace(-2.455e6,-2.375e6,5); yticks=linspace(1.21e6,1.29e6,5);
    xlimits=[xticks(1) xticks(end)]; ylimits = [yticks(1) yticks(end)];
    xlimits2 = -2.4260e6:0.01e6:-2.3960e6; ylimits2 = 1.2353e6:0.01e6:1.2653e6; 
    
% load & display Landsat image 1
    cd([homepath,'../Imagery/LC08_L1GT_218106_20191013_20191018_01_T2']);
    landsat = dir('*B8.TIF');
    [LS.im,LS.R] = readgeoraster(landsat.name); [LS.ny,LS.nx] = size(LS.im);
    % polar stereographic coordinates of image boundaries
    LS.x = linspace(min(LS.R.XWorldLimits),max(LS.R.XWorldLimits),LS.nx); 
    LS.y = linspace(min(LS.R.YWorldLimits),max(LS.R.YWorldLimits),LS.ny);
    % display image on ax(1)
    F1 = figure(1); clf; hold on; 
    set(gcf,'position',[557 75 724 622],'color','w');
    ax(1)=gca;
    im1 = imagesc(ax(1),LS.x,LS.y,flipud(LS.im*1.1)); colormap('gray');
    % set axes properties
    ax(1).XTick=xticks+5e3; ax(1).XTickLabel=string((xticks+5e3)./10^3);
    ax(1).YTick=yticks; ax(1).YTickLabel=string(yticks./10^3);
    set(ax(1),'YDir','normal','XLim',xlimits','YLim',ylimits,'FontSize',14,...
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
    % display contours on ax(1)
    hold on; contour(ax(1),REMA.x,REMA.y,flipud(REMA.im),0:500:2000,'-k','linewidth',2,'ShowText','on');
    
% load & display ITS_LIVE velocity map
    cd([homepath,'data/velocities']);
    v.v = ncread('ANT_G0240_2017.nc','v'); v.v(v.v==-3267)=NaN;
    v.x = ncread('ANT_G0240_2017.nc','x'); v.y = ncread('ANT_G0240_2017.nc','y');
    % display velocity map on ax(2)
    ax(2) = axes; 
    im2=imagesc(ax(2),v.x,v.y,v.v'); colormap(cmocean('haline'));
    im2.AlphaData=0.5; caxis([0 1100]);
    % set axes properties
    ax(2).XLim=xlimits; ax(2).YLim=ylimits; 
    set(ax(2),'YDir','normal','Visible','off','fontsize',12,'XTick',[],'YTick',[]); hold on; 
    % add legend
    l=legend('Position',[0.23 0.8 0.12 0.07],'Color',[175 175 175]/255,...
        'fontname','arial','fontsize',fontsize);
    
% Add lat lon coordinates to top and right axes
    latlon = graticuleps(-66:0.25:-64,-64:0.5:-62,'color',[175 175 175]/255,'linewidth',1,...
        'HandleVisibility','off');
    ax(1).Position = [0.135 0.1 0.73 0.85];
    ax(2).Position=[0.001 0.1 1 0.85];     
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','62.00^oS','textcolor',[175 175 175]/255,...
        'HeadStyle','none','LineStyle', 'none','Position',[.97 .67 0 0],'FontSize',fontsize,'FontName',fontname);     
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','62.50^oS','textcolor',[175 175 175]/255, ...
        'HeadStyle','none','LineStyle', 'none','Position',[.97 .39 0 0],'FontSize',fontsize,'FontName',fontname);     
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','65.25^oW','textcolor',[175 175 175]/255, ...
        'HeadStyle','none','LineStyle', 'none','Position',[.64 .97 0 0],'FontSize',fontsize,'FontName',fontname);     
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','65.00^oW','textcolor',[175 175 175]/255, ...
        'HeadStyle','none','LineStyle', 'none','Position',[.36 .97 0 0],'FontSize',fontsize,'FontName',fontname); 
    
% plot OIB flights
    cd([homepath,'../figures/1_StudyArea']);
    OIBfiles = dir('IRMCR2*.csv');
    for i=1:length(OIBfiles)
        OIB.file = readmatrix(OIBfiles(i).name);
        [OIB.x,OIB.y] = wgs2ps(OIB.file(:,2),OIB.file(:,1),'StandardParallel',-71,'StandardMeridian',0);
        if i==1
            hold on; plot(ax(2),OIB.x,OIB.y,'color','w',...
                'displayname','NASA OIB','linewidth',linewidth-1);
        else
            hold on; plot(ax(2),OIB.x,OIB.y,'color','w',...
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

% plot (b) as outline polygon

% add colorbar
    c = colorbar('Position',[0.29 0.58 0.02 0.15],'fontsize',fontsize,'fontweight',...
        'bold','color','w','fontname',fontname);
    set(get(c,'Title'),'String','Speed (m a^{-1})','color','w','fontname',...
        fontname,'Position',[5 100 0]);

% insert text labels
    % Crane
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Crane', ...
       'HeadStyle','none','LineStyle', 'none', 'TextRotation',70,'Position',[.52 0.65 0 0],...
       'FontSize',fontsize-1,'FontName',fontname,'fontweight','bold','TextColor',[200 200 200]/255);
    % Jorum
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Jorum', ...
        'HeadStyle','none','LineStyle', 'none', 'TextRotation',35,'Position',[.52 .855 0 0],...
        'FontSize',fontsize-1,'FontName',fontname,'fontweight','bold','TextColor',[200 200 200]/255);    
    % Flask
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Flask', ...
       'HeadStyle','none','LineStyle', 'none', 'TextRotation',55,'Position',[.8 .19 0 0],...
       'FontSize',fontsize-1,'FontName',fontname,'fontweight','bold','TextColor',[200 200 200]/255);    
    % Mapple
    %annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Mapple', ...
    %    'HeadStyle','none','LineStyle', 'none', 'TextRotation',66,'Position',[.65 .7 0 0],...
    %    'FontSize',fontsize-1,'FontName',fontname,'fontweight','bold','TextColor',[200 200 200]/255);    
    % Melville
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Melville', ...
        'HeadStyle','none','LineStyle', 'none', 'TextRotation',70,'Position',[.68 .6 0 0],...
        'FontSize',fontsize-1,'FontName',fontname,'fontweight','bold','TextColor',[200 200 200]/255);           
    % Pequod
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Pequod', ...
        'HeadStyle','none','LineStyle', 'none', 'TextRotation',66,'Position',[.74 .59 0 0],...
        'FontSize',fontsize-1,'FontName',fontname,'fontweight','bold','TextColor',[200 200 200]/255);               
    % Starbuck
    annotation('textarrow',[0.6 0.6],[0.6 0.6],'string','Starbuck', ...
        'HeadStyle','none','LineStyle', 'none', 'TextRotation',38,'Position',[.78 .47 0 0],...
        'FontSize',fontsize-1,'FontName',fontname,'fontweight','bold','TextColor',[200 200 200]/255);   
    % Tributary A
    txt = sprintf('A');
        text(-2.413e6,1.232e6,txt,'color',[200 200 200]/255,'fontsize',fontsize-1,'fontweight','bold','fontname',fontname);       
    % Tributary B
    txt = sprintf('B');
        text(-2.415e6,1.237e6,txt,'color',[200 200 200]/255,'fontsize',fontsize-1,'fontweight','bold','fontname',fontname);    
    % Tributary C
    txt = sprintf('C');
        text(-2.417e6,1.251e6,txt,'color',[200 200 200]/255,'fontsize',fontsize-1,'fontweight','bold','fontname',fontname); 
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.925+min(get(gca,'XLim')),(max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        '(a)','fontsize',fontsize,'linewidth',linewidth-1,'backgroundcolor','w','fontname',fontname);          

    % plot LIMA inset in figure
    ax(4)=axes('pos',[0.2 0.12 0.2 0.2]);
    set(ax(4),'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[]);
    cd([homepath,'../Imagery']);
    L = imread('LIMA.jpg');
    imshow(L);
    hold on; plot(1100,4270,'o','color','y','markersize',10,'linewidth',3);
    
% add Landsat NIR band to demonstrate meltwater pooling on Crane surface
cd([homepath,'../Imagery/LC08_L1GT_217106_20200110_20200114_01_T2']);
    landsat = dir('*B4.TIF');
    [LS2.im,LS2.R] = readgeoraster(landsat.name); [LS2.ny,LS2.nx] = size(LS2.im);
    % polar stereographic coordinates of image boundaries
    LS2.x = linspace(min(LS2.R.XWorldLimits),max(LS2.R.XWorldLimits),LS2.nx); 
    LS2.y = linspace(min(LS2.R.YWorldLimits),max(LS2.R.YWorldLimits),LS2.ny);
    % create new figure
    F2 = figure(2); clf; hold on; F2.Position = F1.Position.*2./3;
    % display image on ax(5) 
    ax(5)=gca; ax(5).Position = [0.15 0.15 0.8 0.8];
    im2 = imagesc(ax(5),LS2.x,LS2.y,flipud(LS2.im*1.1)); colormap('gray');
    % set axis properties
    set(ax(5),'YDir','normal','linewidth',2,'fontsize',fontsize,'fontname',fontname,...
        'XLim',[xlimits2(1) xlimits2(end)],'YLim',[ylimits2(1) ylimits2(end)]); 
    ax(5).XTick=xlimits2+6e3; ax(5).XTickLabel=string((xlimits2+6e3)./10^3);
    ax(5).YTick=ylimits2+4700; ax(5).YTickLabel=string((ylimits2+4700)./10^3);
    xlabel('Easting (km)'); ylabel('Northing (km)');  
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.925+min(get(gca,'XLim')),(max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        '(b)','fontsize',fontsize,'linewidth',linewidth-1,'backgroundcolor','w','fontname',fontname);          
    
% save figure 1
if save_figure || merge_figures
    set(F1,'InvertHardCopy','off'); %set(F2,'InvertHardCopy','off'); % save colors as is
    cd([homepath,'figures/']);
    exportgraphics(F1,'studyArea.png','Resolution',300);
    exportgraphics(F2,'meltwaterPonds.png','Resolution',300);
    cd([homepath,'../write-ups/Thesis/figures']);
    exportgraphics(F1,'studyArea.png','Resolution',300);
    exportgraphics(F2,'meltwaterPonds.png','Resolution',300);    
    disp('figures 1 and 2 saved.');    
end

%% Observed Conditions Time Series

close all;

cd([homepath,'inputs-outputs']);

save_figure = 1; % = 1 to save figure
fontsize = 16; % fontsize for plots
font = 'Arial';
linewidth = 2; % line width for plots
markersize = 10; % marker size for plots
col1 = parula(length(2009:2019)); % color scheme for plotting

% 1. Ice surface elevation
h = load('Crane_surfaceElevationObs.mat').h;
col_h = [1 2 3 3 4 4 4 4 4 5 5 5 6 6 7 7 7 8.*ones(1,11) 9*ones(1,7) 10];

% 2. Glacier bed (OIB)
b = load('Crane_ObservedBed_Tate.mat').hb.hb0; % OIB
b0 = load('Crane_flowlineModelInitialization.mat').hb0;
x0 = load('Crane_flowlineModelInitialization.mat').x0;
bathym = load('Crane_bathymetryData.mat').cl_trough; % Rebesco et al. (2014)

% 3. Glacier terminus position
termX = load('LarsenB_centerline.mat').centerline.termx;
termY = load('LarsenB_centerline.mat').centerline.termy;
termx = cl.xi(dsearchn([cl.Xi cl.Yi],[termX' termY']));
termdate = load('LarsenB_centerline.mat').centerline.termdate;
%col_term = [1 3 4 6 8 13.*ones(1,6) 14.*ones(1,7) 15.*ones(1,6) 16*ones(1,14) 17*ones(1,22) 18]; 

% 4. Width-averaged ice surface speeds
U = load('Crane_centerlineSpeedsWidthAveraged_2007-2018.mat').U_widthavg; 
% cut off at terminus
Iu = [2 4 6 8 9 15:20]; % index of annual velocities to use (2007-2017)
%col_U = [1 1 2 2 3 3 4 4 5 5 5 5 5 6:12]+5;

% 5. Glacier width
W = load('Crane_CalculatedWidth.mat').width.W;

% 6. SMB
SMB = load('Crane_downscaledSMB_2002-2019.mat').SMB;

% Plot
figure(2); clf
set(gcf,'Position',[100 100 600 750]);
ax(1)=axes('Position',[0.12 0.7 0.75 0.27],'linewidth',2,'fontsize',fontsize,'fontname',font); % geometry
    hold on; grid on; ylabel('Elevation (m)'); 
    xlim([0 55]); ylim([-1200 1000]); legend('Location','southwest');
    for i=1:length(h)
        if i==16
        else
            term = dsearchn(cl.xi',termx(dsearchn(ConvertSerialYearToDate(termdate'),datenum(h(i).date))));
            plot(ax(1),cl.xi(1:term)./10^3,h(i).surface(1:term),'linewidth',linewidth,...
                'color',col1(col_h(i),:),'HandleVisibility','off');
            plot(ax(1),cl.xi(term)/10^3*[1,1],[interp1(x0,b0,cl.xi(term)) h(i).surface(term)],'--','linewidth',linewidth,...
                'color',col1(col_h(i),:),'HandleVisibility','off');
        end
    end
    plot(ax(1),cl.xi./10^3,b,'--k','linewidth',linewidth,'displayname','b');
    plot(ax(1),x0/10^3,b0,'-k','linewidth',linewidth,'displayname','b_{\mu}');
    % Add text label
    text(51,(max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        '(a)','fontsize',fontsize,'linewidth',linewidth-1,'backgroundcolor','w','fontname',font);    
ax(2)=axes('Position',[0.12 0.38 0.75 0.27],'linewidth',2,'fontsize',fontsize,'fontname',font); % speed
    hold on; grid on; ylabel('Speed (m a^{-1})');
    xlim([0 55]); ylim([0 1600]);
    for i=3:length(Iu)
        plot(ax(2),cl.xi(1:dsearchn(cl.xi',termx(dsearchn(termdate',2006+i))))./10^3,movmean(U(Iu(i)).speed(1:dsearchn(cl.xi',termx(dsearchn(termdate',2006+i)))),2).*3.1536e7,'linewidth',linewidth,'color',col1(i-2,:));
    end
    % Add text label
    text(51,(max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        '(b)','fontsize',fontsize,'linewidth',linewidth-1,'backgroundcolor','w','fontname',font);      
ax(3)= axes('Position',[0.12 0.07 0.75 0.27],'linewidth',2,'fontsize',fontsize,'fontname',font); % speed
    hold on; grid on; xlabel('Distance Along Centerline (km)'); ylabel('Average Annual SMB (m a^{-1})');
    xlim([0 55]); ylim([0 0.6]);set(gca,'clim',[2009 2019]); legend('Location','southwest');
    % Add colorbar
    colormap(col1); 
    c = colorbar('Limits',[2009 2019],'Ticks',[2009 2014 2019],'Position',[.89 .32 .03 .3410],...
        'fontsize',fontsize); 
    smb = zeros(10,length(cl.xi));
    for i=1:length(SMB)-8
        term = dsearchn(cl.xi',termx(dsearchn(termdate',2008+i)));
        smb(i,:) = SMB(i+7).smb_adj;
        plot(ax(3),cl.xi(1:term)./10^3,smb(i,1:term),'color',col1(i,:),'linewidth',linewidth,'HandleVisibility','off');
    end
    plot(ax(3),cl.xi(1:term)./10^3,nanmean(smb(1:end-2,1:term),1),'-k','linewidth',linewidth+1,'displayname','mean');
    % Add text label
    text(51,(max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        '(c)','fontsize',fontsize,'linewidth',linewidth-1,'backgroundcolor','w','fontname',font);          

% Save figure
if save_figure
    cd([homepath,'../write-ups/Thesis/figures']);
    saveas(gcf,'centerlineObservations.png','png'); 
    cd([homepath,'figures/']);
    saveas(gcf,'centerlineObservations.png','png');     
    disp('figure 2 saved.');
end

%% Cumulative Strain Rates & Rate Factor A

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
set(gcf,'Position',[100 100 1000 700],'defaultAxesColorOrder',[[0 0.4 0.8];[0 0 0]]);
yyaxis left; set(gca,'linewidth',2,'fontsize',fontsize,'fontname',fontname);
    xlabel('Distance Along Centerline (km)'); xlim([0 45]); ylim([-0.36 2.27]); 
    %ylabel('$$ \Sigma ( \dot{\eta} ) (s^{-1})$$','Interpreter','latex','fontname','Arial','fontsize',18);
    ylabel('Cumulative Strain Rate');
    hold on; grid on; legend('Location','northwest');
    years = [2007:2011 2013:2018];
    for i=2:length(eta_dot_cum(:,1))
        plot(cl.xi(1:135)/10^3,movmean(eta_dot_cum(i,1:135),2),'-','color',col1(i,:),'linewidth',linewidth,'displayname',num2str(years(i))); drawnow
    end
    %plot(cl.xi(1:135)/10^3,nanmean(eta_dot_cum(:,1:135),1),'-k','linewidth',linewidth+1,'displayname','mean');
yyaxis right; ylabel('Rate Factor (Pa^{-3} a^{-1})'); ylim([1.03e-17 1.75e-17]);
    %cb = colorbar('position',[0.13 0.5 0.02 0.3],'fontname',fontname,...
    %    'fontsize',fontsize-3,'Ticks',0:0.5:1,'TickLabels',[{'2009'},{'2013'},{'2017'}]);
    hold on; grid on; xlim([0 45]); legend('Location','northwest');
    plot(cl.xi(1:135)/10^3,A(1:135)*3.1536e7,'--k','linewidth',linewidth+0.5,'displayname','A');
    plot(cl.xi(1:135)/10^3,A_adj(1:135)*3.1536e7,'-k','linewidth',linewidth+0.5,'displayname','A_{adj}');
   
% save figure
if save_figure
    figure(3); 
    % save in thesis figures folder
    cd([homepath,'../write-ups/Thesis/figures/']);
    saveas(gcf,'rateFactor.png','png');  
    cd([homepath,'figures/']);
    saveas(gcf,'rateFactor.png','png');      
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
betax = load('Crane_betaSolution.mat').xn;
xcf = load('Crane_betaSolution.mat').xcf;
U = load('Crane_betaSolution.mat').Un;
x0 = load('Crane_flowlineModelInitialization.mat').x0;
c0 = load('Crane_flowlineModelInitialization.mat').c0;
load('Crane_SCL_results.mat'); % stress-coupling length results
U_2018 = load('Crane_centerlineSpeedsWidthAveraged_2007-2018.mat').U_widthavg(20).speed;

% plot results
figure(4); clf
set(gcf,'Position',[272 143 1000 700],'defaultAxesColorOrder',[0 0 1; 0 0 0]);
set(gca,'fontsize',fontsize,'fontname',fontname,'linewidth',2);
hold on; grid on; 
legend('position',[0.78 0.15 0.09 0.14]);
xlim([0 45]); xlabel('Distance Along Centerline (km)'); 
yyaxis left; ylabel('\beta (s^{1/m} m^{-1/m})'); 
plot(x0(1:c0)./10^3,beta(1:c0),'-b','linewidth',2,'displayname','\beta');     
yyaxis right; ylabel('U (m a^{-1})'); ylim([0 600]);
plot(cl.xi(1:135)./10^3,U_2018(1:135).*3.1536e7,'k','displayname','U_{obs}','linewidth',linewidth);
plot(betax(1:dsearchn(betax',x0(c0)))./10^3,U(1:dsearchn(betax',x0(c0))).*3.1536e7,'--k','displayname','U_{mod}','linewidth',2);    

if save_figure
    figure(4);
    % save in figures folder
    cd([homepath,'figures/']);
    saveas(gcf,'betaSolution.png','png');  
    % save in thesis figures folder
    cd([homepath,'../write-ups/Thesis/figures/']);
    saveas(gcf,'betaSolution.png','png');  
    disp('figure 4 saved.');    
end
    
%% Sensitivity Tests

close all;

save_figures = 0;    % = 1 to save figure
fontsize = 15;      % font size
fontname = 'Arial'; % font name
linewidth = 2;      % line width
markersize = 10;    % marker size

% load sensitivity test no change variables
cd([homepath,'inputs-outputs']);
load('2100_noChange.mat');
% load observations/parameters
x_cl = load('Crane_centerline.mat').x; y_cl = load('Crane_centerline.mat').y; 
    % define x as distance along centerline
    x0 = zeros(1,length(y_cl));
    for i=2:length(x0)
        x0(i) = sqrt((x_cl(i)-x_cl(i-1))^2+(y_cl(i)-y_cl(i-1))^2)+x0(i-1);
    end
load('flowlineModelInitialization.mat')
smb_mean = load('2100_SMB_mean.mat').smb_mean;

% define time stepping (s)
dt = 0.01*3.1536e7;
t_start = 0*3.1536e7;
t_end = 91*3.1536e7;
t = (t_start:dt:t_end);
clear dt t_start t_end

% Observed discharge & calving front positions
    % Rignot et al. (2004) discharge observations 
    % [deciyear discharge (km^3/a)]
    F_obs = [1996 2.6; 2003+2.75/12 5.6; 2003+3.2/12 6.8; 2003+3.5/12 7.6];
    % convert to Gt/a
    % (km^3/a)*(10^9 m^3/km^3)*(917 kg/m^3)*(10^-3 ton/kg)*(10^-9 Gt/ton) = Gt/a
    F_obs(:,2) = F_obs(:,2).*10^9.*917.*10^-3.*10^-9; % Gt/a
    % load observed speed, width, and modeled thickness
    os.years = [2009 2010 2011 2014 2016 2017];
    os.gl = 134;%[134 140 130 141 153 151]; % index
    W = load('calculatedWidth.mat').width.W;
    os.W = W(os.gl); % m
    h = load('surfaceElevationObs.mat').h;
    os.h = 67;%[h(1).surface(os.gl(1)) h(2).surface(os.gl(2)) h(3).surface(os.gl(3)) h(14).surface(os.gl(4)) h(28).surface(os.gl(5)) h(35).surface(os.gl(6))];
    os.H = (1000/(1000-917)).*os.h; % flotation thickness (m)
    U = load('centerlineSpeedsWidthAveraged_2007-2018.mat').U_widthavg;
    os.U = [U(5).speed(os.gl) U(8).speed(os.gl) U(11).speed(os.gl) U(16).speed(os.gl) U(18).speed(os.gl) U(19).speed(os.gl)];
    % calculate F_obs, convert to Gt/a, and append to vector array
    % (m^3/s)*(3.1536e7 s/a)*(917 kg/m^3)*(10^-3 ton/kg)*(10^-9 Gt/ton)
    % use the mean of a rectangular and ellipsoidal bed
    F_obs = [F_obs; os.years' (os.W.*os.H.*os.U.*3.1536e7*917.*10^-3.*10^-9*pi*1/4)']; % Gt/a;

    % Calving front positions
    termdate = load('LarsenB_centerline.mat').centerline.termdate;
    termX = load('LarsenB_centerline.mat').centerline.termx;
    termY = load('LarsenB_centerline.mat').centerline.termy;
    termx = cl.x(dsearchn([cl.X' cl.Y'],[termX' termY']));

% define color schemes
col1=cmocean('thermal',14); % color scheme for figures 5 & 6 
col1(1,:) = []; col1(end-1:end,:)=[]; % don't use end member colors
colF = [140,81,10; 191,129,45; 178,171,210; 128,115,172; 84,39,136]./255;

% set up axes
loop=1; % loop to minimize plotting commands
while loop==1
    figure(5); clf 
    set(gcf,'Position',[100 200 1000 1000]);
    % SMB
    axA=axes('position',[0.08 0.73 0.25 0.25]); hold on; grid on; % geometry
        set(gca,'fontsize',fontsize,'YTick',-1200:400:1200,'XTickLabel',[]);
        ylabel('Elevation (m)');         
        xlim([25 90]); ylim([-800 600]);
        plot(x1./10^3,hb1,'-k','linewidth',2);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.9+min(get(gca,'YLim')),...
            '(a)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1); 
        % add colorbar
        colormap(col1);
        cb2=colorbar('northoutside'); set(cb2,'fontname',fontname,'fontsize',fontsize-3,'Ticks',0:1/5:1,'TickLabels',string(0:-2:-10));
        set(get(cb2,'label'),'String','\DeltaSMB (m a^{-1})','fontsize',fontsize-3);
    axD=axes('position',[0.08 0.51 0.25 0.2]); hold on; grid on; % thickness
        set(gca,'fontsize',fontsize,'YTick',0:300:1200,'XTickLabel',[]);
        ylabel('Thickness (m)');        
        xlim([25 90]); ylim([0 1000]);
        plot(x1./10^3,hb1,'-k','linewidth',2);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.9+min(get(gca,'YLim')),...
            '(d)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1);         
    axG=axes('position',[0.08 0.29 0.25 0.2]); hold on; grid on; % speed
        set(gca,'fontsize',fontsize,'YTick',0:200:1000,'XTickLabel',[]);        
        ylabel('Speed (m a^{-1})');  
        xlim([25 90]); ylim([175 750]);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.919+min(get(gca,'YLim')),...
            '(g)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1); 
    axJ=axes('position',[0.08 0.07 0.25 0.2]); hold on; grid on; % xgl/xcf
        set(gca,'fontsize',fontsize);
        ylabel('SMB_{mean} (m a^{-1})'); xlabel('Distance Along Centerline (km)');
        xlim([25 90]); ylim([smb_mean(end)*3.1536e7-0.5 smb_mean(1)*3.1536e7+0.5]);              
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.9+min(get(gca,'YLim')),...
            '(j)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1);     
    % TF   
    axB=axes('position',[0.4 0.73 0.25 0.25]); hold on; grid on; % geometry
        set(gca,'fontsize',fontsize,'YTick',-1500:500:1500,'XTickLabel',[]);
        xlim([25 90]); ylim([-800 600]);
        plot(x1./10^3,hb1,'-k','linewidth',2);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.9+min(get(gca,'YLim')),...
            '(b)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1); 
        % add colorbar
        colormap(col1);
        cb1=colorbar('northoutside'); set(cb1,'fontname',fontname,'fontsize',fontsize-3,'Ticks',0:1/5:1,'TickLabels',string(0:0.2:1));
        set(get(cb1,'label'),'String','\DeltaTF (^oC)','fontsize',fontsize-3);
    axE=axes('position',[0.4 0.51 0.25 0.2]); hold on; grid on; % thickness
        set(gca,'fontsize',fontsize,'YTick',0:300:1200,'XTickLabel',[]);
        xlim([25 90]); ylim([0 1000]);               
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.9+min(get(gca,'YLim')),...
            '(e)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1);      
    axH=axes('position',[0.4 0.29 0.25 0.2]); hold on; grid on; % speed
        set(gca,'fontsize',fontsize,'YTick',0:200:1000,'XTickLabel',[]);        
        xlim([25 90]); ylim([175 750]);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.9+min(get(gca,'YLim')),...
            '(h)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1); 
    axK=axes('position',[0.4 0.07 0.25 0.2]); hold on; grid on; % gl/cf
         set(gca,'fontsize',fontsize);
        ylabel('TF (^oC)'); xlabel('Distance Along Centerline (km)');
        xlim([25 90]); ylim([0.1 1.3]);               
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.9+min(get(gca,'YLim')),...
            '(k)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1);         
    % d_fw    
    axC=axes('position',[0.72 0.73 0.25 0.25]); hold on; grid on; % geometry
        set(gca,'fontsize',fontsize,'YTick',-1500:500:1500,'XTickLabel',[]);
        xlim([25 90]); ylim([-800 600]);
        plot(x1./10^3,hb1,'-k','linewidth',2);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.9+min(get(gca,'YLim')),...
            '(c)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1); 
        % add colorbar
        colormap(col1);
        cb3=colorbar('northoutside'); set(cb3,'fontname',fontname,'fontsize',fontsize-3,'Ticks',0:1/5:1,'TickLabels',string(0:2:10));
        set(get(cb3,'label'),'String','\Deltad_{fw} (m)','fontsize',fontsize-3);
    axF=axes('position',[0.72 0.51 0.25 0.2]); hold on; grid on; % thickness
        set(gca,'fontsize',fontsize,'YTick',0:300:1200,'XTickLabel',[]);
        xlim([25 90]); ylim([0 1000]);
        plot(x1./10^3,hb1,'-k','linewidth',2);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.9+min(get(gca,'YLim')),...
            '(f)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1);         
    axI=axes('position',[0.72 0.29 0.25 0.2]); hold on; grid on; % speed
        set(gca,'fontsize',fontsize,'YTick',0:200:1000,'XTickLabel',[]);        
        xlim([25 90]); ylim([175 750]);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.9+min(get(gca,'YLim')),...
            '(i)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1); 
    axL=axes('position',[0.72 0.07 0.25 0.2]); hold on; grid on; % xgl/xcf
        set(gca,'fontsize',fontsize);
        ylabel('d_{fw} (m)'); xlabel('Distance Along Centerline (km)');
        xlim([25 90]); ylim([fwd0-1 fwd0+11]);              
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.9+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.9+min(get(gca,'YLim')),...
            '(l)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1); 
    % axes for other perturbations
    figure(6); clf 
    set(gcf,'Position',[150 200 700 1000]);
    % SMB_enh  
    axAA=axes('position',[0.1 0.73 0.38 0.24]); hold on; grid on; % geometry
        set(gca,'fontsize',fontsize,'YTick',-1500:500:1500,'XTickLabel',[]);
        ylabel('Elevation (m)'); 
        xlim([25 90]); ylim([-800 600]);
        plot(x1./10^3,hb1,'-k','linewidth',2);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.89+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.885+min(get(gca,'YLim')),...
            '(a)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1); 
        % add colorbar
        colormap(col1);
        cb1=colorbar('northoutside'); set(cb1,'fontname',fontname,'fontsize',fontsize-3,'Ticks',0:1/5:1,'TickLabels',string(0:-2:-10));
        set(get(cb1,'label'),'String','\DeltaSMB_{enh} (m a^{-1})','fontsize',fontsize-3);
    axCC=axes('position',[0.1 0.51 0.38 0.2]); hold on; grid on; % thickness
        set(gca,'fontsize',fontsize,'YTick',0:300:1200,'XTickLabel',[]);
        ylabel('Thickness (m)');
        xlim([25 90]); ylim([0 1000]);               
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.89+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.91+min(get(gca,'YLim')),...
            '(c)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1);      
    axEE=axes('position',[0.1 0.29 0.38 0.2]); hold on; grid on; % speed
        set(gca,'fontsize',fontsize,'YTick',0:200:1000,'XTickLabel',[]);        
        ylabel('Speed (m a^{-1})');  
        xlim([25 90]); ylim([175 750]);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.89+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.91+min(get(gca,'YLim')),...
            '(e)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1); 
    axGG=axes('position',[0.1 0.08 0.38 0.2]); hold on; grid on; % gl/cf
         set(gca,'fontsize',fontsize);
        ylabel('SMB_{mean} (m a^{-1})'); xlabel('Distance Along Centerline (km)');
        xlim([25 90]); ylim([smb_mean(end)*3.1536e7-1 smb_mean(1)*3.1536e7+1]);               
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.89+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.91+min(get(gca,'YLim')),...
            '(g)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1);         
    % SMB_enh + TF
    axBB=axes('position',[0.6 0.73 0.38 0.24]); hold on; grid on; % geometry
        set(gca,'fontsize',fontsize,'YTick',-1500:500:1500,'XTickLabel',[]);
        xlim([25 90]); ylim([-800 600]);
        plot(x1./10^3,hb1,'-k','linewidth',2);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.89+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.91+min(get(gca,'YLim')),...
            '(b)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1); 
        % add colorbar
        colormap(col1);
        cb2=colorbar('northoutside'); set(cb2,'fontname',fontname,'fontsize',fontsize-3,'Ticks',0:1/5:1,'TickLabels',string(0:2:10));
        set(get(cb2,'label'),'String','-\DeltaSMB_{enh} (m a^{-1}) and \DeltaTF (^oC)','fontsize',fontsize-3);
    axDD=axes('position',[0.6 0.51 0.38 0.2]); hold on; grid on; % thickness
        set(gca,'fontsize',fontsize,'YTick',0:200:1000,'XTickLabel',[]);
        xlim([25 90]); ylim([0 1000]);
        plot(x1./10^3,hb1,'-k','linewidth',2);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.89+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.91+min(get(gca,'YLim')),...
            '(d)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1);         
    axFF=axes('position',[0.6 0.29 0.38 0.2]); hold on; grid on; % speed
        set(gca,'fontsize',fontsize,'YTick',0:200:1000,'XTickLabel',[]);        
        xlim([25 90]); ylim([175 750]);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.89+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.91+min(get(gca,'YLim')),...
            '(f)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1); 
    axHH=axes('position',[0.6 0.07 0.38 0.2]); hold on; grid on; % xgl/xcf
        set(gca,'fontsize',fontsize);
        xlabel('Distance Along Centerline (km)'); ylabel('TF_{max}');
        xlim([25 90]); ylim([-0.1 1.1]);              
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.89+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.91+min(get(gca,'YLim')),...
            '(h)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1); 
    % ice mass discharge across grounding line
    figure(7); clf;
    set(gcf,'position',[50 300 1000 800]); 
    ax10 = axes('position',[0.08 0.57 0.27 0.39]); hold on; % SMB
        set(gca,'fontsize',fontsize); grid on;
        xlabel('SMB_{mean} (m a^{-1})'); ylabel('F_{gl} (Gt a^{-1})');
        xlim([smb_mean(end)*3.1536e7-0.5 smb_mean(1)*3.1536e7+0.5]); ylim([0.4 1.4]);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.02+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.91+min(get(gca,'YLim')),...
            '(a)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1); 
        % add colorbar
        colormap(col1);
        cb5=colorbar('northoutside','fontname',fontname,'fontsize',fontsize-2,'Ticks',0:1/5:1,'TickLabels',string(0:-2:-10));
        set(get(cb5,'label'),'String','\DeltaSMB (m a^{-1})','fontsize',fontsize-2);
    ax11 = axes('position',[0.39 0.57 0.27 0.39]); hold on; % TF
        set(gca,'fontsize',fontsize); grid on; 
        xlabel('TF (^oC)');  
        xlim([0.1 1.3]); ylim([0.5 1.4]); 
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.02+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.91+min(get(gca,'YLim')),...
            '(b)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1); 
        % add colorbar
        colormap(col1); 
        cb4=colorbar('northoutside','fontname',fontname,'fontsize',fontsize-2,'Ticks',0:1/5:1,'TickLabels',string(0:0.2:1));
        set(get(cb4,'label'),'String','\DeltaTF (^oC)','fontsize',fontsize-2);
    ax12 = axes('position',[0.7 0.57 0.27 0.39]); hold on; % d_fw
        set(gca,'fontsize',fontsize); grid on;
        xlabel('d_{fw} (m)'); 
        xlim([FWD0-0.5 FWD0+10.5]); ylim([0.5 1.4]);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.02+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.91+min(get(gca,'YLim')),...
            '(c)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1); 
         % add colorbar
        colormap(col1);
        cb6=colorbar('northoutside','fontname',fontname,'fontsize',fontsize-2,'Ticks',0:1/5:1,'TickLabels',string(0:2:10)); 
        set(get(cb6,'label'),'String','\Deltad_{fw} (m)','fontsize',fontsize-2);
    ax13 = axes('position',[0.17 0.1 0.3 0.35]); hold on; 
        set(gca,'fontsize',fontsize); grid on; 
        xlabel('SMB_{mean} (m a^{-1})'); ylabel('F_{gl} (Gt a^{-1})'); 
        xlim([smb_mean(end)*3.1536e7-0.5 smb_mean(1)*3.1536e7+0.5]); ylim([0.5 1.4]); 
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.02+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.91+min(get(gca,'YLim')),...
            '(d)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1); 
        % add colorbar
        colormap(col1); 
        cb4=colorbar('eastoutside','fontname',fontname,'fontsize',fontsize-2,'Ticks',0:1/5:1,'TickLabels',string(0:-2:-10));
        set(get(cb4,'label'),'String','\DeltaSMB_{enh} (m a^{-1})','fontsize',fontsize-2);
    ax14 = axes('position',[0.58 0.1 0.3 0.35]); hold on;
        set(gca,'fontsize',fontsize); grid on;
        xlabel('\DeltaTF (^oC)');
        xlim([-0.1 1.1]); ylim([0.5 1.4]);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.02+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.91+min(get(gca,'YLim')),...
            '(e)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1); 
        % add colorbar
        colormap(col1);
        cb5=colorbar('eastoutside','fontname',fontname,'fontsize',fontsize-2,'Ticks',0:1/5:1,'TickLabels',string(0:-2:-10));
        set(get(cb5,'label'),'String','\DeltaSMB_{enh} (m a^{-1})','fontsize',fontsize-2);
    % calving front position and discharge over time
    figure(8); clf;  
    set(gcf,'Position',[300 200 1000 1000]); 
    ax16 = axes('position',[0.06 0.57 0.7 0.4]); hold on;
        set(gca,'fontsize',fontsize); grid on;
        xlim([1995 2100]);
        ylabel('Calving Front Position (km)'); ylim([25 85]);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.02+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.1+min(get(gca,'YLim')),...
            '(a)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1);
    ax17 = axes('position',[0.06 0.09 0.7 0.4]); hold on; 
        set(gca,'fontsize',fontsize); grid on; 
        leg = legend('position',[0.8 0.4 0.17 0.24]);
        xlim([1995 2100]); ylim([0 7.25]);        
        xlabel('Year'); ylabel('Grounding Line Discharge (Gt a^{-1})'); ylim([0 7]);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.02+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.1+min(get(gca,'YLim')),...
            '(b)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1);
    loop=loop+1; % exit loop        
end

% (1) loop through results, split into SMB and FWD change scenarios
cd([homepath,'scripts/3_sensitivityTests/results/1_SMB_FWD/']);
files = dir('*geom.mat');
for i=1:length(files)
    if contains(files(i).name,'SMB0') && contains(files(i).name,'_geom.mat')
        files(i).change = 0;
        files(i).changeIn = 'SMB';
        files(i).name = files(i).name; 
        files(length(files)+1).change = 0;
        files(length(files)).changeIn = 'FWD';
        files(length(files)).name = files(i).name;                        
    elseif contains(files(i).name,'SMR0') && contains(files(i).name,'_geom.mat')
        files(i).change = str2double(files(i).name(regexp(files(i).name,'B')+1:...
            regexp(files(i).name,'_geom')-1)).*3.1536e7; 
        files(i).changeIn = 'SMB';  
    elseif contains(files(i).name,'FWD') && contains (files(i).name,'_geom.mat')
        files(i).change = str2double(files(i).name(regexp(files(i).name,'D_')+2:...
            regexp(files(i).name,'m')-1))-FWD0; % m 
        files(i).changeIn = 'FWD';         
    end
end
    
% sort by changeIn and change
files = struct2table(files);
files = sortrows(files,[8,7]);
% save indices for SMB & SMR
ISMB = flipud(find(strcmp('SMB',table2array(files(:,8)))));
IFWD = find(strcmp('FWD',table2array(files(:,8))));
files = table2struct(files);

% SMB
dsmb = zeros(length(ISMB),1); % change in SMB
dHgl = zeros(length(ISMB),1); % change in thickness at the grounding line
dUgl = zeros(length(ISMB),1); % change in ice speed at the grounding line
dL = zeros(length(ISMB),1); % change in length
dgl = zeros(length(ISMB),1); % change in grounding line position
Ugl = zeros(length(ISMB),1); % speed at grounding line 
Fsmb = zeros(length(ISMB),1); % grounding line discharge
for i=1:length(ISMB)
    load(files(ISMB(i)).name);
    % calculate grounding line discharge
    W = interp1(x0,W0,x2); % interpolate width on spatial grid
    % F (Gt/a) = (U m/s)*(A m^2)*(917 kg/m^3)*(3.1536e7 s/a)*(1e-12 Gt/kg)
    % use the mean of a rectangular and an ellipsoidal bed
    Fsmb(i) = (H2(gl2)*U2(gl2)*W(gl2)*pi*1/4)*917*1e-12*3.1536e7; % Gt/a
    figure(5); hold on;
        % ice surface
        plot(axA,x2(1:c2)./10^3,h2(1:c2),'color',col1(i,:),'linewidth',linewidth-0.5);
        % calving front
        plot(axA,[x2(c2) x2(c2)]./10^3,[h2(c2) h2(c2)-H2(c2)],'color',col1(i,:),'linewidth',linewidth-0.5);
        % floating bed
        plot(axA,x2(gl2:c2)./10^3,h2(gl2:c2)-H2(gl2:c2),'color',col1(i,:),'linewidth',linewidth-0.5);
        % H
        plot(axD,x2(1:c2)/10^3,H2(1:c2),'color',col1(i,:),'linewidth',linewidth-0.5);
        % ice surface speed
        plot(axG,x2(1:c2)./10^3,U2(1:c2).*3.1536e7,'color',col1(i,:),'linewidth',linewidth-0.5);
        % dH
        dHgl(i) = H2(gl2)-H1(gl1); 
        % cf and gl positions
        plot(axJ,x2(gl2)/1e3,smb_mean(i)*3.1536e7,'x','markersize',10,'linewidth',linewidth,'color',col1(i,:));
        plot(axJ,x2(c2)/1e3,smb_mean(i)*3.1536e7,'o','markersize',10,'linewidth',linewidth,'color',col1(i,:));        
    % ice mass discharge
    plot(ax10,smb_mean(i)*3.1536e7,Fsmb(i),'o','markersize',markersize,'linewidth',linewidth,'color',col1(i,:));
    if i==1
        figure(8);
        plot(ax16,termdate,termx/10^3,'xk','markersize',markersize-1,'linewidth',linewidth,'displayName','observed');        
        plot(ax16,t/3.1536e7+2009,movmean(XCF2/10^3,100),'--k','linewidth',linewidth+0.5,'displayName','unperturbed');
        plot(ax17,F_obs(:,1),F_obs(:,2),'xk','markersize',markersize-1,'linewidth',linewidth,'displayName','observed');                
        plot(ax17,t/3.1536e7+2009,movmean(Fgl2,100),'--k','linewidth',linewidth+0.5,'displayName','unperturbed');         
    elseif i==length(ISMB)
        figure(8);
        plot(ax16,t/3.1536e7+2009,movmean(XCF2/10^3,100),'-','color',colF(1,:),'linewidth',linewidth+0.5,'displayName','\DeltaSMB_{max}');        
        plot(ax17,t/3.1536e7+2009,movmean(Fgl2,100),'-','color',colF(1,:),'linewidth',linewidth+0.5,'displayName','\DeltaSMB_{max}'); 
    end
    % save results 
    dsmb(i) = files(ISMB(i)).change; % m/a
    dUgl(i) = (U2(gl2)-U1(gl1)).*3.1536e7; % m/a     
    dL(i) = x2(c2)-x1(c1); % m   
    dgl(i) = x2(gl2)-x1(gl1); % m
    Ugl(i) = U2(gl2)*3.1536e7; % m/a
end
% create table to store result quantities
varNames = {'dsmb','dL','dxgl','dHgl','dUgl','Fsmb'};
T_smb = table(round(dsmb),round(dL)/10^3,round(dgl)/10^3,round(dHgl),round(dUgl),Fsmb,'VariableNames',varNames);

% TF
cd([homepath,'scripts/3_sensitivityTests/results/2_TF/']);
dTF = zeros(length(IFWD),1); % change in TF
dsmr = zeros(length(IFWD),1); % change in smr due to TF
dHgl = zeros(length(IFWD),1); % change in mean thickness
dUgl = zeros(length(IFWD),1); % change in mean ice speed
dL = zeros(length(IFWD),1); % change in length
dgl = zeros(length(IFWD),1); % change in grounding line position
Ugl = zeros(length(IFWD),1); % speed at grounding line 
F = zeros(length(IFWD),1); % grounding line discharge
% define thermal forcing
TF0 = 0.2; % ^oC - estimated from Larsen B icebergs
clear f; f = dir('*geom.mat');
% sort files by perturbation magnitude
for i=1:length(f)
    f(i).change = str2double(f(i).name(regexp(f(i).name,'F')+1:...
                  regexp(f(i).name,'_geom')-1)); % ^oC
end
f = struct2table(f);
f = sortrows(f,7,'ascend');
f = table2struct(f); 
for i=1:length(f)
    load(f(i).name);
    W = interp1(x0,W0,x2); % interpolate width on spatial grid
    % F (Gt/a) = (U m/s)*(A m^2)*(917 kg/m^3)*(3.1536e7 s/a)*(1e-12 Gt/kg)
    % use the mean of a rectangular and an ellipsoidal bed
    F(i) = (H2(gl2)*U2(gl2)*W(gl2)*pi*1/4)*917*1e-12*3.1536e7; % Gt/a
    % extract change in TF
    dTF(i) = f(i).change; % m/a
    % estimate initial melt rate using Eqn from Slater et al. (2020):
    mdot0 = (3*10^-4*-hb0(gl1)*((sum(RO0(1:gl1)))*86400)^0.39 + 0.15)*(TF0^1.18)/86400; % m/s
    figure(5); hold on;
        % ice surface
        plot(axB,x2(1:c2)./10^3,h2(1:c2),'color',col1(i,:),'linewidth',linewidth-0.5);
        % calving front
        plot(axB,[x2(c2) x2(c2)]./10^3,[h2(c2) h2(c2)-H2(c2)],'color',col1(i,:),'linewidth',linewidth-0.5);
        % floating bed
        plot(axB,x2(gl2:c2)./10^3,h2(gl2:c2)-H2(gl2:c2),'color',col1(i,:),'linewidth',linewidth-0.5);
        % H
        plot(axE,x2(1:c2)/10^3,H2(1:c2),'color',col1(i,:),'linewidth',linewidth-0.5);
        % ice surface speed
        plot(axH,x2(1:c2)./10^3,U2(1:c2).*3.1536e7,'color',col1(i,:),'linewidth',linewidth-0.5);
        % dH
        dHgl(i) = H2(gl2)-H1(gl1); 
        % cf and gl positions
        plot(axK,x2(gl2)/1e3,TF0+dTF(i),'x','markersize',markersize,...
            'linewidth',linewidth,'color',col1(i,:));
        plot(axK,x2(c2)/1e3,TF0+dTF(i),'o','markersize',markersize,...
            'linewidth',linewidth,'color',col1(i,:));   
    % calving front position and ice mass discharge
    figure(7);
    plot(ax11,TF0+dTF(i),F(i),'o','markersize',markersize,'linewidth',linewidth,'color',col1(i,:));
    if i==length(f)
        figure(8);
        plot(ax16,t/3.1536e7+2009,movmean(XCF2/10^3,100),'-','color',colF(2,:),'linewidth',linewidth+0.5,'displayName','\DeltaTF');        
        plot(ax17,t/3.1536e7+2009,movmean(Fgl2,100),'-','color',colF(2,:),'linewidth',linewidth+0.5,'displayName','\DeltaTF_{max}');
    end
    % save results for data table
    % calculate additional melt due to the increase in subglacial discharge
    dsmr(i) = ((3*10^-4*(h2(gl2)-H2(gl2))*((sum(RO0)*86400)^0.39) +0.15)*((TF0+dTF(i))^1.18)-mdot0)*3.1536e7; % m/a
    dUgl(i) = (U2(gl2)-U1(gl1)).*3.1536e7; % m/a     
    dL(i) = x2(c2)-x1(c1); % m   
    dgl(i) = x2(gl2)-x1(gl1); % m
    Ugl(i) = U2(gl2)*3.1536e7; % m/a
end
varNames = {'dTF','dSMR','dL','dgl','dHgl','dUgl','Fgl'};
T_TF = table(dTF,dsmr,round(dL)/10^3,round(dgl)/10^3,round(dHgl),round(dUgl),F,'VariableNames',varNames);

% FWD
cd([homepath,'scripts/3_sensitivityTests/results/1_SMB_FWD/']);
dfwd = zeros(length(IFWD),1); % change in FWD
dHgl = zeros(length(IFWD),1); % change in mean thickness
dUgl = zeros(length(IFWD),1); % change in mean ice speed
dL = zeros(length(IFWD),1); % change in length
dgl = zeros(length(IFWD),1); % change in grounding line position
Ugl = zeros(length(IFWD),1); % speed at grounding line 
Ffwd = zeros(length(IFWD),1); % grounding line discharge
col1 = cmocean('thermal',length(IFWD)+3); col1(1,:)=[];col1(end-1:end,:)=[];
for i=1:length(IFWD)
    % load file
    load(files(IFWD(i)).name);
    % calculate grounding line discharge
    W = interp1(x0,W0,x2); % interpolate width on spatial grid
    % F (Gt/a) = (U m/s)*(A m^2)*(917 kg/m^3)*(3.1536e7 s/a)*(1e-12 Gt/kg)
    % use the mean of a rectangular and an ellipsoidal bed
    Ffwd(i) = (H2(gl2)*U2(gl2)*W(gl2))*pi*1/4*917*1e-12*3.1536e7; % Gt/a    
    % plot
    figure(5); hold on;
        % ice surface
        plot(axC,x2(1:c2)./10^3,h2(1:c2),'color',col1(i,:),'linewidth',linewidth-0.5);
        % calving front
        plot(axC,[x2(c2) x2(c2)]./10^3,[h2(c2)-H2(c2) h2(c2)],'color',col1(i,:),'linewidth',linewidth-0.5);
        % floating bed
        plot(axC,x2(gl2:c2)./10^3,h2(gl2:c2)-H2(gl2:c2),'color',col1(i,:),'linewidth',linewidth-0.5);
        % thickness
        plot(axF,x2(1:c2)/10^3,H2(1:c2),'color',col1(i,:),'linewidth',linewidth-0.5);
        % ice surface speed
        plot(axI,x2(1:c2)./10^3,U2(1:c2).*3.1536e7,'color',col1(i,:),'linewidth',linewidth-0.5);
        % dH
        dHgl(i) = H2(gl2)-H1(gl1); % m
        % cf and gl positions
        plot(axL,x2(gl2)/10^3,files(IFWD(i)).change+FWD0,'x',...
            'markersize',10,'linewidth',linewidth,'color',col1(i,:));
        plot(axL,x2(c2)/10^3,files(IFWD(i)).change+FWD0,'o',...
            'markersize',10,'linewidth',linewidth,'color',col1(i,:));        
    % ice mass discharge
    figure(7); 
    plot(ax12,files(IFWD(i)).change+FWD0,Ffwd(i),'o','markersize',markersize,'linewidth',linewidth,'color',col1(i,:));    
    % calving front position and grounding line discharge over time
    if i==length(IFWD)
        figure(8);
        plot(ax16,t/3.1536e7+2009,movmean(XCF2/10^3,100),'-','color',colF(3,:),'linewidth',linewidth+0.5,'displayName','\Deltad_{fw,max}');        
        plot(ax17,t/3.1536e7+2009,movmean(Fgl2,100),'-','color',colF(3,:),'linewidth',linewidth+0.5,'displayName','\Deltad_{fw,max}'); 
    end
    % save results 
    dfwd(i) = files(IFWD(i)).change; % m/a
    dUgl(i) = (U2(gl2)-U1(gl1)).*3.1536e7; % m/a     
    dL(i) = x2(c2)-x1(c1); % m   
    dgl(i) = x2(gl2)-x1(gl1); % m
    Ugl(i) = U2(gl2)*3.1536e7; % m/a
end
% create table to store result quantities
varNames = {'dfwd','dL','dxgl','dHgl','dUgl','Ffwd'};
T_fwd = table(dfwd,round(dL)/10^3,round(dgl)/10^3,round(dHgl),round(dUgl),Ffwd,'VariableNames',varNames);

% (3) SMB_enh
cd([homepath,'scripts/3_sensitivityTests/results/3_SMB_enh/']);
dsmb_enh = zeros(length(IFWD),1); % change in SMB_enh
dHgl = zeros(length(IFWD),1); % change in mean thickness
dUgl = zeros(length(IFWD),1); % change in mean ice speed
dL = zeros(length(IFWD),1); % change in length
dgl = zeros(length(IFWD),1); % change in grounding line position
Ugl = zeros(length(IFWD),1); % speed at grounding line 
F = zeros(length(IFWD),1); % grounding line discharge
% define thermal forcing
TF0 = 0.2; % ^oC - estimated from Larsen B icebergs
clear f; f = dir('*geom.mat');
% sort files by perturbation magnitude
for i=1:length(f)
    f(i).change = str2double(f(i).name(regexp(f(i).name,'B')+1:...
                  regexp(f(i).name,'_geom')-1)).*3.1536e7; % m/a
end
f = struct2table(f);
f = sortrows(f,7,'descend');
f = table2struct(f); 
for i=1:length(f)
    load(f(i).name);
    W = interp1(x0,W0,x2); % interpolate width on spatial grid
    % F (Gt/a) = (U m/s)*(A m^2)*(917 kg/m^3)*(3.1536e7 s/a)*(1e-12 Gt/kg)
    % use the mean of a rectangular and an ellipsoidal bed
    F(i) = (H2(gl2)*U2(gl2)*W(gl2)*pi*1/4)*917*1e-12*3.1536e7; % Gt/a
    % estimate initial melt rate using Eqn from Slater et al. (2020):
    mdot0 = (3*10^-4*-hb0(gl1)*((sum(RO0(1:gl1)))*86400)^0.39 + 0.15)*TF0^1.18/86400; % m/s
    figure(6); hold on;
        % ice surface
        plot(axAA,x2(1:c2)./10^3,h2(1:c2),'color',col1(i,:),'linewidth',linewidth-0.5);
        % calving front
        plot(axAA,[x2(c2) x2(c2)]./10^3,[h2(c2) h2(c2)-H2(c2)],'color',col1(i,:),'linewidth',linewidth-0.5);
        % floating bed
        plot(axAA,x2(gl2:c2)./10^3,h2(gl2:c2)-H2(gl2:c2),'color',col1(i,:),'linewidth',linewidth-0.5);
        % H
        plot(axCC,x2(1:c2)/10^3,H2(1:c2),'color',col1(i,:),'linewidth',linewidth-0.5);
        % ice surface speed
        plot(axEE,x2(1:c2)./10^3,U2(1:c2).*3.1536e7,'color',col1(i,:),'linewidth',linewidth-0.5);
        % dH
        dHgl(i) = H2(gl2)-H1(gl1); 
        % cf and gl positions
        plot(axGG,x2(gl2)/1e3,smb_mean(i)*3.1536e7,'x','markersize',markersize,...
            'linewidth',linewidth,'color',col1(i,:));
        plot(axGG,x2(c2)/1e3,smb_mean(i)*3.1536e7,'o','markersize',markersize,...
            'linewidth',linewidth,'color',col1(i,:));   
    % calving front position and ice mass discharge
    figure(7);
    plot(ax13,smb_mean(i)*3.1536e7,F(i),'o','markersize',markersize,'linewidth',linewidth,'color',col1(i,:));
    if i==length(f)
        figure(8);
        plot(ax16,t/3.1536e7+2009,movmean(XCF2/10^3,100),'-','color',colF(4,:),'linewidth',linewidth+0.5,'displayName','\DeltaSMB_{enh,max}');        
        plot(ax17,t/3.1536e7+2009,movmean(Fgl2,100),'-','color',colF(4,:),'linewidth',linewidth+0.5,'displayName','\DeltaSMB_{enh,max}');
    end
    % save results 
    dsmb_enh(i) = f(i).change; % m/a
    % calculate additional melt due to the increase in subglacial discharge
    dsmr(i) = ((3*10^-4*(h2(gl2)-H2(gl2))*((-dsmb_enh(i)*86400)^0.39) +0.15)*(TF0^1.18)-mdot0)*3.1536e7; % m/a
    dUgl(i) = (U2(gl2)-U1(gl1)).*3.1536e7; % m/a     
    dL(i) = x2(c2)-x1(c1); % m   
    dgl(i) = x2(gl2)-x1(gl1); % m
    Ugl(i) = U2(gl2)*3.1536e7; % m/a
end
varNames = {'dSMB_enh','dSMR','dL','dgl','dHgl','dUgl','Fgl'};
T_smb_enh = table(dsmb_enh,dsmr,round(dL)/10^3,round(dgl)/10^3,round(dHgl),round(dUgl),F,'VariableNames',varNames);

% (4) SMB_enh + TF
cd([homepath,'scripts/3_sensitivityTests/results/4_SMB_enh+TF/']);
dsmb_enh = zeros(length(IFWD),1); % change in SMB_enh
dTF = zeros(length(IFWD),1); % change in TF
dsmr = zeros(length(IFWD),1); % change in SMR
dHgl = zeros(length(IFWD),1); % change in mean thickness
dUgl = zeros(length(IFWD),1); % change in mean ice speed
dL = zeros(length(IFWD),1); % change in length
dgl = zeros(length(IFWD),1); % change in grounding line position
Ugl = zeros(length(IFWD),1); % speed at grounding line 
F = zeros(length(IFWD),1); % grounding line discharge
clear f; f = dir('*geom.mat');
% sort files by perturbation magnitude
for i=1:length(f)
    f(i).change_TF = str2double(f(i).name(regexp(f(i).name,'F')+1:...
                  regexp(f(i).name,'_SMB')-1)); % ^oC
    f(i).change_SMB = str2double(f(i).name(regexp(f(i).name,'B')+1:...
                  regexp(f(i).name,'_geom')-1))*3.1536e7; % m/a
end
f = struct2table(f);
f = sortrows(f,7,'ascend');
f = table2struct(f); 
for i=1:length(f)
    load(f(i).name);
    W = interp1(x0,W0,x2); % interpolate width on spatial grid
    % F (Gt/a) = (U m/s)*(A m^2)*(917 kg/m^3)*(3.1536e7 s/a)*(1e-12 Gt/kg)
    % use the mean of a rectangular and an ellipsoidal bed
    F(i) = (H2(gl2)*U2(gl2)*W(gl2)*pi*1/4)*917*1e-12*3.1536e7; % Gt/a
    % estimate initial melt rate using Eqn from Slater et al. (2020):
    mdot0 = (3*10^-4*-hb0(gl1)*((sum(RO0(1:gl1)))*86400)^0.39 + 0.15)*TF0^1.18/86400; % m/s
    figure(6); hold on;
        % ice surface
        plot(axBB,x2(1:c2)./10^3,h2(1:c2),'color',col1(i,:),'linewidth',linewidth-0.5);
        % calving front
        plot(axBB,[x2(c2) x2(c2)]./10^3,[h2(c2) h2(c2)-H2(c2)],'color',col1(i,:),'linewidth',linewidth-0.5);
        % floating bed
        plot(axBB,x2(gl2:c2)./10^3,h2(gl2:c2)-H2(gl2:c2),'color',col1(i,:),'linewidth',linewidth-0.5);
        % H
        plot(axDD,x2(1:c2)/10^3,H2(1:c2),'color',col1(i,:),'linewidth',linewidth-0.5);
        % ice surface speed
        plot(axFF,x2(1:c2)./10^3,U2(1:c2).*3.1536e7,'color',col1(i,:),'linewidth',linewidth-0.5); 
        % cf and gl positions
        plot(axHH,x2(gl2)/1e3,f(i).change_TF,'x',...
            'markersize',markersize,'linewidth',linewidth,'color',col1(i,:));
        plot(axHH,x2(c2)/1e3,f(i).change_TF,'o',...
            'markersize',markersize,'linewidth',linewidth,'color',col1(i,:));        
    % ice mass discharge
    figure(7);
    plot(ax14,f(i).change_TF,F(i),'o','markersize',markersize,'linewidth',linewidth,'color',col1(i,:));
    if i==length(f)
        figure(8);
        plot(ax16,t/3.1536e7+2009,movmean(XCF2/10^3,100),'-','color',colF(5,:),'linewidth',linewidth+0.5,'displayName','\DeltaSMB_{enh,max} & \DeltaTF_{max}');                
        plot(ax17,t/3.1536e7+2009,movmean(Fgl2,100),'-','color',colF(5,:),'linewidth',linewidth+0.5,'displayName','\DeltaSMB_{enh,max} & \DeltaTF_{max}');
    end  
    % save results 
    dsmb_enh(i) = f(i).change_SMB; % m/a
    dTF(i) = f(i).change_TF; % m/a
    % calculate additional melt due to the increase in subglacial discharge
    dsmr(i) = ((3*10^-4*(h2(gl2)-H2(gl2))*((-dsmb_enh(i)*86400)^0.39) +0.15)*((TF0+dTF(i))^(1.18))-mdot0)*3.1536e7; % m/a
    dUgl(i) = (U2(gl2)-U1(gl1)).*3.1536e7; % m/a     
    dL(i) = x2(c2)-x1(c1); % m   
    dgl(i) = x2(gl2)-x1(gl1); % m
    Ugl(i) = U2(gl2)*3.1536e7; % m/a
end
varNames = {'dSMB_enh','dTF','dSMR','dL','dgl','dHgl','dUgl','Fgl'};
T_smb_enh_TF = table(dsmb_enh,dTF,dsmr,round(dL)/10^3,round(dgl)/10^3,round(dHgl),round(dUgl),F,'VariableNames',varNames);

% save figures 
if save_figures
    cd([homepath,'../write-ups/JGlacPaper/']);
    figure(5); saveas(gcf,'sensitivityTests_geom+speed_independent.png','png');
    figure(6); saveas(gcf,'sensitivityTests_geom+speed_enhanced.png','png');
    figure(7); saveas(gcf,'sensitivityTests_discharge.png','png');
    figure(8); saveas(gcf,'sensitivityTests_FglXcf.png','png');    
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

%% Width segments

save_figure = 1; % = 1 to save figure
linewidth = 2; 
fontsize = 18;
fontname = 'Arial';

% Load Landsat image, width, width segments, and glacier extent polygon
cd([homepath,'data/Imagery/']);
ls = dir('LC08*20201104_01_T2_B8.TIF');
[LS.im,R] = readgeoraster(ls.name); [LS.ny,LS.nx] = size(LS.im);
% Polar stereographic coordinates of image boundaries
LS.x = linspace(min(R.XWorldLimits),max(R.XWorldLimits),LS.nx);
LS.y = linspace(min(R.YWorldLimits),max(R.YWorldLimits),LS.ny);
cd([homepath,'inputs-outputs/']);
extx = load('Crane_calculatedWidth.mat').width.extx;
exty = load('Crane_calculatedWidth.mat').width.exty;
W = load('Crane_calculatedWidth.mat').width.W;
ol = load('Crane_glacierOutline.mat').ol;
        
% Plot
col = flipud(cmocean('ice',5)); % color scheme for potting
figure(9); clf
set(gcf,'units','pixels','position',[200 200 1000 800]);
ax1 = axes('position',[0.08 0.1 0.6 0.85]);
    hold on; imagesc(LS.x/10^3,LS.y/10^3,flipud(LS.im)); colormap("gray");
    set(gca,'fontsize',fontsize,'linewidth',linewidth); 
    xlabel('Easting (km)'); ylabel('Northing (km)'); 
    legend('Location','east','color',[0.8,0.8,0.8]);
    xlim([-2.43e3 -2.385e3]); ylim([1.21e3 1.285e3]); 
    fill(ol.x/10^3,ol.y/10^3,col(1,:),'displayname','glacier extent');
    for i=1:length(extx)
        if i==1
            plot(extx(i,:)/10^3,exty(i,:)/10^3,'color',col(2,:),'linewidth',linewidth,'displayname','width segments');
        else
            plot(extx(i,:)/10^3,exty(i,:)/10^3,'color',col(2,:),'linewidth',linewidth,'HandleVisibility','off');        
        end
    end
    plot(cl.X/10^3,cl.Y/10^3,'color',col(3,:),'linewidth',linewidth,'displayname','centerline');
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.02+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.05+min(get(gca,'YLim')),...
            '(a)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1); 
ax2 = axes('position',[0.73 0.2 0.2 0.65]);
    hold on; set(gca,'fontsize',fontsize,'linewidth',linewidth,'YTick',[],'XDir','reverse'); 
    xlabel('Width (km)'); yyaxis right; ylabel('Distance Along Centerline (km)'); 
    plot(W/10^3,cl.xi/10^3,'-k','linewidth',linewidth); grid on; 
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.95+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.05+min(get(gca,'YLim')),...
            '(b)','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1); 
        
% save figure
if save_figure
    cd([homepath,'figures/']);
    saveas(gcf,'Crane_widthSegments.png','png');
    cd([homepath,'../write-ups/Thesis/figures/']);
    saveas(gcf,'Crane_widthSegments.png','png');
    disp('figure 9 saved');
end

%% 2018 Model Misfits
    
% Note: Must rerun first section before running this section each time

save_figure = 0; % = 1 to save figure
plotTimeSteps = 1; % = 1 to plot model time steps

% define home path in directory and add necessary paths
homepath = '/Users/raineyaberle/Desktop/Research/CraneModeling/CraneGlacier_flowlinemodeling/';
addpath([homepath,'scripts/3_sensitivityTests/']); % add path to U_convergence
addpath([homepath,'matlabFunctions/cmocean_v2.0/cmocean']);
addpath([homepath,'inputs-outputs/']);

dx0 = 200; % grid spacing (m)

% Load Crane Glacier initialization variables
b=1; % use a while loop to make this section collapsible
while b==1
    load('flowlineModelInitialization.mat');
    F0=0; % inner boundary flux (m^3/s)

    % Centerline
    cd([homepath,'inputs-outputs']);
    cl.X = load('Crane_centerline.mat').x; cl.Y = load('Crane_centerline.mat').y;
    % define x as distance along centerline
    cl.x = zeros(1,length(cl.X));
    for i=2:(length(cl.X))
        cl.x(i)=sqrt((cl.X(i)-cl.X(i-1))^2+(cl.Y(i)-cl.Y(i-1))^2)+cl.x(i-1);
    end

    % Load observed conditions
    % ice surface
    h_obs = load('surfaceElevationObs.mat').h;
    % terminus position 
    term = load('terminusPositions_2002-2019.mat').term;
    for i=1:length(term)
        termx_obs(i) = term(i).x;
        termDate_obs(i) = term(i).decidate;
    end
    % fit a quadratic function to the terminus positions to smooth seasonal variations
    termx_obs = feval(fit(termDate_obs',termx_obs','poly2'),termDate_obs');
    term_obs = interp1(termDate_obs',termx_obs,2009:2017);
    clear term 
    % ice speed
    U_obsi = load('centerlineSpeedsWidthAveraged_2007-2018.mat').U_widthavg;
    u = [6 8 9 14 15:20]; % indices of speeds to use annually (2009-2017)
    for i=1:length(u)
        U_obs(i).U = U_obsi(u(i)).speed;
        U_obs(i).date = U_obsi(u(i)).date;
    end
    clear U_obsi u 

    % densities and g
    rho_i = 917; % ice density (kg m^-3)
    rho_sw = 1028; % ocean water density (kg m^-3)
    rho_fw = 1000; % fresh water density (kg m^-3)
    g = 9.81; % acceleration (m s^-2)

    % stress parameters (unitless)
    m = 3; % basal sliding exponent
    n = 3; % flow law exponent
    E = 1; % enhancement factor

    % calving parameters
    Hc = 400; % m -> set the calving front to a default minimum ice thickness value if calving criteria is not met
    sigma_b = 0; % back pressure (Pa) - similar to that employed at Helheim Glacier (Nick et al., 2009)

    % maximum & minimum thickness & speed cut-off to check for instability
    H_max = 2000; % maximum thickness (m)
    H_min = 100;  % minimum thickness (m)
    U_min = 100;  % minimum mean speed (m a^-1)

    % initial conditions
    H0 = h0-hb0; % ice thickness (m)
    dUdx0 = [(U0(1:end-1)-U0(2:end))./(x0(1:end-1)-x0(2:end)) 0]; % strain rate (1/s)
    % find the location of the grounding line and the end of the ice-covered domain
    Hf = -(rho_sw./rho_i).*hb0; % flotation thickness (m)
    gl0 = find(Hf-H0>0,1,'first')-1; % grounding line location 
    H0(gl0+1:length(x0))=h0(gl0+1:end)*rho_sw/(rho_sw-rho_i); % buoyant thickness using surface
    H0(H0>=(h0-hb0))=h0(H0>=(h0-hb0))-hb0(H0>=(h0-hb0)); % thickness can't go beneath bed elevation
    H0(c0+1:end) = 0; % zero ice thickness past calving front

    % extend variables to model length if necessary (set extended points to last value)
    A0(find(isnan(A0),1,'first'):end) = A0(find(isnan(A0),1,'first')-1);
    U0(find(isnan(U0),1,'first'):end) = U0(find(isnan(U0),1,'first')-1);
    dUdx0(find(isnan(dUdx0),1,'first'):end) = dUdx0(find(isnan(dUdx0),1,'first')-1);
    beta0(find(isnan(beta0),1,'first'):end) = beta0(find(isnan(beta0),1,'first')-1);
    hb0(find(isnan(hb0),1,'first'):end) = hb0(find(isnan(hb0),1,'first')-1);
    W0(find(isnan(W0),1,'first'):end) = W0(find(isnan(W0),1,'first')-1);    

    % define time stepping (s)
    dt = 0.01*3.1536e7;
    t_start = 0*3.1536e7;
    t_end = 9*3.1536e7;
    t = (t_start:dt:t_end);

    % initialize ice flux
    Fgl=zeros(1,length(t)); % ice mass flux across the grounding line
    
    b=b+1;
end
    
% initialize variables
x=x0; U=U0; W=W0; gl=gl0; dUdx=dUdx0; A=A0; h=h0; hb=hb0; H=H0; fwd=FWD0; beta=beta0;

% run flowline model
for i=1:length(t)

    % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
    Rxx = 2*nthroot(dUdx./(E.*A),n); % resistive stress (Pa)
    crev_s = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*fwd); % surface crevasse penetration depth (m)
    Hab = H+rho_sw/rho_i*(hb); % height above buoyancy (m)
    crev_b = rho_i/(rho_sw-rho_i).*(Rxx./(rho_i*g)-Hab); % basal crevasse depth (m)
    % calving front located where the inland-most crevasse intersects sea level
    if i==1 % use observed calving front position for first iteration
        xcf = x0(c0);
    else
        if length(h)>=find(h-crev_s<0,1,'first')+1
            xcf_s = interp1(h(find(h-crev_s<0,1,'first')-1:find(h-crev_s<0,1,'first')+1)...
                -crev_s(find(h-crev_s<0,1,'first')-1:find(h-crev_s<0,1,'first')+1),...
                x(find(h-crev_s<0,1,'first')-1:find(h-crev_s<0,1,'first')+1),0,'linear','extrap'); % (m along centerline)
        else
            xcf_s = interp1(h-crev_s,x,0,'linear','extrap');
        end
        if length(h)>=find(h-crev_b<0,1,'first')+1
            xcf_b = interp1(h(find(h-crev_b<0,1,'first')-1:find(h-crev_b<0,1,'first')+1)...
                -crev_b(find(h-crev_b<0,1,'first')-1:find(h-crev_b<0,1,'first')+1),...
                x(find(h-crev_b<0,1,'first')-1:find(h-crev_b<0,1,'first')+1),0,'linear','extrap'); % (m along centerline)
        else
            xcf_b = interp1(h-crev_b,x,0,'linear','extrap');
        end
        if xcf_s<0; xcf_s=NaN; end
        if xcf_b<0; xcf_b=NaN; end
        % calving front = whichever calving criteria occurs the
        % furthest inland
        if xcf_s<xcf_b
            xcf = xcf_s;
        else
            xcf = xcf_b;
        end
        xcf=xcf_s;
        if xcf<20e3 || xcf > 70e3 || isnan(xcf)
            xcf = interp1(feval(fit(x',(h-crev_s)','poly1'),x),x,0,'linear','extrap');
            %xcf = x0(c0);
        end
    end

    % calculate the thickness required to remain grounded at each grid cell
    Hf = -(rho_sw./rho_i).*hb; % flotation thickness (m)
    % find the location of the grounding line and use a floating
    % geometry from the grounding line to the calving front
    if ~isempty(find(Hf-H>0,1,'first'))
        if length(Hf)>=find(Hf-H>0,1,'first')+1
            xgl = interp1(Hf(find(Hf-H>0,1,'first')-1:find(Hf-H>0,1,'first')+1)...
                -H(find(Hf-H>0,1,'first')-1:find(Hf-H>0,1,'first')+1),...
                x(find(Hf-H>0,1,'first')-1:find(Hf-H>0,1,'first')+1),0,'linear','extrap'); % (m along centerline)
        else
            xgl = x(find(Hf-H>0,1,'first')-1);
        end
    else
        xgl=xcf;
    end
    if xgl>xcf % grounding line can't be past calving front
        xgl=xcf;
    end
    XGL(i) = xgl; % save grounding line position over time

    % create coordinate system that hits cf and gl exactly
    % has resolution dxmax near the ice divide
    % has resolution dxmin from gl to c
    % and has smooth variation between
    xl = round(xgl/dx0); %number of ideal grid spaces needed to reach the grounding line
    dx = xgl/xl; %new grid spacing (should be ~dx0)
    xn = 0:dx:xgl; %new distance vector
    if xcf-xgl > 0
        xl = round((xcf-xgl)/dx0);
        dx = (xcf-xgl)/xl;
        xn = [xn xn(end)+dx:dx:xcf];
    end
    clear dx; dxn = [xn(2:end)-xn(1:end-1) xn(end)-xn(end-1)];

    % get geometry on new coordinates
    c = length(xn); gl = dsearchn(xn',xgl); % indeces for cf and gl
    %if the crevasses never intersect sea level
    if isempty(c) == 1 % set the calving front to a default minimum ice thickness value
        c = find(H<Hc,1,'first');
    end
    if isempty(c)==1 % set c to length of x if criteria still not met
        c=length(x);
        disp('calving criteria not met');
    end

    hb = interp1(x0,hb0,xn,'linear','extrap');
    W = interp1(x0,W0,xn,'linear','extrap');
    H = interp1(x,H,xn,'linear','extrap');
    U = interp1(x,U,xn,'linear','extrap');
    A = interp1(x0,A0,xn,'linear','extrap');
    beta = interp1(x0,beta0,xn,'linear','extrap');
    x = xn; dx = dxn;

    % calculate surface elevation
    h = hb+H; % surface elevation (m a.s.l.)
    h(gl+1:c) = (1-rho_i/rho_sw).*H(gl+1:c); %adjust the surface elevation of ungrounded ice to account for buoyancy
    H(h<0)=0-hb(h<0); h(h<0)=0; % surface cannot go below sea level
    h(h-H<hb) = hb(h-H<hb)+H(h-H<hb); % thickness cannot go beneath bed elevation

    % plot geometry, speed, & grounding line and calving front positions
    col = parula(length(t)+20); % color scheme for plots
    if plotTimeSteps
        if t(i)==t_start
            figure(1); clf
            set(gcf,'Position',[0 100 1300 400]);
            ax1 = axes('Position',[0.06 0.12 0.27 0.78]); % glacier geometry
            hold on; grid on;
            set(gca,'FontSize',14,'linewidth',2);
            %title('Glacier Geometry');
            legend('Location','northeast');
            xlim([0 75]); ylim([min(hb)-100 max(h)+200]);
            xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)');
            % ice surface
            plot(x(1:c)./10^3,h(1:c),'color',col(i,:),'linewidth',2,'displayname','2009');
            % calving front
            plot(x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
            % floating bed
            plot(x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off');
            % bed elevation
            plot(x0./10^3,hb0,'k','linewidth',2,'HandleVisibility','off');
            % mean sea level
            plot([x(1),x(end)]/10^3,[0,0],'k--','HandleVisibility','off');
            ax2 = axes('Position',[0.38 0.1 0.28 0.8]); % ice speed
            hold on; grid on;
            set(gca,'FontSize',14,'linewidth',2);
            title('Ice Speed'); legend('Location','northeast');
            xlim([0 65]); ylim([0 1500]);
            xlabel('Distance Along Centerline (km)'); ylabel('Speed (m yr^{-1})');
            % ice speed
            plot(x(1:c)./10^3,U(1:c).*3.1536e7,'color',col(i,:),'linewidth',2,'displayname','2009');
            ax3 = axes('Position',[0.7 0.1 0.28 0.8]); % calving front position
            hold on; grid on;
            set(gca,'FontSize',14,'linewidth',2);
            title('Terminus & Grounding Line Positions'); legend('Location','best');
            xlim([30 65]); ylim([0 t_end./3.1536e7]);
            xlabel('Distance Along Centerline (km)'); ylabel('Year');
            % terminus & grounding line positions
            plot(x(c)/10^3,t(i)./3.1536e7,'.','markersize',15,'color',col(i,:),'displayname','2009');
            plot(ax3,x(gl)./10^3,t(i)./3.1536e7,'x','Color',col(i,:),'markersize',10,'linewidth',2,'HandleVisibility','off');
        elseif mod(i-1,round(length(t)/10))==0 % display every length(t)/10
            figure(1);
            % ice surface
            plot(ax1,x(1:c)/10^3,h(1:c),'-','color',col(i,:),'linewidth',2,'displayname',num2str(round(t(i)./3.1536e7)+2009));
            plot(ax1,x(gl:c)/10^3,h(gl:c)-H(gl:c),'-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
            plot(ax1,[x(c);x(c)]/10^3,[h(c);h(c)-H(c)],'-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
            % calving front
            plot(ax1,x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
            % floating bed (gl:c)
            plot(ax1,x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off');
            % ice speed
            plot(ax2,x(1:c)/10^3,U(1:c).*3.1536e7,'-','Color',col(i,:),'linewidth',2,'DisplayName',num2str(round(t(i)./3.1536e7)+2009));
            % calving front position
            plot(ax3,x(c)/10^3,t(i)/3.1536e7,'.','Color',col(i,:),'markersize',15,'displayname',num2str(round(t(i)./3.1536e7)+2009)); hold on;
            % grounding line position
            plot(ax3,x(gl)/10^3,t(i)/3.1536e7,'x','Color',col(i,:),'markersize',10,'linewidth',2,'HandleVisibility','off'); hold on;
        end
    end
    
    % calculate the effective pressure (ice overburden pressure minus water
    % pressure) assuming an easy & open connection between the ocean and
    % ice-bed interface
    sl = find(hb<=0,1,'first'); % find where the glacier base first drops below sea level
    N_ground = rho_i*g*H(1:sl); % effective pressure where the bed is above sea level (Pa)
    N_marine = rho_i*g*H(sl+1:length(x))+(rho_sw*g*hb(sl+1:length(x))); % effective pressure where the bed is below sea level (Pa)
    N = [N_ground N_marine];
    N(N<0)=0; % cannot have negative values

    % Solve for new velocity
    [U,dUdx,Td,Tlatb,Tlon,vm,G_minus,G,G_plus,Hm,gamma,eta,Am] = U_convergence(x,U,U0,dUdx,H,h,A,E,N,W,dx,c,n,m,beta,rho_i,rho_sw,g,sigma_b,i);
   
    % calculate ice flux
    F = U.*H.*W; % ice flux (m^3 s^-1)
    F(isnan(F))=0;
    F(1)=F(2)+F0;

    % implement SMB, SMR, delta_SMB, & delta_SMR
    if t(i)/3.1536e7<10 % use original SMB & SMR for first 10 model years
        smr = zeros(1,c);
        for k=gl+1:c
            smr(k) = smr0-0.001*(smr0)*(k-gl+1);
        end
        smb = interp1(x0,SMB0+Q0,x);
        % plot
%                 if i==1
%                     figure(2); clf
%                     set(gcf,'position',[200 300 1000 400]);
%                     subplot(1,3,1); hold on; 
%                         set(gca,'fontsize',12,'fontweight','bold','linewidth',2);
%                         plot(x/10^3,smb.*3.1536e7,'color',col(i,:),'linewidth',2);
%                         xlabel('km along centerline'); ylabel('m a^{-1}'); grid on;
%                         title('SMB');
%                     subplot(1,3,2); hold on; 
%                         set(gca,'fontsize',12,'fontweight','bold','linewidth',2);
%                         plot(x/10^3,smr.*3.1536e7,'color',col(i,:),'linewidth',2);
%                         xlabel('km along centerline'); ylabel('m a^{-1}'); grid on;
%                         title('SMR');
%                     subplot(1,3,3); hold on;
%                         set(gca,'fontsize',12,'fontweight','bold','linewidth',2);
%                         plot(t(i)/3.1536e7+2009,fwd,'.','markersize',10,'color',col(i,:));
%                         xlabel('Year'); ylabel('m'); grid on;
%                         title('d_{fw}');
%                 end        
    elseif t(i)/3.1536e7>=10 % implement changes after 10 model years
        smr = zeros(1,c);
        % increase SMR linearly at each time increment to reach delta_smr by 2100
        delta_smri = delta_smr/(2100-2019)*(t(i)/3.1536e7-10); % total increase in smr from 2019 rate
        for k=gl+1:c
            smr(k) = smr0+delta_smri-0.001*(smr0+delta_smri)*(k-gl+1);
        end
        delta_smbi = delta_smb/(2100-2019)*(t(i)/3.1536e7-10); % total increase in smb from 2019 rate 
        smb = interp1(x0,smb0+Q0,x);
        for k=1:c
            smb(k) = smb(k)+delta_smbi*(h0(1)-h(k))/(h0(1)-h0(c0)); 
        end
        % plot
%                 if mod(t(i)/3.1536e7,10)==0
%                     figure(2);
%                     subplot(1,3,1);
%                         plot(x/10^3,smb.*3.1536e7,'color',col(i,:),'linewidth',2);
%                     subplot(1,3,2);
%                         plot(x/10^3,smr.*3.1536e7,'color',col(i,:),'linewidth',2);
%                     subplot(1,3,3);
%                         plot(t(i)/3.1536e7+2009,fwd,'.','markersize',10,'color',col(i,:));
%                 end
    end

    % calculate the  change in ice thickness from continuity
    clearvars dHdt
    dHdt(1) = (-1/W(1))*(F(1)-F(2))/(x(1)-x(2)); % forward difference
    dHdt(2:c-1) = (-1./W(2:c-1)).*(F(1:c-2)-F(3:c))./(x(1:c-2)-x(3:c)); % central difference
    dHdt(c:length(x)) = (-1./W(c:length(x))).*(F(c-1:length(x)-1)-F(c:length(x)))./(x(c-1:length(x)-1)-x(c:length(x))); % backward difference
    dH = dHdt.*dt;

    % new thickness (change from dynamics, SMB, & SMR)
    Hn = H+dH+(smb.*dt)+(smr.*dt);
    Hn(Hn < 0) = 0; % remove negative values
    H = Hn; % set as the new thickness value

    % stop the model if it behaves unstably (monitored by ice thickness and speed)
    if max(H) > H_max
        disp(['Adjust dt']);
        break;
    end
    if mean(U) < U_min/3.1536e7
        disp('Too slow!');
        break;
    end
    if any(~isfinite(H(1:c))) || any(~isfinite(U(1:c))) || any(~isfinite(h(1:c)))
        disp('non finite values');
        break;
    end

end

% calculate and plot misfits
figure(2); clf;
set(gcf,'position',[200 200 1000 700],'defaultAxesColorOrder',[[0 0 0];[0.8 0.1 0.1]]);
ax1 = axes('position',[0.08 0.67 0.36 0.3]); hold on; grid on;
    set(gca,'fontsize',18,'fontname','Arial','linewidth',2); 
    xlim([0 52]); ylabel('Misfit (m)'); 
    plot(x/10^3,h-interp1(cl.x,h_obs(36).surface,x),'-','linewidth',2,...
        'color',[0.8 0.1 0.1],'HandleVisibility','off');
    plot(x/10^3,nanmean(h-interp1(cl.x,h_obs(36).surface,x))*ones(1,length(x)),...
        '--','linewidth',2,'color',[0.8 0.1 0.1],'HandleVisibility','off');
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.92+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.93+min(get(gca,'YLim')),...
            '(a)','backgroundcolor','w','fontsize',18,'linewidth',1);  
ax2 = axes('position',[0.56 0.67 0.36 0.3]); hold on; grid on;
    set(gca,'fontsize',18,'fontname','Arial','linewidth',2); 
    xlim([0 52]); ylabel('Misfit (m a^{-1})');
    plot(x/10^3,(U-interp1(cl.x,U_obs(10).U,x))*3.1536e7,'-','linewidth',2,...
        'color',[0.8 0.1 0.1],'HandleVisibility','off');
    plot(x/10^3,nanmean(U-interp1(cl.x,U_obs(10).U,x))*3.1536e7*ones(1,length(x)),...
        '--','linewidth',2,'color',[0.8 0.1 0.1],'HandleVisibility','off');
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.92+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.93+min(get(gca,'YLim')),...
            '(b)','backgroundcolor','w','fontsize',18,'linewidth',1); 
set(gcf,'position',[200 200 1000 700],'defaultAxesColorOrder',[[0 0 0];[0 0.4 0.8]]);
ax3 = axes('position',[0.08 0.1 0.36 0.5]); hold on; grid on; legend('Location','west'); 
    set(gca,'fontsize',18,'fontname','Arial','linewidth',2); 
    xlim([0 52]); ylabel('Elevation (m)'); xlabel('Distance Along Centerline (km)'); 
    plot(x/10^3,movmean(h,20),'-k','linewidth',2,'displayname','h_{mod}');
    plot(cl.x(1:150)/10^3,h_obs(36).surface(1:150),'--k','linewidth',2,'displayname','h_{obs}');
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.92+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.93+min(get(gca,'YLim')),...
            '(c)','backgroundcolor','w','fontsize',18,'linewidth',1);     
    %yyaxis right; set(ax3,'YTick',[],'YTickLabel',[]);
ax4 = axes('position',[0.56 0.1 0.36 0.5]); hold on; grid on; legend('Location','west');
    set(gca,'fontsize',18,'fontname','Arial','linewidth',2); 
    xlim([0 52]); xlabel('Distance Along Centerline (km)'); 
    yyaxis left; ylabel('Speed (m a^{-1})'); ylim([100 850]);
        plot(x/10^3,U*3.1536e7,'-k','linewidth',2,'displayname','U_{mod}');
        plot(cl.x(1:145)/10^3,U_obs(10).U(1:145)*3.1536e7,'--k','linewidth',2,'displayname','U_{obs}');
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.92+min(get(gca,'XLim')),...
                (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.93+min(get(gca,'YLim')),...
                '(d)','backgroundcolor','w','fontsize',18,'linewidth',1);
    yyaxis right; ylabel('Basal Roughness Factor (s^{1/m} m^{-1/m})'); ylim([0 1.4]);
        plot(x/10^3,movmean(beta,5),'-','linewidth',2,'color',[0 0.4 0.8],'displayname','\beta');

% Plot stresses
figure(3); clf; 
set(gcf,'position',[324    83   913   369]);
col2 = parula(4); 
hold on; grid on; legend; 
set(gca,'fontsize',14,'linewidth',2); xlabel('x (km)'); 
yyaxis left; ylabel('Td (Pa)');plot(x/10^3,Td,'linewidth',2,'color',col2(1,:),'displayname','T_d');
yyaxis right; ylabel('Pa'); plot(x/10^3,Tlon,'linewidth',2,'color',col2(2,:),'displayname','T_{lon}');
plot(x/10^3,Tlatb,'linewidth',2,'color',col2(3,:),'displayname','T_{lat,b}');

% Plot viscosity, A, U and dUdx
figure(4); clf;
set(gcf,'position',[322   331   920   366]);
subplot(1,2,1);
hold on; grid on; legend;
set(gca,'fontsize',14,'linewidth',2); xlabel('x (km)'); 
yyaxis left; ylabel('vm (Pa s)'); plot(x/10^3,vm,'-m','linewidth',2,'displayname','vm');
yyaxis right; ylabel('A (Pa^{-3} a^{-1})'); plot(x/10^3,A*3.1536e7,'-b','linewidth',2,'displayname','A');
plot(x/10^3,5.6e-17*ones(1,length(x)),'--b','displayname','A (Nick et al., 2010)','linewidth',2);
subplot(1,2,2);
hold on; grid on; legend;
set(gca,'fontsize',14,'linewidth',2); xlabel('x (km)'); 
yyaxis left; ylabel('m/y'); plot(x/10^3,U*3.1536e7,'-k','linewidth',2,'displayname','U');
yyaxis right; ylabel('1/y'); plot(x/10^3,dUdx*3.1536e7,'--k','linewidth',2,'displayname','dUdx');

% plot SMB, SMR, H
figure(5); clf; 
set(gca,'fontsize',14,'linewidth',2);
hold on; grid on; legend; xlabel('x (km)');
yyaxis left; ylabel('(m/yr)');
plot(x/10^3,smb.*3.1536e7,'-k','displayname','SMB','linewidth',2);
plot(x/10^3,smr.*3.1536e7,'--k','displayname','SMR','linewidth',2);
yyaxis right; ylabel('H (m)');
plot(x/10^3,H,'-b','displayname','H','linewidth',2);
plot(x/10^3,dH,'--b','displayname','dH','linewidth',2);

% display grounding line misfit
disp('Mean misfits at the grounding line:');
disp(['    Speed: ',num2str((U(gl)-interp1(cl.x,U_obs(10).U,x(gl)))*3.1536e7),' m/yr']);
disp(['    Elevation: ',num2str((h(gl)-interp1(cl.x,h_obs(36).surface,x(gl)))),' m']);

% save figure
if save_figure
    cd([homepath,'figures/']);
    saveas(gcf,'misfits2018.png','png');
    cd([homepath,'../write-ups/Thesis/figures/']);
    saveas(gcf,'misfits2018.png','png');
    disp('figure saved');
end
