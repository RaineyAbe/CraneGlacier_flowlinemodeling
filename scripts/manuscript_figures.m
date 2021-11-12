%% Create Figures for MS Thesis and JGlac Submission
% RKA 2021

% Figures:
%   - Map of the study area
%   - Observed Conditions Time Series
%   - Cumulative Strain Rates / A adjustment
%   - Beta Solution
%   - Sensitivity Tests
%   - Regional glacier terminus time series
%   - 2018 Model Misfits

close all; clear all;

% Define home path and add paths to necessary functions
homepath = '/Users/raineyaberle/Desktop/Research/CraneModeling/CraneGlacier_flowlinemodeling/';
addpath([homepath,'matlabFunctions/']);
addpath([homepath,'matlabFunctions/gridLegend_v1.4/']);
addpath([homepath,'matlabFunctions/cmocean_v2.0/cmocean/']);
addpath([homepath,'inputs-outputs/']);
addpath([homepath,'scripts/']);
addpath([homepath,'scripts/1_synthesize_datasets/OIBPicking/functions/']);
addpath([homepath,'data/OIB/L1B']);

% Load Crane centerline
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
fontsize = 16;      % font size
fontname = 'Arial'; % font name
markersize = 8;    % marker size
linewidth = 1.5;      % line width

% Add paths to functions
    addpath([homepath,'matlabFunctions/']);
    addpath([homepath,'matlabFunctions/cmocean_v2.0/cmocean/']);
    addpath([homepath,'matlabFunctions/AntarcticMappingTools/']);

% Specify xticks/yticks coordinates and axes limits
    xticks = linspace(-2.455e6,-2.375e6,5); yticks=linspace(1.21e6,1.29e6,5);
    xlimits=[xticks(1) xticks(end)]; ylimits = [yticks(1) yticks(end)];
    xlimits2 = -2.4196e6:0.005e6:-2.4046e6; ylimits2 = 1.2460e6:0.005e6:1.2610e6; 
    
% load & display Landsat image 1
    cd([homepath,'data/Imagery/']);
    landsat = dir('LC08_L1GT_217106_20200110_20200114_01_T2_B8.TIF');
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
    xlabel('Easting [km]'); ylabel('Northing [km]');
    
% load & display REMA as contours
    cd([homepath,'data/Imagery/']);
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
    set(ax(2),'YDir','normal','Visible','off','fontsize',fontsize,'XTick',[],'YTick',[]); hold on; 
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
    cd([homepath,'data/OIB/L2/']);
    OIBfiles = dir('*IRMCR2*.csv');
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
    
% plot (b) as outline polygon
    r = rectangle('Position',[xlimits2(1) ylimits2(1) xlimits2(end)-xlimits2(1) ylimits2(end)-ylimits2(1)],...
        'EdgeColor','k','linewidth',linewidth+0.5);
    
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
        'bold','color','w','fontname',fontname);
    set(get(c,'Title'),'String','Speed [m a^{-1}]','color','w','fontname',...
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
    % a - panel label
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.925+min(get(gca,'XLim')),(max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' a ','fontsize',fontsize,'linewidth',linewidth-1,'backgroundcolor','w','fontname',fontname,'fontweight','bold');          

    % plot LIMA inset in figure
    ax(4)=axes('pos',[0.2 0.12 0.2 0.2]);
    set(ax(4),'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[]);
    cd([homepath,'data/Imagery/']);
    L = imread('LIMA.jpg');
    imshow(L);
    hold on; plot(1100,4270,'o','color','y','markersize',10,'linewidth',3);
    
% add zoomed in Landsat B8 figure to demonstrate meltwater pooling on Crane surface
    % create new figure
    F2 = figure(2); clf; hold on; F2.Position = [200 200 450 390];
    % display image on ax(5) 
    ax(5)=gca; ax(5).Position = [0.15 0.15 0.8 0.8];
    im2 = imagesc(ax(5),LS.x,LS.y,flipud(LS.im*1.1)); colormap('gray');
    % set axis properties
    set(ax(5),'YDir','normal','linewidth',2,'fontsize',fontsize,'fontname',fontname,...
        'XLim',[xlimits2(1) xlimits2(end)],'YLim',[ylimits2(1) ylimits2(end)]); 
    ax(5).XTick=xlimits2-0.4e3; ax(5).XTickLabel=string((xlimits2-0.4e3)./10^3);
    ax(5).YTick=ylimits2+4e3; ax(5).YTickLabel=string((ylimits2+4e3)./10^3);
    xlabel('Easting [km]'); ylabel('Northing [km]');  
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.925+min(get(gca,'XLim')),(max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        ' b ','fontsize',fontsize,'linewidth',linewidth-1,'backgroundcolor','w','fontname',fontname,'fontweight','bold');          
    % Draw boxes over melt features        
    rectangle(ax(5),'Position',[-2.4109e6 1.2565e6 0.0009e6 0.002e6],...
        'EdgeColor','w','linewidth',linewidth);
    rectangle(ax(5),'Position',[-2.4149e6 1.2469e6 0.0055e6 0.0032e6],...
        'EdgeColor','w','linewidth',linewidth);
    rectangle(ax(5),'Position',[-2.4117 1.2588 0.0005 0.0007]*1e6,...
        'EdgeColor','w','linewidth',linewidth);
    
% save figures 1 and 2
if save_figure
    set(F1,'InvertHardCopy','off'); %set(F2,'InvertHardCopy','off'); % save colors as is
    cd([homepath,'figures/']);
    exportgraphics(F1,'studyArea.png','Resolution',600);
    exportgraphics(F2,'meltwaterPonds.png','Resolution',600);   
    disp('figures 1 and 2 saved.');    
end

%% Centerline Observations Time Series

close all;

cd([homepath,'inputs-outputs']);

save_figure = 1; % = 1 to save figure
fontsize = 16; % fontsize for plots
font = 'Arial';
linewidth = 1.5; % line width for plots
markersize = 10; % marker size for plots
col1 = parula(length(2009:2019)); % color scheme for plotting

% 1. Ice surface elevation
h = load('surfaceElevationObs.mat').h;
col_h = [1 2 3 3 4 4 4 4 4 5 5 5 6 6 7 7 7 8.*ones(1,11) 9*ones(1,7) 10];

% 3. Glacier terminus position
termX = load('LarsenB_centerline.mat').centerline.termx;
termY = load('LarsenB_centerline.mat').centerline.termy;
termx = cl.xi(dsearchn([cl.Xi cl.Yi],[termX' termY']));
termdate = load('LarsenB_centerline.mat').centerline.termdate;
%col_term = [1 3 4 6 8 13.*ones(1,6) 14.*ones(1,7) 15.*ones(1,6) 16*ones(1,14) 17*ones(1,22) 18]; 

% 2. Glacier bed (OIB)
b = load('observedBed.mat').HB.hb0; % OIB
b0 = load('flowlineModelInitialization.mat').hb0;
x0 = load('flowlineModelInitialization.mat').x0;
W0 = load('flowlineModelInitialization.mat').W0;
bathym = load('bathymetryData.mat').cl_trough; % Rebesco et al. (2014)
% Model floating bed using observed surface and width-averaged bed
    rho_sw = 1000; % density of sea water (kg/m^3)
    rho_i = 917; % density of ice (kg/m^3)
    % calculate the thickness required to remain grounded at each grid cell
    Hf = -(rho_sw./rho_i).*b0; % flotation thickness (m)
    % Loop through observed surface elevations
    Hi = zeros(length(h),length(x0)); hi = zeros(length(h),length(x0));
    c = zeros(1,length(h)); gl = zeros(1,length(h));
    figure(10); clf; hold on; grid on; col = parula(length(h)+4); legend; 
    plot(x0/10^3,b0,'-k','linewidth',2,'displayname','b_0');
    for i=1:length(h)
        hi(i,:) = interp1(h(i).x(~isnan(h(i).surface)),h(i).surface(~isnan(h(i).surface)),x0);
        c(i) = dsearchn(x0',termx(dsearchn(ConvertSerialYearToDate(termdate'),datenum(h(i).date))));
        hi(i,c(i):end) = 0;
        Hi(i,:) = hi(i,:) - b0;
        % find the location of the grounding line and use a floating
        % geometry from the grounding linU_widthavge to the calving front
        if any(find(Hf(~isnan(Hi(i,:)))-Hi(i,(~isnan(Hi(i,:))))>0,1,'first'))
            gl(i) = find(Hf(~isnan(Hi(i,:)))-Hi(i,(~isnan(Hi(i,:))))>0,1,'first'); % (m along centerline)
        else
            gl(i) = NaN;
        end
        %if gl(i)<200
        %    gl(i) = NaN;
        %end
        % calculate floating thickness using surface elevation, 
        % adjust thickness to account for width-averaging assuming a
        % parabolic bed shape
        if ~isnan(gl(i))
            % calculate floating bed using surface elevation
            Hi(i,gl(i)+1:c(i)) = (rho_sw/(rho_sw-rho_i)).*hi(i,gl(i)+1:c(i)); 
            % calculate parabolic cross-sectional area 
            % y = a(x-h)^2+k --> Hn = a(x-(W/2))^2-H, where a = 4H/(W^2)
            for j=gl(i)+1:c(i)
                xj = 0:10:W0(j);
                a = 4*Hi(i,j)/(W0(j)^2);
                Hi(i,j) = nanmean(-(a.*(xj-W0(i)/2).^2-Hi(i)));
            end
            Hi(i,c(i)+1:end) = 0;
            plot(x0(gl(i):c(i))/10^3,hi(i,gl(i):c(i))-Hi(i,gl(i):c(i)),'--','color',col(i,:),'linewidth',2,'HandleVisibility','off');
        else
            Hi(i,:) = NaN;
        end
        plot(x0/10^3,hi(i,:),'-','displayname',num2str(i),'linewidth',2,'color',col(i,:));
    end

% 4. Width-averaged ice surface speeds
U = load('centerlineSpeedsWidthAveraged_2007-2018.mat').U_widthavg; 
% cut off at terminus
Iu = [2 4 6 8 9 15:20]; % index of annual velocities to use (2007-2017)
%col_U = [1 1 2 2 3 3 4 4 5 5 5 5 5 6:12]+5;

% 5. Glacier width
W = load('calculatedWidth.mat').width.W;

% 6. SMB
SMB = load('downscaledClimateVariables_2009-2019.mat').SMB;
RO = load('flowlineModelInitialization.mat').RO0;

% Plot
figure(2); clf
set(gcf,'Position',[-1000 100 600 750]);
ax(1)=axes('Position',[0.12 0.7 0.75 0.27],'linewidth',2,'fontsize',fontsize,'fontname',font); % geometry
    hold on; grid on; ylabel('Elevation [m]'); 
    xlim([0 55]); ylim([-1200 1000]); legend('Location','southwest');
    for i=1:length(h)
        if i~=16
            plot(ax(1),x0(1:c(i))./10^3,hi(i,1:c(i)),'linewidth',linewidth,...
                'color',col1(col_h(i),:),'HandleVisibility','off');
            plot(ax(1),x0(c(i))/10^3*[1,1],[b0(c(i)) 50],'--','linewidth',linewidth,...
                'color',col1(col_h(i),:),'HandleVisibility','off');
        end
    end
    plot(ax(1),cl.xi./10^3,b,'--k','linewidth',linewidth,'displayname','b');
    plot(ax(1),x0/10^3,b0,'-k','linewidth',linewidth,'displayname','b_{\mu}');
    % Add text label
    text(51,(max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        'a','fontsize',fontsize,'linewidth',linewidth-1,'backgroundcolor','w','fontname',font,'fontweight','bold');    
ax(2)=axes('Position',[0.12 0.38 0.75 0.27],'linewidth',2,'fontsize',fontsize,'fontname',font); % speed
    hold on; grid on; ylabel('Speed [m a^{-1}]');
    xlim([0 55]); ylim([0 1600]);
    for i=3:length(Iu)
        plot(ax(2),cl.xi(1:dsearchn(cl.xi',termx(dsearchn(termdate',2006+i))))./10^3,movmean(U(Iu(i)).speed(1:dsearchn(cl.xi',termx(dsearchn(termdate',2006+i)))),2).*3.1536e7,'linewidth',linewidth,'color',col1(i-2,:));
    end
    % Add text label
    text(51,(max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        'b','fontsize',fontsize,'linewidth',linewidth-1,'backgroundcolor','w','fontname',font,'fontweight','bold');      
ax(3)= axes('Position',[0.12 0.07 0.75 0.27],'linewidth',2,'fontsize',fontsize,'fontname',font); % SMB
    hold on; grid on; xlabel('Distance Along Centerline [km]'); ylabel('Average Annual SMB [m a^{-1}]');
    xlim([0 55]); ylim([0.3 0.7]);set(gca,'clim',[2009 2019]); legend('Location','northwest');
    % Add colorbar
    colormap(col1); 
    colorbar('Limits',[2009 2019],'Ticks',[2009.5 2014 2018.5],'TickLabels',string(2009:5:2019),'Position',[.89 .32 .03 .3410],...
        'fontsize',fontsize); 
    smb = zeros(8,length(cl.xi));
    for i=1:length(smb(:,1))
        term = dsearchn(cl.xi',termx(dsearchn(termdate',2008+i)));
        smb(i,:) = SMB.linear(i,:);
        plot(ax(3),cl.xi(1:term)./10^3,smb(i,1:term),'color',col1(i,:),'linewidth',linewidth,'HandleVisibility','off');
    end
    plot(ax(3),cl.xi(1:term)./10^3,SMB.downscaled_average_linear(1:term),'--k','linewidth',linewidth,'displayname','SMB_{\mu}');
    plot(ax(3),x0/10^3,interp1(cl.xi(1:term),SMB.downscaled_average_linear(1:term),x0) - RO0.*3.1536e7,'-k','linewidth',linewidth+0.5,'displayname','SMB_{\mu} - RO');
    % Add text label
    text(51,(max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.925+min(get(gca,'YLim')),...
        'c','fontsize',fontsize,'linewidth',linewidth-1,'backgroundcolor','w','fontname',font,'fontweight','bold');          

% Save figure
if save_figure
    figure(2);
    cd([homepath,'figures/']);
    exportgraphics(gcf,'centerlineObservations.png','Resolution',600);     
    disp('figure 2 saved.');
end

%% Regional Glacier Termini Time Series

close all;

cd([homepath,'data/terminus/regional/']);

save_figure = 1; % = 1 to save figure
fontsize = 16; 
fontname = 'Arial';
linewidth = 2; 

% add path to required functions
addpath([homepath,'matlabFunctions/']);
addpath([homepath,'matlabFunctions/line2arrow-kakearney-pkg-8aead6f/']);
addpath([homepath,'matlabFunctions/line2arrow-kakearney-pkg-8aead6f/line2arrow/']);
addpath([homepath,'matlabFunctions/line2arrow-kakearney-pkg-8aead6f/axescoord2figurecoord/']);
addpath([homepath,'matlabFunctions/line2arrow-kakearney-pkg-8aead6f/parsepv']);

% load Landsat images
cd([homepath,'data/Imagery/']);
% LSA 
landsatA = dir('*217105*B8.TIF');
    [LSA.im,LSA.R] = readgeoraster(landsatA.name); [LSA.ny,LSA.nx] = size(LSA.im);
    % polar stereographic coordinates of image boundaries
    LSA.x = linspace(min(LSA.R.XWorldLimits),max(LSA.R.XWorldLimits),LSA.nx); 
    LSA.y = linspace(min(LSA.R.YWorldLimits),max(LSA.R.YWorldLimits),LSA.ny);
% LSB
landsatB = dir('*20211013_01_T2_B8.tif');
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
figure(3); clf; 
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
    text(940,3950,' a ','edgecolor','k','fontsize',fontsize-3,'linewidth',1,'backgroundcolor','w','fontweight','bold');
    text(1000,4120,' b ','edgecolor','k','fontsize',fontsize-3,'linewidth',1,'backgroundcolor','w','fontweight','bold');
    text(1013,4371,' c ','edgecolor','k','fontsize',fontsize-3,'linewidth',1,'backgroundcolor','w','fontweight','bold');    
    text(1080,4460,' d ','edgecolor','k','fontsize',fontsize-3,'linewidth',1,'backgroundcolor','w','fontweight','bold');    
    text(1085,4550,' e ','edgecolor','k','fontsize',fontsize-3,'linewidth',1,'backgroundcolor','w','fontweight','bold');    
    % add colorbar
    colormap(ax1,'parula'); 
    colorbar('Limits',[0 1],'Ticks',[0 0.5 1],'TickLabels',[{'2014'},{'2017'},{'2021'}],...
        'Location','southoutside','fontsize',fontsize-1,'FontName',fontname); 
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
clear f; files(1).lon = NaN; files(1).lat = NaN; files(1).date = NaN;
% loop through files to grab terminus positions
for j=1:length(files)
    if ~contains(files(j).name,'CL')
        file = shaperead(files(j).name);
        % loop through file
        for k=1:length(file)
            files(length(files)+1).lon = file(k).X; % Lon
            files(length(files)).lat = file(k).Y; % Lat
            % convert to polar stereographic coordinates
            [files(length(files)).X,files(length(files)).Y] = wgs2ps(files(length(files)).lon,files(length(files)).lat,'StandardParallel',-71,'StandardMeridian',0);
            files(length(files)).x = x(dsearchn([cl.X' cl.Y'],[nanmean(files(length(files)).X) nanmean(files(length(files)).Y)])); % x
        end
    end 
end
% plot
col1 = parula(length(files)); % color scheme for plotting
ax2=axes('position',[0.39 0.58 0.23 0.4]); hold on; grid on;      
    set(gca,'fontname',fontname,'fontsize',fontsize);
    xlabel('Easting [km]'); ylabel('Northing [km]');
    % plot LIMA
    colormap(ax2,'gray');
    imagesc(LSA.x/10^3,LSA.y/10^3,flipud(LSA.im)); 
    ax2.XLim=[-2.4508e3 -2.443e3]; ax2.YLim=[1.4138e3 1.4204e3]; 
    % plot centerline
    l1=line(cl.X(3:8)/10^3,cl.Y(3:8)/10^3,'color','k','linewidth',linewidth);
    % arrow
    line2arrow(l1,'color','k','linewidth',linewidth,'headwidth',20,'headlength',20);
    % plot terminus positions
    for j=1:length(files)
        plot(files(j).X/10^3,files(j).Y/10^3,'color',col1(j,:),'linewidth',linewidth);
    end
    % plot label
    text((ax2.XLim(2)-ax2.XLim(1))*0.88+ax2.XLim(1),(max(ax2.YLim)-min(ax2.YLim))*0.925+min(ax2.YLim),...
        ' a ','edgecolor','k','fontsize',fontsize,'linewidth',1,'backgroundcolor','w','fontweight','bold'); 
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
clear f; files(1).lon = NaN; files(1).lat = NaN; files(1).date = NaN;
% loop through files to grab terminus positions
for j=1:length(files)
    if ~contains(files(j).name,'CL')
        file = shaperead(files(j).name);
        % loop through file
        for k=1:length(file)
            files(length(files)+1).lon = file(k).X; % Lon
            files(length(files)).lat = file(k).Y; % Lat
            % convert to polar stereographic coordinates
            [files(length(files)).X,files(length(files)).Y] = wgs2ps(files(length(files)).lon,files(length(files)).lat,'StandardParallel',-71,'StandardMeridian',0);
            files(length(files)).x = x(dsearchn([cl.X' cl.Y'],[nanmean(files(length(files)).X) nanmean(files(length(files)).Y)])); % x
        end
    end 
end
% plot
col1 = parula(length(files)); % color scheme for plotting
ax3=axes('position',[0.72 0.58 0.24 0.4]); hold on; grid on;      
    set(gca,'fontname',fontname,'fontsize',fontsize);
    xlabel('Easting [km]'); ylabel('Northing [km]');
    % plot LIMA
    colormap(ax3,'gray');
    imagesc(LSA.x/10^3,LSA.y/10^3,flipud(LSA.im)); 
    ax3.XLim = [-2.4416e3 -2.4247e3]; ax3.YLim = [1.3552e3 1.3721e3];        
    % plot centerline
    l3=line(cl.X(4:7)/10^3,cl.Y(4:7)/10^3,'color','k','linewidth',linewidth);
    % arrow
    line2arrow(l3,'color','k','linewidth',linewidth,'headwidth',20,'headlength',20);        
    % plot terminus positions
    for j=1:length(files)
        plot(files(j).X/10^3,files(j).Y/10^3,'color',col1(j,:),'linewidth',linewidth);
    end
    % plot label
    text((ax3.XLim(2)-ax3.XLim(1))*0.88+ax3.XLim(1),(max(ax3.YLim)-min(ax3.YLim))*0.925+min(ax3.YLim),...
        ' b ','edgecolor','k','fontsize',fontsize,'linewidth',1,'backgroundcolor','w','fontweight','bold'); 
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
clear f; files(1).lon = NaN; files(1).lat = NaN; files(1).date = NaN;
% loop through files to grab terminus positions
for j=1:length(files)
    if ~contains(files(j).name,'CL')
        file = shaperead(files(j).name);
        % loop through file
        for k=1:length(file)
            files(length(files)+1).lon = file(k).X; % Lon
            files(length(files)).lat = file(k).Y; % Lat
            % convert to polar stereographic coordinates
            [files(length(files)).X,files(length(files)).Y] = wgs2ps(files(length(files)).lon,files(length(files)).lat,'StandardParallel',-71,'StandardMeridian',0);
            files(length(files)).x = x(dsearchn([cl.X' cl.Y'],[nanmean(files(length(files)).X) nanmean(files(length(files)).Y)])); % x
        end
    end 
end
% plot
col1 = parula(length(files)); % color scheme for plotting
ax4=axes('position',[0.08 0.08 0.23 0.4]); hold on; grid on;      
    set(gca,'fontname',fontname,'fontsize',fontsize);
    xlabel('Easting [km]'); ylabel('Northing [km]');
    % plot landsat image
    colormap(ax4,'gray');
    imagesc(LSB.x/10^3,LSB.y/10^3,flipud(LSB.im)); 
    ax4.XLim = [-2.4364e3 -2.4164e3]; ax4.YLim = [1.3002e3 1.3201e3];        
    % plot centerline
    l4=line(cl.X(2:8)/10^3,cl.Y(2:8)/10^3,'color','k','linewidth',linewidth);
    % arrow
    line2arrow(l4,'color','k','linewidth',linewidth,'headwidth',20,'headlength',20);        
    % plot terminus positions
    for j=1:length(files)
        plot(files(j).X/10^3,files(j).Y/10^3,'color',col1(j,:),'linewidth',linewidth);
    end
    text((ax4.XLim(2)-ax4.XLim(1))*0.88+ax4.XLim(1),(max(ax4.YLim)-min(ax4.YLim))*0.925+min(ax4.YLim),...
        ' c ','edgecolor','k','fontsize',fontsize,'linewidth',1,'backgroundcolor','w','fontweight','bold'); 
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
clear f; files(1).lon = NaN; files(1).lat = NaN; files(1).date = NaN;
% loop through files to grab terminus positions
for j=1:length(files)
    if ~contains(files(j).name,'CL')
        file = shaperead(files(j).name);
        % loop through file
        for k=1:length(file)
            files(length(files)+1).lon = file(k).X; % Lon
            files(length(files)).lat = file(k).Y; % Lat
            % convert to polar stereographic coordinates
            [files(length(files)).X,files(length(files)).Y] = wgs2ps(files(length(files)).lon,files(length(files)).lat,'StandardParallel',-71,'StandardMeridian',0);
            files(length(files)).x = x(dsearchn([cl.X' cl.Y'],[nanmean(files(length(files)).X) nanmean(files(length(files)).Y)])); % x
        end
    end 
end
% plot
col1 = parula(length(files)); % color scheme for plotting
ax5=axes('position',[0.39 0.08 0.23 0.4]); hold on; grid on;      
    set(gca,'fontname',fontname,'fontsize',fontsize);
    xlabel('Easting [km]'); ylabel('Northing [km]');
    % plot landsat image
    colormap(ax5,'gray');
    imagesc(LSB.x/10^3,LSB.y/10^3,flipud(LSB.im)); 
    ax5.XLim = [-2.4193e3 -2.4134e3]; ax5.YLim = [1.2762e3 1.2821e3];                
    % plot centerline
    l5=line(cl.X(5:10)/10^3,cl.Y(5:10)/10^3,'color','k','linewidth',linewidth);
    % arrow
    line2arrow(l5,'color','k','linewidth',linewidth,'headwidth',20,'headlength',20);        
    % plot terminus positions
    for j=1:length(files)
        plot(files(j).X/10^3,files(j).Y/10^3,'color',col1(j,:),'linewidth',linewidth);
    end
    text((ax5.XLim(2)-ax5.XLim(1))*0.88+ax5.XLim(1),(max(ax5.YLim)-min(ax5.YLim))*0.925+min(ax5.YLim),...
        ' d ','edgecolor','k','fontsize',fontsize,'linewidth',1,'backgroundcolor','w','fontweight','bold'); 
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
clear f; files(1).lon = NaN; files(1).lat = NaN; files(1).date = NaN;
% loop through files to grab terminus positions
for j=1:length(files)
    if ~contains(files(j).name,'CL')
        file = shaperead(files(j).name);
        % loop through file
        for k=1:length(file)
            files(length(files)+1).lon = file(k).X; % Lon
            files(length(files)).lat = file(k).Y; % Lat
            % convert to polar stereographic coordinates
            [files(length(files)).X,files(length(files)).Y] = wgs2ps(files(length(files)).lon,files(length(files)).lat,'StandardParallel',-71,'StandardMeridian',0);
            files(length(files)).x = x(dsearchn([cl.X' cl.Y'],[nanmean(files(length(files)).X) nanmean(files(length(files)).Y)])); % x
        end
    end 
end
% plot
col1 = parula(length(files)); % color scheme for plotting
ax6=axes('position',[0.72 0.08 0.23 0.4]); hold on; grid on;      
    set(gca,'fontname',fontname,'fontsize',fontsize);
    xlabel('Easting [km]'); ylabel('Northing [km]');
    % plot landsat image
    colormap(ax6,'gray');
    imagesc(LSB.x/10^3,LSB.y/10^3,flipud(LSB.im)); 
    ax6.XLim=[-2.4132e3 -2.4005e3]; ax6.YLim=[1.2635e3 1.2762e3];                
    % plot centerline
    l6=line(cl.X(6:13)/10^3,cl.Y(6:13)/10^3,'color','k','linewidth',linewidth);
    % arrow
    line2arrow(l6,'color','k','linewidth',linewidth,'headwidth',20,'headlength',20);        
    % plot terminus positions
    for j=1:length(files)
        plot(files(j).X/10^3,files(j).Y/10^3,'color',col1(j,:),'linewidth',linewidth);
    end
    text((ax6.XLim(2)-ax6.XLim(1))*0.88+ax6.XLim(1),(max(ax6.YLim)-min(ax6.YLim))*0.925+min(ax6.YLim),...
        ' e ','edgecolor','k','fontsize',fontsize,'linewidth',1,'backgroundcolor','w','fontweight','bold'); 

% save figure
if save_figure
    cd([homepath,'figures/']);
    exportgraphics(gcf,'regionalTermini.png','Resolution',600);
    disp('figure 3 saved');
end

%% Width segments

save_figure = 1; % = 1 to save figure
linewidth = 2; 
fontsize = 16;
fontname = 'Arial';

% Load Landsat image, width, width segments, and glacier extent polygon
cd([homepath,'data/Imagery/']);
ls = dir('LC08*20201104_01_T2_B8.TIF');
[LS.im,R] = readgeoraster(ls.name); [LS.ny,LS.nx] = size(LS.im);
% Polar stereographic coordinates of image boundaries
LS.x = linspace(min(R.XWorldLimits),max(R.XWorldLimits),LS.nx);
LS.y = linspace(min(R.YWorldLimits),max(R.YWorldLimits),LS.ny);
cd([homepath,'inputs-outputs/']);
cl.X = load('Crane_centerline.mat').x; cl.Y = load('Crane_centerline.mat').y;
extx = load('calculatedWidth.mat').width.extx;
exty = load('calculatedWidth.mat').width.exty;
W = load('calculatedWidth.mat').width.W;
ol = load('glacierOutline.mat').ol;
        
% Plot
col = flipud(cmocean('ice',5)); % color scheme for potting
figure(4); clf
set(gcf,'units','pixels','position',[200 200 1000 800],'defaultAxesColorOrder',[[1 1 1];[0 0 0]]);
ax1 = axes('position',[0.08 0.1 0.6 0.85]);
    hold on; imagesc(LS.x/10^3,LS.y/10^3,flipud(LS.im)); colormap("gray");
    set(gca,'fontsize',fontsize,'linewidth',linewidth); 
    xlabel('Easting [km]'); ylabel('Northing [km]'); 
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
            ' a ','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1,'fontweight','bold'); 
ax2 = axes('position',[0.73 0.2 0.2 0.65]);
    hold on; set(gca,'fontsize',fontsize,'linewidth',linewidth,'YTick',[],'XDir','reverse'); 
    xlabel('Width [km]'); yyaxis right; ylabel('Distance Along Centerline [km]'); 
    plot(W/10^3,cl.xi/10^3,'-k','linewidth',linewidth); grid on; 
    text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.95+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.05+min(get(gca,'YLim')),...
            ' b ','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1,'fontweight','bold'); 
        
% save figure
if save_figure
    cd([homepath,'figures/']);
    exportgraphics(gcf,'widthSegments.png','Resolution',600);
    disp('figure 4 saved');
end

%% Delineated Bed Elevation Profile

save_figure = 1; % = 1 to save figure
fontsize = 16; % font size for plot

% load raw picks
rawX = load('OIBpicks_raw.mat').RawX;
rawY = load('OIBpicks_raw.mat').RawY;
mdata_WGS84 = load('OIBpicks_raw.mat').mdata_WGS84;

% load radar echogram (using workflow from main_L1B.m)
fn = 'IRMCR1B_20181016_01_005.nc';
mdata = load_L1B(fn);
param.ylims_bins = [-inf inf];
good_bins = round(max(1,min(mdata.Surface)+param.ylims_bins(1)) : min(max(mdata.Surface)+param.ylims_bins(2),size(mdata.Data,1)));
% Apply AGC 
isGain = 1;
if isGain
    mdata.gainData = AGCgain(mdata.Data(good_bins,:),50,2);
end
% Elevation Correction
param = [];
if isGain
param.gain = true;
end
param.update_surf = false;
param.filter_surf = true;
param.er_ice = 3.15; % electric permitivity
param.depth = '[min(Surface_Elev)-2000 max(Surface_Elev)+25]';
% param.depth = '[min(Surface_Elev) max(Surface_Elev)]';
[mdata_WGS84,depth_good_idxs] = elevation_compensation(mdata,param);
% Compute Distance Along Fly Line
[mdata_WGS84.Easting,mdata_WGS84.Northing] = ll2utm(mdata_WGS84.Latitude, mdata_WGS84.Longitude);
UTMe = mdata_WGS84.Easting;
UTMn = mdata_WGS84.Northing;
Elvis = mdata_WGS84.Elevation;
% Euclidian Distance
for jj = 1:length(UTMe)-1
    D(jj) = sqrt( (UTMe(jj+1)-UTMe(jj))^2 + (UTMn(jj+1)-UTMn(jj))^2 + (Elvis(jj+1)-Elvis(jj))^2);
end
D = [0,D];
GPSlength = cumsum(D);
% Store in Data Structure
mdata_WGS84.Distance = GPSlength;
[mdata_WGS84.Easting,mdata_WGS84.Northing] = ...
    wgs2ps(mdata_WGS84.Longitude,mdata_WGS84.Latitude,...
    'StandardParallel',-71,'StandardMeridian',0);
clear('UTMe','UTMn','Elvis','D')
% plot vs. WGS84 elevation
figure(3); clf;
set(gcf,'position',[100 100 1000 700]);
imagesc(mdata_WGS84.Distance/1000,mdata_WGS84.Elevation_Fasttime(depth_good_idxs),mdata_WGS84.Data(depth_good_idxs,:));
    xlabel('Distance [km]');
    ylabel('WGS-84 Elevation [m]')
    set(gca,'YDir','normal')
    set(gca,'fontsize',18,'fontname','Arial')
    colormap('bone')
    hold on
    hold off;
    isContrast = 1;
    if isContrast
        keyboard
        % Run the Following line in Debug mode, then press F5 to continue...
        hcontrast = imcontrast;
    end

% Plot the Surface and Bed Picks
figure(3); hold on;
plot(mdata_WGS84.Distance/1000,mdata_WGS84.Surface_Elev,'-b','linewidth',2,'displayname','Surface Picks');
plot(rawX,rawY,'-','color',[230,245,152]./255,'linewidth',1.5,'displayname','Bed Picks');
legend('Location','southwest'); 

% Add inset plot: Landsat image, flight path, centerline
cd([homepath,'data/Imagery/']);
landsat = dir('LC08_L1GT_218106_20201015_20201104_01_T2_B8.TIF');
    [LS.im,LS.R] = readgeoraster(landsat.name); [LS.ny,LS.nx] = size(LS.im);
    % polar stereographic coordinates of image boundaries
    LS.x = linspace(min(LS.R.XWorldLimits),max(LS.R.XWorldLimits),LS.nx); 
    LS.y = linspace(min(LS.R.YWorldLimits),max(LS.R.YWorldLimits),LS.ny);
ax2 = axes('pos',[0.62 0.14 0.25 0.33],'fontname','Arial','fontsize',fontsize,...
    'XTickLabels',[],'YTickLabels',[],'linewidth',2);
    hold on; 
    imagesc(LS.x/1e3,LS.y/1e3,flipud(LS.im*1.1)); colormap('gray');
    plot(cl.X/1e3,cl.Y/1e3,'-.m','linewidth',2,'displayname','Centerline');
    plot(mdata_WGS84.Easting/1e3,mdata_WGS84.Northing/1e3,'-b','linewidth',2,'displayname','OIB flight path');
    legend('Location','east');
    xlim([-2.4293 -2.3574].*1e3); ylim([1.2119 1.2812].*1e3);

if save_figure
    cd([homepath,'figures/']);
    exportgraphics(gcf,'bedPicks.png','Resolution',600);
    disp('figure 3 saved');
end

%% Cumulative Strain Rates & Rate Factor A

close all; 

save_figure = 1; 
fontsize = 16;
fontname = 'Arial';
linewidth = 2;

% Load A, A_adj, and eta_dot_cum
cd([homepath,'inputs-outputs']);
A = load('rateFactorA.mat').A; % rate factor
A_adj = load('adjustedRateFactor.mat').A_adj; % adjusted rate factor
eta_dot_cum = load('adjustedRateFactor.mat').eta_dot_cum; % annual cumulative strain rates

% plot
col1 = parula(length(eta_dot_cum(:,1))); % color scheme for plotting

figure(5); clf; 
set(gcf,'Position',[100 100 1000 700],'defaultAxesColorOrder',[[0 0.4 0.8];[0 0 0]]);
yyaxis left; set(gca,'linewidth',2,'fontsize',fontsize,'fontname',fontname);
    xlabel('Distance Along Centerline [km]'); xlim([0 45]); ylim([-0.36 2.27]); 
    %ylabel('$$ \Sigma ( \dot{\eta} ) (s^{-1})$$','Interpreter','latex','fontname','Arial','fontsize',18);
    ylabel('Cumulative Strain Rate');
    hold on; grid on; 
    years = [2007:2011 2013:2018];
    for i=2:length(eta_dot_cum(:,1))
        plot(cl.xi(1:135)/10^3,movmean(eta_dot_cum(i,1:135),2),'-','color',col1(i,:),'linewidth',linewidth,'HandleVisibility','off');drawnow %,'displayname',num2str(years(i)));
    end
    %plot(cl.xi(1:135)/10^3,nanmean(eta_dot_cum(:,1:135),1),'-k','linewidth',linewidth+1,'displayname','mean');
yyaxis right; ylabel('Rate Factor [Pa^{-3} a^{-1}]'); ylim([1.03e-17 1.75e-17]);
    colormap(parula(9));
    cb = colorbar('position',[0.2 0.4 0.02 0.3],'fontname',fontname,...
        'fontsize',fontsize-2,'Ticks',[0.0556 0.5 0.9444],'TickLabels',[{'2009'},{'2013'},{'2017'}]);
    hold on; grid on; xlim([0 45]); legend('Position',[0.15 0.75 0.1 0.08]);
    plot(cl.xi(1:135)/10^3,A(1:135)*3.1536e7,'--k','linewidth',linewidth+0.5,'displayname','A');
    plot(cl.xi(1:135)/10^3,A_adj(1:135)*3.1536e7,'-k','linewidth',linewidth+0.5,'displayname','A_{adj}');
   
% save figure
if save_figure
    figure(5);  
    cd([homepath,'figures/']);
    exportgraphics(gcf,'rateFactor.png','Resolution',600);      
    disp('figure 5 saved.'); 
end
        
%% Beta Solution

close all;

save_figure = 1;                % = 1 to save figure
fontsize = 16;                   % font size
fontname = 'Arial';             % font name
linewidth = 2.5;                % line width
markersize = 25;                % marker size

% load variables
cd([homepath,'inputs-outputs/']);
beta = load('betaSolution.mat').beta;
betax = load('betaSolution.mat').xn;
xcf = load('betaSolution.mat').xcf;
U = load('betaSolution.mat').Un;
x0 = load('flowlineModelInitialization.mat').x0;
c0 = load('flowlineModelInitialization.mat').c0;
load('SCL_results.mat'); % stress-coupling length results
U_2018 = load('centerlineSpeedsWidthAveraged_2007-2018.mat').U_widthavg(20).speed;

% plot results
figure(6); clf
set(gcf,'Position',[272 143 1000 700],'defaultAxesColorOrder',[0 0 1; 0 0 0]);
set(gca,'fontsize',fontsize,'fontname',fontname,'linewidth',2);
hold on; grid on; 
legend('position',[0.78 0.15 0.09 0.14]);
xlim([0 45]); xlabel('Distance Along Centerline [km]'); 
yyaxis left; ylabel('\beta [s^{1/m} m^{-1/m}]'); 
plot(x0(1:c0)./10^3,beta(1:c0),'-b','linewidth',2,'displayname','\beta');     
yyaxis right; ylabel('U [m a^{-1}]'); ylim([0 650]);
plot(cl.xi(1:135)./10^3,U_2018(1:135).*3.1536e7,'k','displayname','U_{obs}','linewidth',linewidth);
plot(betax(1:dsearchn(betax',x0(c0)))./10^3,U(1:dsearchn(betax',x0(c0))).*3.1536e7,'--k','displayname','U_{mod}','linewidth',2);    

if save_figure
    figure(6);
    cd([homepath,'figures/']);
    exportgraphics(gcf,'betaSolution.png','Resolution',600);  
    disp('figure 6 saved.');    
end

%% 2018 Model Misfits
    
% Note: Must rerun first section before running this section each time

save_figure = 1;       % = 1 to save figure
save_geometry = 0;     % = 1 to save modeled 2018 geometry
plotTimeSteps = 1;     % = 1 to plot geometry, speed, cf/gl positions every decade
plotMisfits = 1;       % = 1 to plot misfit with 2018 conditions
plotClimateParams = 0; % = 1 to plot climate parameters
SMB_enhance = 0;       % = 1 to increase SMR due to decreased SMB    

% Load observed conditions
% ice surface
h_obs = load('surfaceElevationObs.mat').h;
% terminus position 
clear term termx_obs termDate_obs term_obs
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

% define time stepping (s)
dt = 0.01*3.1536e7;
t_start = 0*3.1536e7;
t_end = 9*3.1536e7;

% run the flowline model
%beta0 = load('flowlineModelInitialization.mat').beta0;
load('flowlineModelInitialization.mat', 'x0','beta0','DFW0');
[x,U,h,hb,H,gl,c,xcf,dUdx,Fgl,XCF,XGL,~] = flowlineModel(homepath,plotTimeSteps,plotMisfits,plotClimateParams,dt,t_start,t_end,beta0,DFW0,0,0,0,SMB_enhance);

% save figure
if save_figure
    cd([homepath,'figures/']);
    exportgraphics(gcf,'misfits2018.png','Resolution',600);
    disp('figure 2 saved');
end

% save geometry
if save_geometry
    cd([homepath,'inputs-outputs/']);
    save('2018_modeledConditions.mat','x','U','h','hb','H','gl','c');
    disp('final geometry saved.');
end

%% Sensitivity Tests

close all;

save_figures = 0;    % = 1 to save figure
fontsize = 16;      % font size
fontname = 'Arial'; % font name
linewidth = 2;      % line width
markersize = 10;    % marker size

% load sensitivity test no change variables
% load observations/parameters
load('flowlineModelInitialization.mat')
load('2100_noChange.mat');
% modeled 2018 conditions
H18 = load('2018_modeledConditions.mat').H;
U18 = load('2018_modeledConditions.mat').U;
c18 = load('2018_modeledConditions.mat').c;
gl18 = load('2018_modeledConditions.mat').gl;
h18 = load('2018_modeledConditions.mat').h;
hb18 = load('2018_modeledConditions.mat').hb;
x18 = load('2018_modeledConditions.mat').x;
% mean SMB over time (unperturbed)
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

    % model thickness in 2018 using surface and bed
    % calculate the thickness required to remain grounded at each grid cell
    rho_sw = 1000; % density of sea water (kg/m^3)
    rho_i = 917; % density of ice (kg/m^3)
    Hf = -(rho_sw./rho_i).*hb0; % flotation thickness (m)
    h_2018 = interp1(h(36).x,h(36).surface,x0);
    H_2018 = h_2018 - hb0;
    % find the location of the grounding line and use a floating
    % geometry from the grounding linU_widthavge to the calving front
    if length(Hf)>=find(Hf-H_2018>0,1,'first')+1
        xgl = interp1(Hf(find(Hf-H_2018>0,1,'first')-1:find(Hf-H_2018>0,1,'first')+1)...
            -H_2018(find(Hf-H_2018>0,1,'first')-1:find(Hf-H_2018>0,1,'first')+1),...
            x0(find(Hf-H_2018>0,1,'first')-1:find(Hf-H_2018>0,1,'first')+1),0,'linear','extrap'); % (m along centerline)
    else
        xgl = x0(find(Hf-H_2018>0,1,'first')-1);
    end
    gl_2018 = dsearchn(x0',xgl);
    h_2018(gl_2018+1:c0) = (1-rho_i/rho_sw).*H_2018(gl_2018+1:c0); %adjust the surface elevation of ungrounded ice to account for buoyancy
    H_2018(h_2018<0)=0-hb0(h_2018<0); h_2018(h_2018<0)=0; % surface cannot go below sea level
    h_2018(h_2018-H_2018<hb0) = hb0(h_2018-H_2018<hb0)+H_2018(h_2018-H_2018<hb0); % thickness cannot go beneath bed elevation
    
    % Calving front positions
    termdate = load('LarsenB_centerline.mat').centerline.termdate;
    termX = load('LarsenB_centerline.mat').centerline.termx;
    termY = load('LarsenB_centerline.mat').centerline.termy;
    termx = cl.x(dsearchn([cl.X' cl.Y'],[termX' termY']));

% define color schemes
% figures 5 & 6
col1=cmocean('thermal',14);  
col1(1,:) = []; col1(end-1:end,:)=[]; % don't use end member colors
% figure 8
colF = [197,27,125;222,119,174;241,182,218; 127,188,65;77,146,33]./255;

% set up axes
loop=1; % loop to minimize plotting commands
while loop==1
    figure(5); clf 
    set(gcf,'Position',[-1700 100 1000 800]);
    % (1) 
    %   (a) SMB
    axA=axes('position',[0.07 0.64 0.25 0.25]); hold on; grid on; % geometry
        set(gca,'fontsize',fontsize,'YTick',-1200:400:1200,'XTickLabel',[]);
        ylabel('Elevation [m]');         
        xlim([25 95]); ylim([-600 500]);
        plot(x1./10^3,hb1,'-k','linewidth',2);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.92+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.89+min(get(gca,'YLim')),...
            'a','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1,'fontweight','bold'); 
        % add colorbar
        colormap(col1);
        cb2=colorbar('horiz'); set(cb2,'fontname',fontname,'fontsize',fontsize-3,...
            'Ticks',0:0.5:1,'TickLabels',string(0:-5:-10),'position',[0.07 0.9 0.25 0.02]);
        set(get(cb2,'label'),'String','\DeltaSMB_{max} [m a^{-1}]','fontsize',fontsize-3);   
    axD=axes('position',[0.07 0.36 0.25 0.25]); hold on; grid on; % speed
        set(gca,'fontsize',fontsize,'YTick',0:300:1200,'XTickLabel',[]);        
        ylabel('Speed [m a^{-1}]');  
        xlim([25 95]); ylim([0 1050]);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.92+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.89+min(get(gca,'YLim')),...
            'd','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1,'fontweight','bold'); 
    axG=axes('position',[0.07 0.08 0.25 0.25]); hold on; grid on; % xgl/xcf
        set(gca,'fontsize',fontsize,'YTick',-9:3:9);
        ylabel('SMB_{mean} [m a^{-1}]'); xlabel('Distance Along Centerline [km]');
        xlim([25 95]); ylim([smb_mean(end)*3.1536e7-1 smb_mean(1)*3.1536e7+2]);              
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.92+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.89+min(get(gca,'YLim')),...
            'g','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1,'fontweight','bold');     
    %   (b) F_T   
    axB=axes('position',[0.4 0.64 0.25 0.25]); hold on; grid on; % geometry
        set(gca,'fontsize',fontsize,'YTick',-1500:500:1500,'XTickLabel',[]);
        xlim([25 95]); ylim([-600 500]);
        plot(x1./10^3,hb1,'-k','linewidth',2);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.92+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.89+min(get(gca,'YLim')),...
            'b','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1,'fontweight','bold'); 
        % add colorbar
        colormap(col1);
        cb3=colorbar('horiz'); set(cb3,'fontname',fontname,'fontsize',fontsize-3,...
            'Ticks',0:0.5:1,'TickLabels',string(0:0.5:1),'position',[0.4 0.9 0.25 0.02]);
        set(get(cb3,'label'),'String','\DeltaF_T [^oC]','fontsize',fontsize-3);
    axE=axes('position',[0.4 0.36 0.25 0.25]); hold on; grid on; % speed
        set(gca,'fontsize',fontsize,'YTick',0:300:1200,'XTickLabel',[]);        
        xlim([25 90]); ylim([0 1050]);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.92+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.89+min(get(gca,'YLim')),...
            'e','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1,'fontweight','bold'); 
    axH=axes('position',[0.4 0.08 0.25 0.25]); hold on; grid on; % gl/cf
         set(gca,'fontsize',fontsize,'YTick',0.25:0.25:1.5);
        ylabel('F_T [^oC]'); xlabel('Distance Along Centerline [km]');
        xlim([25 95]); ylim([0.1 1.3]);               
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.92+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.89+min(get(gca,'YLim')),...
            'h','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1,'fontweight','bold');         
    %   (c) DFW
    axC=axes('position',[0.72 0.64 0.25 0.25]); hold on; grid on; % geometry
        set(gca,'fontsize',fontsize,'YTick',-1500:500:1500,'XTickLabel',[]);
        xlim([25 95]); ylim([-600 500]);
        plot(x1./10^3,hb1,'-k','linewidth',2);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.92+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.89+min(get(gca,'YLim')),...
            'c','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1,'fontweight','bold'); 
         % add colorbar
        colormap(col1);
        cb4=colorbar('horiz'); set(cb4,'fontname',fontname,'fontsize',fontsize-3,...
            'Ticks',0:0.5:1,'TickLabels',string(0:5:10),'position',[0.72 0.9 0.25 0.02]);
        set(get(cb4,'label'),'String','\Deltad_{fw} [m]','fontsize',fontsize-3);
    axF=axes('position',[0.72 0.36 0.25 0.25]); hold on; grid on; % speed
        set(gca,'fontsize',fontsize,'YTick',0:300:1200,'XTickLabel',[]);        
        xlim([25 95]); ylim([0 1050]);
        % add text label             
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.92+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.89+min(get(gca,'YLim')),...
            'f','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1,'fontweight','bold'); 
    axI=axes('position',[0.72 0.08 0.25 0.25]); hold on; grid on; % xgl/xcf
        set(gca,'fontsize',fontsize);
        ylabel('d_{fw} [m]'); xlabel('Distance Along Centerline [km]');
        xlim([25 95]); ylim([DFW0-1 DFW0+11]);              
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.92+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.89+min(get(gca,'YLim')),...
            'i','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1,'fontweight','bold'); 
    % axes for other perturbations
    figure(6); clf 
    set(gcf,'Position',[-1600 200 800 800]);
    % (2) SMB_enhanced  
    axAA=axes('position',[0.1 0.64 0.38 0.25]); hold on; grid on; % geometry
        set(gca,'fontsize',fontsize,'YTick',-1500:500:1500,'XTickLabel',[]);
        ylabel('Elevation [m]'); 
        xlim([25 95]); ylim([-600 500]);
        plot(x1./10^3,hb1,'-k','linewidth',2);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.92+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.89+min(get(gca,'YLim')),...
            'a','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1,'fontweight','bold'); 
        % add colorbar
        colormap(col1);
        cb5=colorbar('horiz'); set(cb5,'fontname',fontname,'fontsize',fontsize-3,...
            'Ticks',0:0.5:1,'TickLabels',string(0:5:10),'position',[0.1 0.9 0.38 0.02]);
        set(get(cb5,'label'),'String','\DeltaSMB_{max,enh} [m a^{-1}]','fontsize',fontsize-3);      
    axCC=axes('position',[0.1 0.36 0.38 0.25]); hold on; grid on; % speed
        set(gca,'fontsize',fontsize,'YTick',0:300:1200,'XTickLabel',[]);        
        ylabel('Speed [m a^{-1}]');  
        xlim([25 95]); ylim([0 1050]);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.92+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.89+min(get(gca,'YLim')),...
            'c','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1,'fontweight','bold'); 
    axEE=axes('position',[0.1 0.08 0.38 0.25]); hold on; grid on; % gl/cf
         set(gca,'fontsize',fontsize);
        ylabel('SMB_{mean} [m a^{-1}]'); xlabel('Distance Along Centerline [km]');
        xlim([25 95]); ylim([smb_mean(end)*3.1536e7-1 smb_mean(1)*3.1536e7+1]);               
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.92+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.89+min(get(gca,'YLim')),...
            'e','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1,'fontweight','bold');         
    % (3) SMB_enh + F_T
    axBB=axes('position',[0.58 0.64 0.38 0.25]); hold on; grid on; % geometry
        set(gca,'fontsize',fontsize,'YTick',-1500:500:1500,'XTickLabel',[]);
        xlim([25 95]); ylim([-600 500]);
        plot(x1./10^3,hb1,'-k','linewidth',2);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.92+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.89+min(get(gca,'YLim')),...
            'b','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1,'fontweight','bold'); 
        % add colorbar
        colormap(col1);
        cb6=colorbar('horiz'); set(cb6,'fontname',fontname,'fontsize',fontsize-3,'Ticks',0:1/2:1,...
            'TickLabels',[{'0 & 0'},{'-5 & +0.5'},{'-10 & +1.0'}],'position',[0.58 0.9 0.38 0.02]);
        set(get(cb6,'label'),'String','\DeltaSMB_{max,enh} (m a^{-1}) & \DeltaF_T [^oC]','fontsize',fontsize-3);        
    axDD=axes('position',[0.58 0.36 0.38 0.25]); hold on; grid on; % speed
        set(gca,'fontsize',fontsize,'YTick',0:500:1600,'XTickLabel',[]);        
        xlim([25 95]); ylim([200 1600]);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.92+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.89+min(get(gca,'YLim')),...
            'd','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1); 
    axFF=axes('position',[0.58 0.08 0.38 0.25]); hold on; grid on; % xgl/xcf
        set(gca,'fontsize',fontsize);
        xlabel('Distance Along Centerline [km]'); ylabel('F_T');
        xlim([25 95]); ylim([-0.1 1.1]);              
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.92+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.89+min(get(gca,'YLim')),...
            'f','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1,'fontweight','bold'); 
    % ice mass discharge across grounding line
    figure(7); clf;
    set(gcf,'position',[-1550 300 1000 700]); 
    ax10 = axes('position',[0.08 0.57 0.27 0.39]); hold on; % SMB
        set(gca,'fontsize',fontsize); grid on;
        xlabel('SMB_{mean} [m a^{-1}]'); ylabel('Q_{gl} [Gt a^{-1}]');
        xlim([smb_mean(end)*3.1536e7-0.5 smb_mean(1)*3.1536e7+0.5]); ylim([0.8 1.8]);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.05+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.91+min(get(gca,'YLim')),...
            'a','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1,'fontweight','bold'); 
        % add colorbar
        colormap(col1);
        cb5=colorbar('northoutside','fontname',fontname,'fontsize',fontsize-2,'Ticks',0:0.5:1,'TickLabels',string(0:-5:-10));
        set(get(cb5,'label'),'String','\DeltaSMB_{max} [m a^{-1}]','fontsize',fontsize-2);
    ax11 = axes('position',[0.39 0.57 0.27 0.39]); hold on; % F_T
        set(gca,'fontsize',fontsize); grid on; 
        xlabel('F_T [^oC]');  
        xlim([0.1 1.3]); ylim([0.8 1.8]); 
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.05+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.91+min(get(gca,'YLim')),...
            'b','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1,'fontweight','bold'); 
        % add colorbar
        colormap(col1); 
        cb4=colorbar('northoutside','fontname',fontname,'fontsize',fontsize-2,'Ticks',0:0.5:1,'TickLabels',string(0:0.5:1));
        set(get(cb4,'label'),'String','\DeltaF_T [^oC]','fontsize',fontsize-2);
    ax12 = axes('position',[0.7 0.57 0.27 0.39]); hold on; % d_fw
        set(gca,'fontsize',fontsize); grid on;
        xlabel('d_{fw} [m]'); 
        xlim([DFW0-0.5 DFW0+10.5]); ylim([0.8 1.8]);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.05+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.91+min(get(gca,'YLim')),...
            'c','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1,'fontweight','bold'); 
         % add colorbar
        colormap(col1);
        cb6=colorbar('northoutside','fontname',fontname,'fontsize',fontsize-2,'Ticks',0:0.5:1,'TickLabels',string(0:5:10)); 
        set(get(cb6,'label'),'String','\Deltad_{fw} [m]','fontsize',fontsize-2);
    ax13 = axes('position',[0.17 0.1 0.3 0.35]); hold on; 
        set(gca,'fontsize',fontsize); grid on; 
        xlabel('SMB_{mean} [m a^{-1}]'); ylabel('Q_{gl} [Gt a^{-1}]'); 
        xlim([smb_mean(end)*3.1536e7-0.5 smb_mean(1)*3.1536e7+0.5]); ylim([0.8 1.8]); 
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.05+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.91+min(get(gca,'YLim')),...
            'd','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1,'fontweight','bold'); 
        % add colorbar
        colormap(col1); 
        cb4=colorbar('eastoutside','fontname',fontname,'fontsize',fontsize-2,'Ticks',0:0.5:1,'TickLabels',string(0:-5:-10));
        set(get(cb4,'label'),'String','\DeltaSMB_{max,enh} [m a^{-1}]','fontsize',fontsize-2);
    ax14 = axes('position',[0.58 0.1 0.3 0.35]); hold on;
        set(gca,'fontsize',fontsize); grid on;
        xlabel('\DeltaF_T [^oC]');
        xlim([-0.1 1.1]); ylim([0.8 1.8]);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.05+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.91+min(get(gca,'YLim')),...
            'e','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1,'fontweight','bold'); 
        % add colorbar
        colormap(col1);
        cb5=colorbar('eastoutside','fontname',fontname,'fontsize',fontsize-2,'Ticks',0:0.5:1,'TickLabels',string(0:-5:-10));
        set(get(cb5,'label'),'String','\DeltaSMB_{max,enh} [m a^{-1}]','fontsize',fontsize-2);
    % calving front position and discharge over time
    figure(8); clf;  
    set(gcf,'Position',[-1300 200 1000 800]); 
    ax16 = axes('position',[0.06 0.57 0.7 0.4]); hold on;
        set(gca,'fontsize',fontsize); grid on;
        xlim([1995 2100]); ylim([25 95]);
        ylabel('Calving Front Position [km]');
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.02+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.9+min(get(gca,'YLim')),...
            'a','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1,'fontweight','bold');
        % add grey rectangle for Larsen B ice shelf collapse
        clps_st = 2002.0849; % Jan 31, 2002
        clps_end = 2002.2822; % April 13, 2002
        rectangle(ax16,'Position',[clps_st 25 clps_end-clps_st 90],'FaceColor',[189,189,189]./255,'EdgeColor',[189,189,189]./255);
    ax17 = axes('position',[0.06 0.09 0.7 0.4]); hold on; 
        set(gca,'fontsize',fontsize); grid on; 
        leg = legend('position',[0.8 0.4 0.17 0.24]);
        xlim([1995 2100]); ylim([0 7.5]);        
        xlabel('Year'); ylabel('Grounding Line Discharge [Gt a^{-1}]'); ylim([0 7]);
        % add text label            
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.02+min(get(gca,'XLim')),...
            (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.9+min(get(gca,'YLim')),...
            'b','backgroundcolor','w','fontsize',fontsize,'linewidth',linewidth-1,'fontweight','bold');
        % add grey rectangle for Larsen B ice shelf collapse
        rectangle(ax17,'Position',[clps_st 0 clps_end-clps_st 7],'FaceColor',[189,189,189]./255,'EdgeColor',[189,189,189]./255);
    loop=loop+1; % exit loop        
end

% loop through results, split into each scenario
cd([homepath,'scripts/3_sensitivity_tests/results/1_SMB_DFW_TF/']);
files = dir('*geom.mat');
for i=1:length(files)
    if contains(files(i).name,'SMB0') && contains(files(i).name,'DFW0') && contains(files(i).name,'TF0_')
        files(i).change = 0;
        files(i).changeIn = 'SMB';
        files(i).name = files(i).name; 
        files(length(files)+1).change = 0;
        files(length(files)).changeIn = 'DFW';
        files(length(files)).name = files(i).name;
        files(length(files)+1).change = 0; 
        files(length(files)).changeIn = 'TF';
        files(length(files)).name = files(i).name;        
    % (a) SMB
    elseif contains(files(i).name,'DFW0m') && contains(files(i).name,'TF0_')
        files(i).change = str2double(files(i).name(regexp(files(i).name,'B')+1:...
            regexp(files(i).name,'_D')-1)); 
        files(i).changeIn = 'SMB';  
        files(i).name = files(i).name;
    % (b) DFW
    elseif contains(files(i).name,'SMB0') && contains (files(i).name,'TF0_')
        files(i).change = str2double(files(i).name(regexp(files(i).name,'W')+1:...
            regexp(files(i).name,'m_')-1))-DFW0; % m 
        files(i).changeIn = 'DFW';  
        files(i).name = files(i).name;
    % (c) F_T
    elseif contains(files(i).name,'SMB0') && contains (files(i).name,'DFW0m')        
        files(i).change = str2double(files(i).name(regexp(files(i).name,'TF')+2:...
            regexp(files(i).name,'_geom')-1)); % m 
        files(i).changeIn = 'TF';  
        files(i).name = files(i).name;        
    end
end
    
% sort by changeIn and change
files = struct2table(files);
files = sortrows(files,[8,7]);
% save indices for SMB & SMR
ISMB = flipud(find(strcmp('SMB',table2array(files(:,8)))));
IDFW = find(strcmp('DFW',table2array(files(:,8))));
ITF = find(strcmp('TF',table2array(files(:,8))));
files = table2struct(files);

% (1) SMB
dSMB = zeros(length(ISMB),1); % change in SMB
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
    Fsmb(i) = (H2(gl2)*U2(gl2)*W(gl2))*917*1e-12*3.1536e7; % Gt/a
    figure(5); hold on;
        % ice surface
        plot(axA,x2(1:c2)./10^3,h2(1:c2),'color',col1(i,:),'linewidth',linewidth-0.5);
        % calving front
        plot(axA,[x2(c2) x2(c2)]./10^3,[h2(c2) h2(c2)-H2(c2)],'color',col1(i,:),'linewidth',linewidth-0.5);
        % floating bed
        plot(axA,x2(gl2:c2)./10^3,h2(gl2:c2)-H2(gl2:c2),'color',col1(i,:),'linewidth',linewidth-0.5);
        % ice surface speed
        plot(axD,x2(1:c2)./10^3,U2(1:c2).*3.1536e7,'color',col1(i,:),'linewidth',linewidth-0.5);
        % dH
        dHgl(i) = H2(gl2)-H1(gl1); 
        % cf and gl positions
        plot(axG,x2(gl2)/1e3,smb_mean(i)*3.1536e7,'x','markersize',10,'linewidth',linewidth,'color',col1(i,:));
        plot(axG,x2(c2)/1e3,smb_mean(i)*3.1536e7,'o','markersize',10,'linewidth',linewidth,'color',col1(i,:));        
    % ice mass discharge
    plot(ax10,smb_mean(i)*3.1536e7,Fsmb(i),'o','markersize',markersize,'linewidth',linewidth,'color',col1(i,:));
    if i==1
        figure(8);
        plot(ax16,termdate,termx/10^3,'xk','markersize',markersize-1,'linewidth',linewidth,'displayName','observed');        
        plot(ax16,t/3.1536e7+2009,movmean(XCF2/10^3,100),'-k','linewidth',linewidth+0.5,'displayName','unperturbed');
        plot(ax17,F_obs(:,1),F_obs(:,2),'xk','markersize',markersize-1,'linewidth',linewidth,'displayName','observed');                
        plot(ax17,t/3.1536e7+2009,movmean(Fgl2,100),'-k','linewidth',linewidth+0.5,'displayName','unperturbed');         
    elseif i==length(ISMB)
        figure(8);
        plot(ax16,t/3.1536e7+2009,movmean(XCF2/10^3,100),'-','color',colF(1,:),'linewidth',linewidth+0.5,'displayName','\DeltaSMB_{max}');        
        plot(ax17,t/3.1536e7+2009,movmean(Fgl2,100),'-','color',colF(1,:),'linewidth',linewidth+0.5,'displayName','\DeltaSMB_{max}'); 
    end
    % save results 
    dSMB(i) = files(ISMB(i)).change; % m/a
    dUgl(i) = (U2(gl2)-U1(gl1)).*3.1536e7; % m/a     
    dL(i) = x2(c2)-x1(c1); % m   
    dgl(i) = x2(gl2)-x1(gl1); % m
    Ugl(i) = U2(gl2)*3.1536e7; % m/a
end
% create table to store result quantities
varNames = {'dsmb','dL','dxgl','dHgl','dUgl','Qsmb'};
T_SMB = table(round(dSMB),round(dL)/10^3,round(dgl)/10^3,round(dHgl),round(dUgl),Fsmb,'VariableNames',varNames);

% (2) F_T
dTF = zeros(length(ITF),1); % change in F_T
dsmr = zeros(length(ITF),1); % change in smr due to F_T
dHgl = zeros(length(ITF),1); % change in mean thickness
dUgl = zeros(length(ITF),1); % change in mean ice speed
dL = zeros(length(ITF),1); % change in length
dgl = zeros(length(ITF),1); % change in grounding line position
Ugl = zeros(length(ITF),1); % speed at grounding line 
F = zeros(length(ITF),1); % grounding line discharge
% define thermal forcing
TF0 = 0.2; % ^oC - estimated from Larsen B icebergs
for i=1:length(ITF)
    load(files(ITF(i)).name);
    W = interp1(x0,W0,x2); % interpolate width on spatial grid
    % F (Gt/a) = (U m/s)*(H m)*(W m)*(917 kg/m^3)*(3.1536e7 s/a)*(1e-12 Gt/kg)
    % use the mean of a rectangular and an ellipsoidal bed
    F(i) = (H2(gl2)*U2(gl2)*W(gl2))*917*1e-12*3.1536e7; % Gt/a
    % extract change in F_T
    dTF(i) = files(ITF(i)).change; % m/a
    figure(5); hold on;
        % ice surface
        plot(axB,x2(1:c2)./10^3,h2(1:c2),'color',col1(i,:),'linewidth',linewidth-0.5);
        % calving front
        plot(axB,[x2(c2) x2(c2)]./10^3,[h2(c2) h2(c2)-H2(c2)],'color',col1(i,:),'linewidth',linewidth-0.5);
        % floating bed
        plot(axB,x2(gl2:c2)./10^3,h2(gl2:c2)-H2(gl2:c2),'color',col1(i,:),'linewidth',linewidth-0.5);
        % ice surface speed
        plot(axE,x2(1:c2)./10^3,U2(1:c2).*3.1536e7,'color',col1(i,:),'linewidth',linewidth-0.5);
        % dH
        dHgl(i) = H2(gl2)-H1(gl1); 
        % cf and gl positions
        plot(axH,x2(gl2)/1e3,TF0+dTF(i),'x','markersize',markersize,...
            'linewidth',linewidth,'color',col1(i,:));
        plot(axH,x2(c2)/1e3,TF0+dTF(i),'o','markersize',markersize,...
            'linewidth',linewidth,'color',col1(i,:));   
    % calving front position and ice mass discharge
    figure(7);
    plot(ax11,TF0+dTF(i),F(i),'o','markersize',markersize,'linewidth',linewidth,'color',col1(i,:));
    if i==length(ITF)
        figure(8);
        plot(ax16,t/3.1536e7+2009,movmean(XCF2/10^3,100),'-','color',colF(2,:),'linewidth',linewidth+0.5,'displayName','\DeltaF_T');        
        plot(ax17,t/3.1536e7+2009,movmean(Fgl2,100),'-','color',colF(2,:),'linewidth',linewidth+0.5,'displayName','\DeltaF_{T,max}');
    end
    % save results for data table
    %dSMR(i) = dSMR_max*3.1536e7;
    dUgl(i) = (U2(gl2)-U1(gl1)).*3.1536e7; % m/a     
    dL(i) = x2(c2)-x1(c1); % m   
    dgl(i) = x2(gl2)-x1(gl1); % m
    Ugl(i) = U2(gl2)*3.1536e7; % m/a
end
varNames = {'dTF','dL','dgl','dHgl','dUgl','Qgl'};
T_TF = table(dTF,round(dL)/10^3,round(dgl)/10^3,round(dHgl),round(dUgl),F,'VariableNames',varNames);

% (3) DFW
dDFW = zeros(length(IDFW),1); % change in DFW
dHgl = zeros(length(IDFW),1); % change in mean thickness
dUgl = zeros(length(IDFW),1); % change in mean ice speed
dL = zeros(length(IDFW),1); % change in length
dgl = zeros(length(IDFW),1); % change in grounding line position
Ugl = zeros(length(IDFW),1); % speed at grounding line 
Fdfw = zeros(length(IDFW),1); % grounding line discharge
col1 = cmocean('thermal',length(IDFW)+3); col1(1,:)=[];col1(end-1:end,:)=[];
for i=1:length(IDFW)
    % load file
    load(files(IDFW(i)).name);
    % calculate grounding line discharge
    W = interp1(x0,W0,x2); % interpolate width on spatial grid
    % F (Gt/a) = (U m/s)*(H m)*(W m)*(917 kg/m^3)*(3.1536e7 s/a)*(1e-12 Gt/kg)
    % use the mean of a rectangular and an ellipsoidal bed
    Fdfw(i) = (H2(gl2)*U2(gl2)*W(gl2))*917*1e-12*3.1536e7; % Gt/a    
    % plot
    figure(5); hold on;
        % ice surface
        plot(axC,x2(1:c2)./10^3,h2(1:c2),'color',col1(i,:),'linewidth',linewidth-0.5);
        % calving front
        plot(axC,[x2(c2) x2(c2)]./10^3,[h2(c2)-H2(c2) h2(c2)],'color',col1(i,:),'linewidth',linewidth-0.5);
        % floating bed
        plot(axC,x2(gl2:c2)./10^3,h2(gl2:c2)-H2(gl2:c2),'color',col1(i,:),'linewidth',linewidth-0.5);
        % ice surface speed
        plot(axF,x2(1:c2)./10^3,U2(1:c2).*3.1536e7,'color',col1(i,:),'linewidth',linewidth-0.5);
        % dH
        dHgl(i) = H2(gl2)-H1(gl1); % m
        % cf and gl positions
        plot(axI,x2(gl2)/10^3,files(IDFW(i)).change+DFW0,'x',...
            'markersize',10,'linewidth',linewidth,'color',col1(i,:));
        plot(axI,x2(c2)/10^3,files(IDFW(i)).change+DFW0,'o',...
            'markersize',10,'linewidth',linewidth,'color',col1(i,:));        
    % ice mass discharge
    figure(7); 
    plot(ax12,files(IDFW(i)).change+DFW0,Fdfw(i),'o','markersize',markersize,'linewidth',linewidth,'color',col1(i,:));    
    % calving front position and grounding line discharge over time
    if i==length(IDFW)
        figure(8);
        plot(ax16,t/3.1536e7+2009,movmean(XCF2/10^3,100),'-','color',colF(3,:),'linewidth',linewidth+0.5,'displayName','\Deltad_{fw,max}');        
        plot(ax17,t/3.1536e7+2009,movmean(Fgl2,100),'-','color',colF(3,:),'linewidth',linewidth+0.5,'displayName','\Deltad_{fw,max}'); 
    end
    % save results 
    dDFW(i) = files(IDFW(i)).change; % m/a
    dUgl(i) = (U2(gl2)-U1(gl1)).*3.1536e7; % m/a     
    dL(i) = x2(c2)-x1(c1); % m   
    dgl(i) = x2(gl2)-x1(gl1); % m
    Ugl(i) = U2(gl2)*3.1536e7; % m/a
end
% create table to store result quantities
varNames = {'dfwd','dL','dxgl','dHgl','dUgl','Qfwd'};
T_DFW = table(dDFW,round(dL)/10^3,round(dgl)/10^3,round(dHgl),round(dUgl),Fdfw,'VariableNames',varNames);

% (4) SMB_enh
cd([homepath,'scripts/3_sensitivity_tests/results/2_SMB_enh/']);
dsmb_enh = zeros(length(IDFW),1); % change in SMB_enh
dHgl = zeros(length(IDFW),1); % change in mean thickness
dUgl = zeros(length(IDFW),1); % change in mean ice speed
dL = zeros(length(IDFW),1); % change in length
dgl = zeros(length(IDFW),1); % change in grounding line position
Ugl = zeros(length(IDFW),1); % speed at grounding line 
F = zeros(length(IDFW),1); % grounding line discharge
% define thermal forcing
TF0 = 0.2; % ^oC - estimated from Larsen B icebergs
clear files; files = dir('*geom.mat');
% sort files by perturbation magnitude
for i=1:length(files)
    files(i).change = str2double(files(i).name(regexp(files(i).name,'B')+1:...
                  regexp(files(i).name,'_enh')-1)); % m/a
end
files = struct2table(files);
files = sortrows(files,7,'descend');
files = table2struct(files); 
for i=1:length(files)
    load(files(i).name);
    W = interp1(x0,W0,x2); % interpolate width on spatial grid
    % F (Gt/a) = (U m/s)*(H m)*(W m)*(917 kg/m^3)*(3.1536e7 s/a)*(1e-12 Gt/kg)
    % use the mean of a rectangular and an ellipsoidal bed
    F(i) = (H2(gl2)*U2(gl2)*W(gl2))*917*1e-12*3.1536e7; % Gt/a
    % estimate initial melt rate using Eqn from Slater et al. (2020):
    mdot0 = (3*10^-4*-hb0(gl1)*((sum(RO0(1:gl1)))*86400)^0.39 + 0.15)*TF0^1.18/86400; % m/s
    figure(6); hold on;
        % ice surface
        plot(axAA,x2(1:c2)./10^3,h2(1:c2),'color',col1(i,:),'linewidth',linewidth-0.5);
        % calving front
        plot(axAA,[x2(c2) x2(c2)]./10^3,[h2(c2) h2(c2)-H2(c2)],'color',col1(i,:),'linewidth',linewidth-0.5);
        % floating bed
        plot(axAA,x2(gl2:c2)./10^3,h2(gl2:c2)-H2(gl2:c2),'color',col1(i,:),'linewidth',linewidth-0.5);
        % ice surface speed
        plot(axCC,x2(1:c2)./10^3,U2(1:c2).*3.1536e7,'color',col1(i,:),'linewidth',linewidth-0.5);
        % dH
        dHgl(i) = H2(gl2)-H1(gl1); 
        % cf and gl positions
        plot(axEE,x2(gl2)/1e3,smb_mean(i)*3.1536e7,'x','markersize',markersize,...
            'linewidth',linewidth,'color',col1(i,:));
        plot(axEE,x2(c2)/1e3,smb_mean(i)*3.1536e7,'o','markersize',markersize,...
            'linewidth',linewidth,'color',col1(i,:));   
    % calving front position and ice mass discharge
    figure(7);
    plot(ax13,smb_mean(i)*3.1536e7,F(i),'o','markersize',markersize,'linewidth',linewidth,'color',col1(i,:));
    if i==length(files)
        figure(8);
        plot(ax16,t/3.1536e7+2009,movmean(XCF2/10^3,100),'-','color',colF(4,:),'linewidth',linewidth+0.5,'displayName','\DeltaSMB_{enh}');        
        plot(ax17,t/3.1536e7+2009,movmean(Fgl2,100),'-','color',colF(4,:),'linewidth',linewidth+0.5,'displayName','\DeltaSMB_{enh}');
    end
    % save results 
    dsmb_enh(i) = files(i).change; % m/a
    %dSMR(i) = dSMR_max*3.1536e7; % m/a
    dUgl(i) = (U2(gl2)-U1(gl1)).*3.1536e7; % m/a     
    dL(i) = x2(c2)-x1(c1); % m   
    dgl(i) = x2(gl2)-x1(gl1); % m
    Ugl(i) = U2(gl2)*3.1536e7; % m/a
end
varNames = {'dSMB_enh','dL','dgl','dHgl','dUgl','Qgl'};
T_smb_enh = table(dsmb_enh,round(dL)/10^3,round(dgl)/10^3,round(dHgl),round(dUgl),F,'VariableNames',varNames);

% (5) SMB_enh + F_T
cd([homepath,'scripts/3_sensitivity_tests/results/3_SMB_enh+TF/']);
dsmb_enh = zeros(length(IDFW),1); % change in SMB_enh
dTF = zeros(length(IDFW),1); % change in F_T
dHgl = zeros(length(IDFW),1); % change in mean thickness
dUgl = zeros(length(IDFW),1); % change in mean ice speed
dL = zeros(length(IDFW),1); % change in length
dgl = zeros(length(IDFW),1); % change in grounding line position
Ugl = zeros(length(IDFW),1); % speed at grounding line 
F = zeros(length(IDFW),1); % grounding line discharge
clear files; files = dir('*geom.mat');
% sort files by perturbation magnitude
for i=1:length(files)
    files(i).change_TF = str2double(files(i).name(regexp(files(i).name,'F')+1:...
                  regexp(files(i).name,'_geom')-1)); % ^oC
    files(i).change_SMB = str2double(files(i).name(regexp(files(i).name,'B')+1:...
                  regexp(files(i).name,'_enh')-1)); % m/a
end
files = struct2table(files);
files = sortrows(files,7,'ascend');
files = table2struct(files); 
for i=1:length(files)
    load(files(i).name);
    W = interp1(x0,W0,x2); % interpolate width on spatial grid
    % F (Gt/a) = (U m/s)*(H m)*(W m)*(917 kg/m^3)*(3.1536e7 s/a)*(1e-12 Gt/kg)
    % use the mean of a rectangular and an ellipsoidal bed
    F(i) = (H2(gl2)*U2(gl2)*W(gl2))*917*1e-12*3.1536e7; % Gt/a
    figure(6); hold on;
        % ice surface
        plot(axBB,x2(1:c2)./10^3,h2(1:c2),'color',col1(i,:),'linewidth',linewidth-0.5);
        % calving front
        plot(axBB,[x2(c2) x2(c2)]./10^3,[h2(c2) h2(c2)-H2(c2)],'color',col1(i,:),'linewidth',linewidth-0.5);
        % floating bed
        plot(axBB,x2(gl2:c2)./10^3,h2(gl2:c2)-H2(gl2:c2),'color',col1(i,:),'linewidth',linewidth-0.5);
        % ice surface speed
        plot(axDD,x2(1:c2)./10^3,U2(1:c2).*3.1536e7,'color',col1(i,:),'linewidth',linewidth-0.5); 
        % cf and gl positions
        plot(axFF,x2(gl2)/1e3,files(i).change_TF,'x',...
            'markersize',markersize,'linewidth',linewidth,'color',col1(i,:));
        plot(axFF,x2(c2)/1e3,files(i).change_TF,'o',...
            'markersize',markersize,'linewidth',linewidth,'color',col1(i,:));        
    % ice mass discharge
    figure(7);
    plot(ax14,files(i).change_TF,F(i),'o','markersize',markersize,'linewidth',linewidth,'color',col1(i,:));
    if i==length(files)
        figure(8);
        plot(ax16,t/3.1536e7+2009,movmean(XCF2/10^3,100),'-','color',colF(5,:),'linewidth',linewidth+0.5,'displayName','\DeltaSMB_{enh,max} & \DeltaF_{T,max}');                
        plot(ax17,t/3.1536e7+2009,movmean(Fgl2,100),'-','color',colF(5,:),'linewidth',linewidth+0.5,'displayName','\DeltaSMB_{enh,max} & \DeltaF_{T,max}');
    end  
    % save results 
    dsmb_enh(i) = files(i).change_SMB; % m/a
    dTF(i) = files(i).change_TF; % m/a
    % dH
    dHgl(i) = H2(gl2)-H1(gl1);
    %dSMR(i) = dSMR_max*3.1536e7; % m/a
    dUgl(i) = (U2(gl2)-U1(gl1)).*3.1536e7; % m/a     
    dL(i) = x2(c2)-x1(c1); % m   
    dgl(i) = x2(gl2)-x1(gl1); % m
    Ugl(i) = U2(gl2)*3.1536e7; % m/a
end
varNames = {'dSMB_enh','dTF','dL','dgl','dHgl','dUgl','Qgl'};
T_smb_enh_TF = table(dsmb_enh,dTF,round(dL)/10^3,round(dgl)/10^3,round(dHgl),round(dUgl),F,'VariableNames',varNames);

% Add observed 2018 conditions to figures 5 and 6
b=1; % create a loop to make section collapsible
while b==1
    % geometry
    plot(axA,x18/10^3,h18,'--k','linewidth',linewidth-0.5);
    plot(axA,x18/10^3,h18-H18,'--k','linewidth',linewidth-0.5);
    plot(axA,[x18(c18)/10^3; x18(c18)/10^3],[h18(c18); h18(c18)-H18(c18)],'--k','linewidth',linewidth-0.5);
    plot(axB,x18/10^3,h18,'--k','linewidth',linewidth-0.5);
    plot(axB,x18/10^3,h18-H18,'--k','linewidth',linewidth-0.5);
    plot(axB,[x18(c18)/10^3; x18(c18)/10^3],[h18(c18)-H18(c18); h18(c18)],'--k','linewidth',linewidth-0.5);
    plot(axC,x18/10^3,h18,'--k','linewidth',linewidth-0.5);
    plot(axC,x18/10^3,h18-H18,'--k','linewidth',linewidth-0.5);
    plot(axC,[x18(c18)/10^3; x18(c18)/10^3],[h18(c18)-H18(c18); h18(c18)],'--k','linewidth',linewidth-0.5);
    plot(axAA,x18/10^3,h18,'--k','linewidth',linewidth-0.5);
    plot(axAA,x18/10^3,h18-H18,'--k','linewidth',linewidth-0.5);
    plot(axAA,[x18(c18)/10^3; x18(c18)/10^3],[h18(c18); h18(c18)-H18(c18)],'--k','linewidth',linewidth-0.5);
    plot(axBB,x18/10^3,h18,'--k','linewidth',linewidth-0.5);
    plot(axBB,x18/10^3,h18-H18,'--k','linewidth',linewidth-0.5);
    plot(axBB,[x18(c18)/10^3; x18(c18)/10^3],[h18(c18)-H18(c18); h18(c18)],'--k','linewidth',linewidth-0.5);
    % speed
    plot(axD,x18/10^3,U18*3.1536e7,'--k','linewidth',linewidth-0.5);
    plot(axE,x18/10^3,U18*3.1536e7,'--k','linewidth',linewidth-0.5);
    plot(axF,x18/10^3,U18*3.1536e7,'--k','linewidth',linewidth-0.5);
    plot(axCC,x18/10^3,U18*3.1536e7,'--k','linewidth',linewidth-0.5);
    plot(axDD,x18/10^3,U18*3.1536e7,'--k','linewidth',linewidth-0.5);
%     % calving front and groundling line positions
%     plot(axJ,x18(c18)/10^3,smb_mean(1),'ok','linewidth',linewidth-0.5,'markersize',markersize);
%     plot(axJ,x18(gl18)/10^3,smb_mean(1),'xk','linewidth',linewidth-0.5,'markersize',markersize);
%     plot(axK,x18(c18)/10^3,TF0,'ok','linewidth',linewidth-0.5,'markersize',markersize);
%     plot(axK,x18(gl18)/10^3,TF0,'xk','linewidth',linewidth-0.5,'markersize',markersize);
%     plot(axL,x18(c18)/10^3,DFW0,'ok','linewidth',linewidth-0.5,'markersize',markersize);
%     plot(axL,x18(gl18)/10^3,DFW0,'xk','linewidth',linewidth-0.5,'markersize',markersize);
%     plot(axGG,x18(c18)/10^3,smb_mean(1),'ok','linewidth',linewidth-0.5,'markersize',markersize);
%     plot(axGG,x18(gl18)/10^3,smb_mean(1),'xk','linewidth',linewidth-0.5,'markersize',markersize);
%     plot(axHH,x18(c18)/10^3,TF0,'ok','linewidth',linewidth-0.5,'markersize',markersize);
%     plot(axHH,x18(gl18)/10^3,TF0,'xk','linewidth',linewidth-0.5,'markersize',markersize);
    
    b=b+1;
end
    
% save figures 
if save_figures
    cd([homepath,'figures/']);
    figure(5); exportgraphics(gcf,'sensitivityTests_geom+speed_independent.png','Resolution',600);
    figure(6); exportgraphics(gcf,'sensitivityTests_geom+speed_enhanced.png','Resolution',600);
    figure(7); exportgraphics(gcf,'sensitivityTests_discharge.png','Resolution',600);
    figure(8); exportgraphics(gcf,'sensitivityTests_QglXcf.png','Resolution',600);
    disp('figures 5-8 saved.');
end

