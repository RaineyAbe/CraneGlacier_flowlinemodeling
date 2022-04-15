%% Script to calculate Crane Glacier width using glacier outline and using the centerline vector to cut off the width.
% Rainey Aberle
% Fall 2020

clear all; close all;

% Define homepath
homepath = '/Users/raineyaberle/Research/MS/CraneGlacier_flowlinemodeling/';
cd([homepath,'inputs-outputs/']);

% Add path to required functions
% addpath([homepath,'functions/OIBPicking/functions/']);

% Load centerline
cl.x = load('Crane_centerline.mat').x; cl.y = load('Crane_centerline.mat').y;

outlineExists = 1;  % Does the glacier outline already exist? (0=no, 1=yes)
save_width = 1;     % = 1 to save width

%% Load and display Landsat image

cd([homepath,'data/imagery/']);

ls = dir('LC08*20200114_01_T2_B8.TIF');
[LS.im,R] = readgeoraster(ls.name); [LS.ny,LS.nx] = size(LS.im);
% Polar stereographic coordinates of image boundaries
LS.x = linspace(min(R.XWorldLimits),max(R.XWorldLimits),LS.nx);
LS.y = linspace(min(R.YWorldLimits),max(R.YWorldLimits),LS.ny);
        
% Plot Landsat image
fig1 = figure(1); clf
ax1 = subplot(1,3,[1 2]); hold on; 
    imagesc(ax1, LS.x,LS.y,flipud(LS.im)); 
    colormap("gray");
    set(gca,'fontsize',14,'linewidth',2); 
    set(gcf,'units','centimeters','position',[13 5 22 25]);
    xlabel('Easting [m]'); ylabel('Northing [m]'); 
    legend('Location','east','color',[0.6,0.6,0.6]);
    % Specify the region of interest (ROI)
    xlim([-2.43e6 -2.385e6]); ylim([1.215e6 1.282e6]); 
set(gcf,'position',[2 3 41 22]); 
col = parula(5); % color scheme for potting

% Plot centerline
plot(ax1, cl.x,cl.y,'-w','linewidth',2,'displayname','glacier centerline');
    
% Manually pick or load glacier outline
if outlineExists==0
    % Make manual picks on figure
    [ol.x,ol.y]=ginput();
        %[ol.xa,ol.ya]=ginput();
        %[ol.xb,ol.yb]=ginput();
        %ol.x = [ol.x; ol.xa; ol.xb]; ol.y = [ol.y; ol.ya; ol.yb];
        
    % Create polygon from picks, plot
    ol.polygon = polyshape(ol.x,ol.y);
    plot(ol.polygon,'displayname','glacier outline','facecolor',col(3,:));
    
    % Save picks
    cd([homepath,'inputs-outputs/']);
    save('glacier_outline_2.mat','ol');
else
    cd([homepath,'inputs-outputs/']);
    ol = load('glacier_outline_2.mat').ol; 
    % Create polygon from picks, plot
    ol.polygon = polyshape(ol.x,ol.y);
    plot(ol.polygon,'displayname','glacier outline','facecolor',col(3,:));
end

%% Draw lines perpendicular to each point along centerline
% Use a staggered grid of centerline coordinates for calculating slope
cl.xm = (cl.x(2:end)+cl.x(1:end-1))./2;
cl.ym = (cl.y(2:end)+cl.y(1:end-1))./2;

% Initialize variables
segs.l = 7000; % length of lines extending from centerline (m)
dx = 10; % distance between points of line (m)
segs.x = zeros(length(cl.x),2); % glacier width segments (m)
segs.y = zeros(length(cl.x),2);
% slopes (unitless)
m = zeros(1,length(cl.x)); 
m_rev = zeros(1,length(cl.x)); 
% width segments with higher spatial resolution (m)
segs.xn = zeros(length(cl.x),segs.l*2/dx+1);
segs.yn = zeros(length(cl.x),segs.l*2/dx+1);    
% sort by x coordinate
% [segs.xn, segs.yn] = sort([segs.xn, segs.yn],1);

for i=2:length(cl.x)
    % Calculate slope angle at each point
    if i==length(cl.x)
        m(i) = m(i-1); 
        m_rev(i) = -1/m(i);
    else
        m(i) = (cl.ym(i)-cl.ym(i-1))/(cl.xm(i)-cl.xm(i-1));
        m_rev(i) = -1/m(i);
    end 

    % Calculate distance to points  
    delta_x = sqrt((segs.l)^2/(m_rev(i)^2+1));
    delta_y = sqrt((segs.l)^2/(m_rev(i)^2+1))*m_rev(i);

    % Add dy/dx to the point on the centerline
    segs.x(i,2) = cl.x(i)+delta_x; segs.x(i,1) = cl.x(i)-delta_x;
    segs.y(i,2) = cl.y(i)+delta_y; segs.y(i,1) = cl.y(i)-delta_y;        

    % Increase number of points in segment
    segs.xn(i,:) = segs.x(i,1):(segs.x(i,2)-segs.x(i,1))*(dx/(segs.l*2)):segs.x(i,2);
    segs.yn(i,:) = segs.y(i,1):(segs.y(i,2)-segs.y(i,1))*(dx/(segs.l*2)):segs.y(i,2);

    % Plot results on figure
    figure(1);
    if i==2
        plot(ax1, segs.xn(i,:),segs.yn(i,:),'-b','linewidth',1,'displayname','perpendicular lines');            
    else
        plot(ax1, segs.xn(i,:),segs.yn(i,:),'-b','linewidth',1,'handlevisibility','off');
    end

    % Once all other points are complete, go to the first point
    % and use the 2nd point slope to calculate segment 
    % (b/c the slope is 0 at i=1)
    if i==length(cl.x)
        m(1) = m(2); m_rev(1) = -1/m(1);
        % Calculate distance to points  
        delta_x = sqrt((segs.l)^2/(m_rev(1)^2+1));
        delta_y = sqrt((segs.l)^2/(m_rev(1)^2+1))*m_rev(1);        
        % Add dy/dx to the point on the centerline
        segs.x(1,2) = cl.x(1)+delta_x; segs.x(1,1) = cl.x(1)-delta_x;
        segs.y(1,2) = cl.y(1)+delta_y; segs.y(1,1) = cl.y(1)-delta_y; 
        % Increase number of points in segment
        segs.xn(1,:) = segs.x(1,1):(segs.x(1,2)-segs.x(1,1))*(dx/(segs.l*2)):segs.x(1,2);
        segs.yn(1,:) = segs.y(1,1):(segs.y(1,2)-segs.y(1,1))*(dx/(segs.l*2)):segs.y(1,2);            
        % Plot results
        plot(ax1, segs.xn(1,:),segs.yn(1,:),'-b','linewidth',1,'handlevisibility','off');
    end

end 

%% Cut off width segments where they intersect with the glacier outline polygon, calculate glacier width.

% initialize variables: 
% segment intersection points
segs.intx = NaN*zeros(length(cl.x),2);
segs.inty = NaN*zeros(length(cl.y),2);
% clipped width segments with ~100 m-spaced points (for width-averaging)
segs.xn_clip = NaN*zeros(length(cl.x),10e3/(dx*10));
segs.yn_clip = NaN*zeros(length(cl.x),10e3/(dx*10));
% loop through each width segment
for i=1:length(segs.x(:,1))
    % if segment intersects polygon...
    if any(inpolygon(segs.xn(i,:),segs.yn(i,:),ol.x,ol.y))
        % find the intersection points 
        [segs.intx(i,:),segs.inty(i,:), I] = polyxpoly(ol.x,ol.y,segs.xn(i,:),segs.yn(i,:));
        % calculate width from these segments
        segs.W(i) = sqrt((segs.intx(i,1)-segs.intx(i,2))^2+(segs.inty(i,1)-segs.inty(i,2))^2);
        % add points for segment that each have ~200 m spacing
        segs.xn_clip(i,1:length(segs.xn(i,min(I(:,2)):10:max(I(:,2))))) = segs.xn(i,min(I(:,2)):10:max(I(:,2)));
        segs.yn_clip(i,1:length(segs.yn(i,min(I(:,2)):10:max(I(:,2))))) = segs.yn(i,min(I(:,2)):10:max(I(:,2)));
    else 
        segs.intx(i,:) = NaN;
        segs.inty(i,:) = NaN;
        segs.W(i) = NaN;
    end
    
    % Plot resulting segments onto figure above
    figure(1); hold on;
    if i==1
        plot(ax1, segs.intx(i,:),segs.inty(i,:),'LineWidth',2,'color',col(3,:),...
            'DisplayName','width segments');
    else
        plot(ax1, segs.intx(i,:),segs.inty(i,:),'LineWidth',2,...
            'color',col(3,:),'HandleVisibility','off');        
    end
end 

% Define x as distance along the centerline
x = zeros(1,length(cl.x));
for i=2:length(cl.x)
    x(i) = sqrt((cl.x(i)-cl.x(i-1))^2+(cl.y(i)-cl.y(i-1)))+x(i-1);
end
    
% Load and plot the 2019 terminus position for reference
centerline = load('LarsenB_centerline.mat').centerline;
cl.termx = centerline.termx(end); cl.termy = centerline.termy(end); % last col = 2019
cl.Iterm = dsearchn([cl.x cl.y],[cl.termx cl.termy]); 
figure(1); plot(ax1, cl.x(cl.Iterm),cl.y(cl.Iterm),'x','linewidth',2,...
    'color',col(1,:),'markersize',10,'DisplayName','2019 terminus position');
% Plot resulting width along centerline
ax2 = subplot(1,3,3); hold on; 
    plot(ax2,segs.W/10^3,x/10^3,'linewidth',2,'color',col(2,:)); grid on;
    plot(ax2,segs.W(cl.Iterm)/10^3, x(cl.Iterm)/10^3, 'x', 'linewidth', 2,...
        'color',col(1,:),'markersize',10,'DisplayName','2019 terminus position');
    set(gca,'fontname','arial','fontsize',14,'linewidth',2);
    ylabel('distance along centerline [km]'); xlabel('glacier width [km]');

%% Save results! 

% Save figures
cd([homepath,'figures/']);
saveas(fig1,'glacier_width.png','png');
disp('figure saved');
% Save calculated width values
if save_width
    width.W = segs.W; width.x = cl.x; width.y = cl.y; 
        width.perpx = segs.x; width.perpy = segs.y;
        width.extx = segs.intx; width.exty = segs.inty;
        width.segsx = segs.xn; width.segsy = segs.yn;
        width.segsx_clip = segs.xn_clip; width.segsy_clip = segs.yn_clip;
    cd([homepath,'inputs-outputs/']);
    save('glacier_width.mat','width');
    disp(['saved width variables in: ',pwd]);
end


