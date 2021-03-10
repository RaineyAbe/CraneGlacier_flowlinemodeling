%% Script to tune the basal roughness factor beta, using 2009 observations of glacier
% speed and geometry.

%% 0. define time and space independent variables
  
clear all; close all;
warning off; % turn off warnings (velocity coefficient matrix is close to singular)
    
% define home path in directory and add necessary paths
homepath = '/Users/raineyaberle/Desktop/Research/CraneGlacier_modeling/';
addpath([homepath,'scripts/tuningBeta/']); % add path to betaSolve.m
addpath([homepath,'../matlabFunctions/']); % add path to bestPolyFit.m
cd([homepath,'inputs-outputs/']);

% grid spacing
dx0 = 1000:100:5000;

% Load Crane Glacier initialization variables
load('Crane_flowline_initialization.mat');
    
% Load 2009 observed ice surface speed
U_obs = load('Crane_CenterlineSpeeds_2007-2017.mat').U(6).speed;
        
% densities and g
rho_i = 917; % ice density (kg m^-3)
rho_sw = 1028; % ocean water density (kg m^-3)
rho_fw = 1000; % fresh water density (kg m^-3)
g = 9.81; % acceleration (m s^-2)

% stress parameters (unitless)
m = 1; % basal sliding exponent
n = 3; % flow law exponent

% calving parameters
Hc = 400; % m -> set the calving front to a default minimum ice thickness value if calving criteria is not met
fwd = 30; % fresh water depth in crevasses (m)
sigma_b = 0; % back pressure (Pa) due to sea ice or sikkusak

% maximum & minimum thickness & speed cut-off to check for instability
H_max = 2000; % maximum thickness (m)
H_min = 100;  % minimum thickness (m)
U_min = 200;  % minimum mean speed (m a^-1)

% initial conditions
H0 = h0-hb0; % ice thickness (m)
dUdx0 = [(U0(1:end-1)-U0(2:end))./(x0(1:end-1)-x0(2:end)) 0]; % strain rate (1/s)
% find the location of the grounding line and the end of the ice-covered domain
Hf = -(rho_sw./rho_i).*hb0; % flotation thickness (m)
gl0 = find(Hf-H0>0,1,'first')-1; % grounding line location 
H0(gl0:length(x0))=h0(gl0:end)*rho_sw/(rho_sw-rho_i); % buoyant thickness using surface
H0(H0>=(h0-hb0))=h0(H0>=(h0-hb0))-hb0(H0>=(h0-hb0)); % thickness can't go beneath bed elevation
H0(c0+1:end) = 0; % zero ice thickness past calving front

% extend variables past where there are NaNs (set to last value)
A0(find(isnan(A0),1,'first'):end) = A0(find(isnan(A0),1,'first')-1);
U0(find(isnan(U0),1,'first'):end) = U0(find(isnan(U0),1,'first')-1);
hb0(find(isnan(hb0),1,'first'):end) = hb0(find(isnan(hb0),1,'first')-1);
W0(find(isnan(W0),1,'first'):end) = W0(find(isnan(W0),1,'first')-1);    
E0 = ones(1,length(x0)); % enhancement factor

%% 1. run the flowline model

save_betabest = 1; % = 1 to save the beta with the lowest U RMSE

RMSE = zeros(1,length(dx0)); % pre-allocate RMSE
clear beta

% loop through grid spatial resolutions
for j=1:length(dx0)

    % save dx
    beta(j).dx = dx0(j);
    
    % initialize variables
    x=x0; h=h0; H=H0; hb=hb0; U=U0; A=A0; E=E0; dUdx=dUdx0;
    
    % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
    Rxx = 2*nthroot(dUdx./(E.*A),n); % resistive stress (Pa)
    crev = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*fwd); % crevasse penetration depth (m)
    % calving front located where the inland-most crevasse intersects sea level
    xcf = interp1(h-crev,x,0,'linear','extrap'); % (m along centerline)

    % calculate the thickness required to remain grounded at each grid cell
    Hf = -(rho_sw./rho_i).*hb; % flotation thickness (m)
    % find the location of the grounding line and use a floating
    % geometry from the grounding line to the calving front
    if ~isempty(find(Hf-H>0,1,'first'))
        xgl = interp1(Hf(find(Hf-H>0,1,'first')-1:find(Hf-H>0,1,'first')+1)...
            -H(find(Hf-H>0,1,'first')-1:find(Hf-H>0,1,'first')+1),...
            x(find(Hf-H>0,1,'first')-1:find(Hf-H>0,1,'first')+1),0,'linear','extrap'); % (m along centerline)
    else
        xgl=xcf;
    end
    if xgl>xcf % grounding line can't be past calving front
        xgl=xcf;
    end 

    % create coordinate system that hits cf and gl exactly
    % has resolution dxmax near the ice divide
    % has resolution dxmin from gl to c
    % and has smooth variation between
    xl = round(xgl/dx0(j)); %number of ideal grid spaces needed to reach the grounding line
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
    if isempty(c) == 1 %set the calving front to a default minimum ice thickness value
        c = find(H<Hc,1,'first'); 
    end
    if isempty(c)==1 % set c to length of x if criteria still not met
        c=length(x);
        disp('calving criteria not met');
    end
    hb = interp1(x0,hb0,xn,'pchip'); 
    W = interp1(x0,W0,xn,'pchip'); 
    H = interp1(x,H,xn,'pchip'); 
    U = interp1(x,U,xn,'pchip'); 
    A = interp1(x0,A0,xn,'pchip'); 
    E = interp1(x,E,xn,'pchip'); 
    x = xn; dx = dxn;

    % calculate surface elevation
    h = hb+H; % surface elevation (m a.s.l.)
    h(gl+1:c) = (1-rho_i/rho_sw).*H(gl+1:c); %adjust the surface elevation of ungrounded ice to account for buoyancy
    H(h<0)=0-hb(h<0); h(h<0)=0; % surface cannot go below sea level

    % calculate the effective pressure (ice overburden pressure minus water
    % pressure) assuming an easy & open connection between the ocean and
    % ice-bed interface
    sl = find(hb<=0,1,'first'); % find where the glacier base first drops below sea level
    N_ground = rho_i*g*H(1:sl); % effective pressure where the bed is above sea level (Pa)
    N_marine = rho_i*g*H(sl+1:length(x))+(rho_sw*g*hb(sl+1:length(x))); % effective pressure where the bed is below sea level (Pa)
    N = [N_ground N_marine];
    N(N<0)=0; % cannot have negative values

    % solve for new beta, save results
    [beta(j).beta,dUdx,Un_rev] = betaSolve(H,c,x,U,n,A,E,m,dx,rho_i,g,h,rho_sw,sigma_b,W,N);    
    beta(j).x = x; beta(j).U = Un_rev; 
    
    % plot geometry, speed, & calving front position
    if j==1
        col = parula(length(dx0)+2); % color scheme for plots
        figure(1); clf
        set(gcf,'Position',[0 100 1300 400]);
        ax1 = axes('Position',[0.05 0.1 0.25 0.8]); % glacier geometry
            hold on; grid on; set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
            title('Glacier Geometry'); xlim([0 50]); ylim([min(hb)-100 max(h)+200]);
            xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)');
            % ice surface
            plot(x(1:c)./10^3,h(1:c),'color',col(j,:),'linewidth',2,'displayname',num2str(dx0(j)));
            % calving front
            plot(x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(j,:),'linewidth',2,'HandleVisibility','off');
            % floating bed
            plot(x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(j,:),'linewidth',2,'HandleVisibility','off');
            % bed elevation
            plot(x0./10^3,hb0,'k','linewidth',2,'HandleVisibility','off');
            % mean sea level
            plot([x(1),x(end)]/10^3,[0,0],'k--','HandleVisibility','off');
            % colorbar
            colormap(parula(length(col))); cb=colorbar('Ticks',[0 0.5 1],...
                'YTickLabel',[{num2str(dx0(1)/10^3)} {num2str(dx0(round(length(dx0)/2))/10^3)} {num2str(dx0(end)/10^3)}],...
                'Position',[.9 .35 .025 .3410],'Fontname','arial','fontweight','bold');  
                cb.Label.String = 'dx (km)';
        ax2 = axes('Position',[0.35 0.1 0.25 0.8]); % ice speed
            hold on; grid on; set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
            title('Ice Speed'); xlim([0 50]); ylim([0 2500]);
            xlabel('Distance Along Centerline (km)'); ylabel('Speed (m yr^{-1})');
            plot(x(1:c)./10^3,Un_rev(1:c).*3.1536e7,'color',col(j,:),'linewidth',2);
        ax3 = axes('Position',[0.7 0.1 0.25 0.8]); % beta
            hold on; grid on; set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
            title('Basal Roughness Factor \beta'); xlim([0 50]);
            xlabel('Distance Along Centerline (km)'); ylabel('\beta (s^{1/m} m^{-1/m})');
            plot(x(1:c-1)./10^3,beta(j).beta(1:c-1),'color',col(j,:),'linewidth',2);
    else
        % ice surface
        plot(ax1,x(1:c)./10^3,h(1:c),'color',col(j,:),'linewidth',2);
        % calving front
        plot(ax1,x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(j,:),'linewidth',2);
        % floating bed
        plot(ax1,x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(j,:),'linewidth',2);        
        % ice speed
        plot(ax2,x(1:c)./10^3,Un_rev(1:c).*3.1536e7,'color',col(j,:),'linewidth',2);
        % beta
        plot(ax3,x(1:c-1)./10^3,beta(j).beta(1:c-1),'color',col(j,:),'linewidth',2);
    end
    if j==length(dx0)
        figure(1);
        % colorbar
        colormap(parula(length(col))); cb=colorbar('Ticks',[0 0.5 1],...
            'YTickLabel',[{num2str(dx0(1)/10^3)} {num2str(dx0(round(length(dx0)/2))/10^3)} {num2str(dx0(end)/10^3)}],...
            'Position',[.95 .35 .025 .3410],'Fontname','arial','fontweight','bold');  
            cb.Label.String = 'dx (km)';        
    end

    % Calculate RMSE of resulting speed using beta solution
    RMSE(j) = sqrt((sum(Un_rev-U).^2./length(U)));

    % display results
    disp(['dx = ',num2str(dx0(j)),' m: ']);
    disp(['RMSE = ',num2str(RMSE(j)),' m/s']);

end 

% plot RMSE results
Ibest = find(RMSE==min(RMSE)); % index of lowest RMSE
if length(Ibest)>1 % if lowest RMSE exists for multiple dx, use lowest dx
    Ibest(2:end)=[];
end
figure(4); clf
set(gcf,'Position',[441   145   894   652]);
subplot(2,2,1:2);
    set(gca,'fontsize',14,'fontname','arial','linewidth',2);
    xlabel('x (km along centerline)'); ylabel('\beta (s^{1/m} m^{-1/m})'); 
    title('\beta_{best}'); hold on; grid on; 
    plot(beta(Ibest).x./10^3,beta(Ibest).beta,'linewidth',2);  
subplot(2,2,3);
    set(gca,'fontsize',14,'fontname','arial','linewidth',2);
    xlabel('dx (km)'); ylabel('RMSE (m/s)'); title('RMSE of \beta Solutions');
    hold on; plot(dx0./10^3,RMSE,'.','markersize',20); grid on;
    plot(dx0(Ibest)/10^3,RMSE(Ibest),'*','markersize',20,'linewidth',2);  
subplot(2,2,4);
    set(gca,'fontsize',14,'fontname','arial','linewidth',2);
    xlabel('x (m along centerline)'); ylabel('U (m a^{-1})'); title('U Solutions');
    hold on; grid on; legend('Location','best');
    plot(x0(1:c0)./10^3,U0(1:c0).*3.1536e7,'b','displayname','observed','linewidth',2);
    plot(beta(Ibest).x./10^3,beta(Ibest).U.*3.1536e7,'displayname','solution','linewidth',2);    

% save beta with lowest U RMSE
if save_betabest
    cd([homepath,'inputs-outputs/']);
    betabest.beta = beta(Ibest).beta; 
    betabest.x = beta(Ibest).x;
    betabest.dx = beta(Ibest).dx;
    save('betabest.mat','betabest');
    disp('betabest saved.');
end


