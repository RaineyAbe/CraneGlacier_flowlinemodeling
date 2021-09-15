%% Script to augment the basal roughness factor (beta) solved by gradient 
% descent (in tuneParameters.m) to further reduce model speed misfits
% using observations for 2009-2018
% Adapted from Ellyn Enderlin's flowline model (Enderlin et al., 2013)
% Rainey Aberle
% 2020-2021

%% 0. define time and space independent variables

clear all; close all;
warning off; % turn off warnings (velocity coefficient matrix is close to singular)
    
% define home path in directory and add necessary paths
homepath = '/Users/raineyaberle/Desktop/Research/CraneModeling/CraneGlacier_flowlinemodeling/';
addpath([homepath,'scripts/3_sensitivityTests/']); % add path to U_convergence
addpath([homepath,'../matlabFunctions/cmocean_v2.0/cmocean']);
addpath([homepath,'inputs-outputs/']);

dx0 = 200; % grid spacing (m)

% Load Crane Glacier initialization variables
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
% submarine melt rate profile fit
SMR_mean_fit = load('LarsenC_MeanMeltRate.mat').mr_mean_fit;
        
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
H0(c0:end) = 0; % no ice thickness past the calving front
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

%% 1. Define basal roughness factor beta

% Note: Run Section 0 before re-running this section each time

% add Gaussian bump where speed is too slow
width = -3e3:dx0:3e3; % x distance
sigma = 8e3; % bumpwidth, kind of
mu = 1; 
gauss = normpdf(width,mu,sigma);
gauss = gauss*25e3; % increase magnitude
% add bump between 5 km and 25km 
beta0 = 0.9*ones(1,length(x0));
beta0(dsearchn(x0',24e3):dsearchn(x0',30e3)) = beta0(dsearchn(x0',24e3):dsearchn(x0',30e3))+gauss;
%        0   5   10  15   20   25  30  35 40  45  50  55  60
% beta = [0.9 1.1 0.6 0.2 0.2 1.0 1.5 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0];
% betax = 0:(x0(end)/(length(beta)-1)):x0(end);
% beta0 = interp1(betax,beta,x0);
% beta0 = movmean(beta0,10);
figure(5); clf; plot(x0,beta0);

%% 2. Run the flowline mode and calculate model misfits
    
% Note: Must rerun previous sections before running this section each time

plot_timeSteps = 1; % = 1 to plot geometry/speed at each model year
save_figure = 0; % = 1 to save final figure

% define time stepping (s)
dt = 0.01*3.1536e7;
t_start = 0*3.1536e7;
t_end = 9*3.1536e7;
t = (t_start:dt:t_end);

% initialize ice flux
Fgl=zeros(1,length(t)); % ice mass flux across the grounding line

% initialize variables
x=x0; U=U0; W=W0; gl=gl0; dUdx=dUdx0; A=A0; h=h0; hb=hb0; H=H0; d_fw=d_fw0; beta = beta0;

% run flowline model
for i=1:length(t)

    % implement change to d_fw after 10 model years
    if t(i)>10*3.1536e7
        % linearly increase d_fw to reach delta_d_fw by 2100
        delta_d_fwi = delta_d_fw/(2100-2019)*(t(i)/3.1536e7-10); % total increase in d_fw at this time increment
        d_fw = d_fw0 + delta_d_fwi; % m
    end 

    % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
    Rxx = 2*nthroot(dUdx./(E.*A),n); % resistive stress (Pa)
    crev_s = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*d_fw); % surface crevasse penetration depth (m)
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
        if xcf<20e3 || xcf > 100e3 || isnan(xcf)
            xcf = x(dsearchn(x',x(c)));
            %xcf = interp1(feval(fit(x',(h-crev_s)','poly1'),x),x,0,'linear','extrap');
            %xcf = x0(c0);
        end
    end
    XCF(i) = xcf; % save calving front position over time

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
    if isempty(c) == 1 %set the calving front to a default minimum ice thickness value
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
    if plot_timeSteps
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
    [U,dUdx,Td,Tlatb,Tlon,vm] = U_convergence(x,U,U0,dUdx,H,h,A,E,N,W,dx,c,n,m,beta,rho_i,rho_sw,g,sigma_b,i);
    if t(i)==t_end
        figure(9); hold on;
        plot(x,Tlatb,'-b','displayname','T_{latb}'); 
        plot(x,Tlon,'-k','displayname','T_{lon}');
    end 

    % calculate ice flux
    F = U.*H.*W; % ice flux (m^3 s^-1)
    F(isnan(F))=0;
    F(1)=F(2)+F0;

    % implement SMB, SMR, & delta_SMB
    SMR = zeros(1,c);
    % use the Larsen C mean melt rate profile to scale smr
    % using the max initial smr
    SMR(gl+1:c) = SMR0/(SMR_mean_fit.a+1)*feval(SMR_mean_fit,x(gl+1:c)-x(gl));
    % Define SMB and estimated runoff
    SMB = interp1(x0,SMB0,x);
    RO = interp1(x0,RO0,x);                      

    % calculate the  change in ice thickness from continuity
    clearvars dHdt
    dHdt(1) = (-1/W(1))*(F(1)-F(2))/(x(1)-x(2)); % forward difference
    dHdt(2:c-1) = (-1./W(2:c-1)).*(F(1:c-2)-F(3:c))./(x(1:c-2)-x(3:c)); % central difference
    dHdt(c:length(x)) = (-1./W(c:length(x))).*(F(c-1:length(x)-1)-F(c:length(x)))./(x(c-1:length(x)-1)-x(c:length(x))); % backward difference
    dH = dHdt.*dt;

    % new thickness (change from dynamics, SMB, & SMR)
    Hn = H+dH+(SMB.*dt)+(SMR.*dt)-(interp1(x0,RO0,x)*dt)+(interp1(x0,Q0,x)*dt);
    Hn(Hn < 0) = 0; % remove negative values
    H = Hn; % set as the new thickness value

    Fgl(i) = F(gl)*pi*1/4*917*1e-12*3.1536e7; % Gt/a

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
b=1;
while b==1
    figure(2); clf;
    set(gcf,'position',[200 200 1000 700],'defaultAxesColorOrder',[[0 0 0];[0.8 0.1 0.1]]);
    ax1 = axes('position',[0.08 0.67 0.36 0.3]); hold on; grid on;
        set(gca,'fontsize',18,'fontname','Arial','linewidth',2); 
        xlim([0 52]); ylabel('Misfit (m)'); 
        plot(x/10^3,h-interp1(cl.x,h_obs(36).surface,x),'-','linewidth',2,...
            'color',[0.8 0.1 0.1],'HandleVisibility','off');
        plot(x/10^3,mean(h-interp1(cl.x,h_obs(36).surface,x),'omitnan')*ones(1,length(x)),...
            '--','linewidth',2,'color',[0.8 0.1 0.1],'HandleVisibility','off');
        text((max(get(gca,'XLim'))-min(get(gca,'XLim')))*0.92+min(get(gca,'XLim')),...
                (max(get(gca,'YLim'))-min(get(gca,'YLim')))*0.93+min(get(gca,'YLim')),...
                '(a)','backgroundcolor','w','fontsize',18,'linewidth',1);  
    ax2 = axes('position',[0.56 0.67 0.36 0.3]); hold on; grid on;
        set(gca,'fontsize',18,'fontname','Arial','linewidth',2); 
        xlim([0 52]); ylabel('Misfit (m a^{-1})');
        plot(x/10^3,(U-interp1(cl.x,U_obs(10).U,x))*3.1536e7,'-','linewidth',2,...
            'color',[0.8 0.1 0.1],'HandleVisibility','off');
        plot(x/10^3,mean(U-interp1(cl.x,U_obs(10).U,x),'omitnan')*3.1536e7*ones(1,length(x)),...
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
        yyaxis right; ylabel('Basal Roughness Factor (s^{1/m} m^{-1/m})'); ylim([0 max(beta)+0.2]);
            plot(x/10^3,beta,'-','linewidth',2,'color',[0 0.4 0.8],'displayname','\beta');

    % display grounding line misfit
    disp('Mean misfits at the grounding line:');
    disp(['    Speed: ',num2str((U(gl)-interp1(cl.x,U_obs(10).U,x(gl)))*3.1536e7),' m/yr']);
    disp(['    Elevation: ',num2str((h(gl)-interp1(cl.x,h_obs(36).surface,x(gl)))),' m']);
    b=b+1; % exit while loop 
end 

% save figure
if save_figure
    cd([homepath,'figures/']);
    saveas(gcf,'misfits2018.png','png');
    cd([homepath,'../write-ups/Thesis/figures/']);
    saveas(gcf,'misfits2018.png','png');
    disp('figure saved');
end
