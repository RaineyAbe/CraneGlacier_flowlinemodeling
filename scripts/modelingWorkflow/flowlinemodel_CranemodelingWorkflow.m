%% Glacier Flowline Model: Modeling workflow for Crane Glacier, Antarctic Peninsula
% Last edited: Spring 2021
% Adapted by Rainey Aberle from Ellyn Enderlin's flowline model package 
% (Enderlin et al., 2013)
%
% Workflow: Run step 0 before any other step to initialize variables.
%   0. Define time and space independent variables by loading the flowline 
%       initialization file and regridding to the desired grid spacing. 
%   1. Tune the fresh water depth in crevasses fwd to minimize the misfit
%       between modeled and observed terminus positions 2009-2019.
%   2. Tune the enhancement factor E along the centerline to minimize
%       misfit between modeled and observed ice surface elevation.
%   3. Tune the backstress sigma_b using the resistive stress term to 
%       minimize misfit between modeled and observed terminus position 
%       to account for sea ice or sikkusak (Nick et al., 2010)
%   4. Run sensitivity tests for surface mass balance (SMB) & submarine
%       melting rate (SMR). 

clear all; close all;
warning off; % turn off warnings (velocity coefficient matrix is close to singular)
    
% define home path in directory and add necessary paths
homepath = '/Users/raineyaberle/Desktop/Research/CraneGlacier_flowlinemodeling/';
addpath([homepath,'scripts/modelingWorkflow/']); % add path to U_convergence
cd([homepath,'inputs-outputs/']);

%% 0. define time and space independent variables

dx0 = 200; % grid spacing (m)

% Load Crane Glacier initialization variables
load('Crane_flowlineModelInitialization.mat');
    
% Load observed conditions
% dH
dH_obs = load('Crane_dHdt_2009-2018.mat').dHdt.dH_total; % (m) total change in thickness 2009-2018
% terminus position 
term = load('Crane_terminusPositions_2002-2019.mat').term;
for i=1:length(term)
    termx_obs(i) = term(i).x;
    termDate_obs(i) = term(i).decidate;
end
% fit a quadratic to the terminus positions to smooth seasonal variations
termx_obs = feval(fit(termDate_obs',termx_obs','poly2'),termDate_obs');
term_obs = interp1(termDate_obs',termx_obs,2009:2017);
clear term 
% ice speed
U_obsi = load('Crane_centerlineSpeeds_2007-2017.mat').U;
u = [6 8 9 14 15:19]; % indices of speeds to use annually (2009-2017)
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

% calving parameters
Hc = 400; % m -> set the calving front to a default minimum ice thickness value if calving criteria is not met
fwd = 5; % fresh water depth in crevasses (m) - similar to RACMO modeled and downscaled snowfall at the terminus for 2016
sigma_b = 0; % back pressure (Pa) - similar to that employed at Helheim Glacier (Nick et al., 2009)

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

% extend variables to model length if necessary (set extended points to last value)

A0(find(isnan(A0),1,'first'):end) = A0(find(isnan(A0),1,'first')-1);
U0(find(isnan(U0),1,'first'):end) = U0(find(isnan(U0),1,'first')-1);
beta0(find(isnan(beta0),1,'first'):end) = beta0(find(isnan(beta0),1,'first')-1);
hb0(find(isnan(hb0),1,'first'):end) = hb0(find(isnan(hb0),1,'first')-1);
W0(find(isnan(W0),1,'first'):end) = W0(find(isnan(W0),1,'first')-1);    
E0 = ones(1,length(x0)); % enhancement factor

%% 1. Tune enhancement factor E (via U & h)

close all; 

save_optE = 1; % = 1 to save Ebest
   
Efit = (25:30)'.*ones(1,length(x0));

% pre-allocate misfit variables
dH_tot = NaN.*zeros(length(Efit(:,1)),length(x0)); % total dH for each point along centerline
dH_RMSE = NaN.*zeros(length(Efit(:,1)),1); % misfit of total dH
dU_RMSE = NaN.*zeros(length(Efit(:,1)),1);

% loop through all Efit values
for f=1:length(Efit(:,1))
    
    % define time stepping (s)
    dt = 0.01*3.1536e7;
    t_start = 0*3.1536e7;
    t_end = 10*3.1536e7;
    t = (t_start:dt:t_end);
    
    % re-initialize variables
    x=x0; H=H0; U=U0; dUdx=dUdx0; A=A0; h=h0; hb=hb0; gl=gl0; c=c0;
    
    % define E
    E = Efit(f,:);
    sigma_b=0;
    
    % run flowline model
    for i=1:length(t)
        
        % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
        Rxx = 2*nthroot(dUdx./(E.*A),n); % resistive stress (Pa)
        crev = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*fwd); % crevasse penetration depth (m)
        % calving front located where the inland-most crevasse intersects sea level
        xcf = interp1(h-crev,x,0,'linear','extrap'); % (m along centerline)
        if i==1 % use observed calving front position for first iteration
            xcf = x0(c0);
        end
        
        % calculate the thickness required to remain grounded at each grid cell
        Hf = -(rho_sw./rho_i).*hb; % flotation thickness (m)
        % find the location of the grounding line and use a floating
        % geometry from the grounding line to the calving front
        if ~isempty(find(Hf-H>0,1,'first'))
            %xgl = x(find(Hf-H>0,1,'first')-1);
            xgl = interp1(Hf(find(Hf-H>0,1,'first')-1:find(Hf-H>0,1,'first')+1)...
               -H(find(Hf-H>0,1,'first')-1:find(Hf-H>0,1,'first')+1),...
               x(find(Hf-H>0,1,'first')-1:find(Hf-H>0,1,'first')+1),0,'linear','extrap'); % (m along centerline)
        else
            xgl=xcf;
        end
        if xgl>xcf % grounding line can't be past calving front
            xgl=xcf;
        elseif xgl<0 % grounding line can't be negative
            xgl=xcf;
        end
        
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
        E = interp1(x,E,xn,'linear','extrap');
        beta = interp1(x0,beta0,xn,'linear','extrap');
        x = xn; dx = dxn;
        
        % calculate surface elevation
        h = hb+H; % surface elevation (m a.s.l.)
        h(gl+1:c) = (1-rho_i/rho_sw).*H(gl+1:c); %adjust the surface elevation of ungrounded ice to account for buoyancy
        H(h<0)=0-hb(h<0); h(h<0)=0; % surface cannot go below sea level
        h(h-H<hb) = hb(h-H<hb)+H(h-H<hb); % thickness cannot go beneath bed elevation
        
        % plot geometry, speed, & calving front position
        if t(i)==t_start
            col = parula(length(t)+20); % color scheme for plots
            figure(1); clf
            set(gcf,'Position',[0 100 1300 400]);
            ax1 = axes('Position',[0.05 0.1 0.28 0.8]); % glacier geometry
            hold on; grid on;
            set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
            title('Glacier Geometry'); legend('Location','northeast');
            xlim([0 65]); ylim([min(hb)-100 max(h)+200]);
            xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)');
            % ice surface
            plot(x(1:c)./10^3,h(1:c),'color',col(i,:),'linewidth',2,'displayname','0');
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
            set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
            title('Ice Speed'); legend('Location','northeast');
            xlim([0 65]); ylim([0 2500]);
            xlabel('Distance Along Centerline (km)'); ylabel('Speed (m yr^{-1})');
            % ice speed
            plot(x(1:c)./10^3,U(1:c).*3.1536e7,'color',col(i,:),'linewidth',2,'displayname','0');
            ax3 = axes('Position',[0.7 0.1 0.28 0.8]); % calving front position
            hold on; grid on;
            set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
            title('Calving Front Position'); legend('Location','best');
            xlim([30 65]); ylim([0 10]);
            xlabel('Distance Along Centerline (km)'); ylabel('Year');
            % calving front position
            plot(x(c)/10^3,t(i)./3.1536e7,'.','markersize',15,'color',col(i,:),'displayname','0');
        elseif mod(i-1,round(length(t)/50))==0 % display every length(t)/50
            figure(1);
            if mod(i-1,round(length(t)/10))==0 % display in legend every length(t)/10
                % ice surface
                plot(ax1,x(1:c)/10^3,h(1:c),'-','color',col(i,:),'linewidth',2,'displayname',num2str(t(i)./3.1536e7));
                % ice speed
                plot(ax2,x(1:c)/10^3,U(1:c).*3.1536e7,'-','Color',col(i,:),'linewidth',2,'DisplayName',num2str(t(i)./3.1536e7));
                % calving front position
                plot(ax3,x(c)./10^3,t(i)./3.1536e7,'.','Color',col(i,:),'markersize',15,'displayname',num2str(t(i)./3.1536e7)); hold on;
            else
                % ice surface
                plot(ax1,x(1:c)/10^3,h(1:c),'-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
                % ice speed
                plot(ax2,x(1:c)/10^3,U(1:c).*3.1536e7,'-','Color',col(i,:),'linewidth',2,'HandleVisibility','off'); hold on;
                % calving front position
                plot(ax3,x(c)./10^3,t(i)./3.1536e7,'.','Color',col(i,:),'markersize',15,'HandleVisibility','off'); hold on;
            end
            % calving front
            plot(ax1,x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
            % floating bed (gl:c)
            plot(ax1,x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off');
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
        [U,dUdx,vm,Un,dUndx,M,T] = U_convergence(x,U,dUdx,H,h,A,E,N,W,dx,c,n,m,beta,rho_i,rho_sw,g,sigma_b,i);
        
        % calculate ice flux
        F = U.*H.*W; % ice flux (m^3 s^-1)
        F(isnan(F))=0;
        
        % implement SMB & SMR
        smr = zeros(1,length(x));
        smr(gl+1:length(x)) = feval(fit([x(gl);x(gl+1);x(c)],[0;smr0;0.95*smr0],'pchip'),x(gl+1:end));
        smb = interp1(x0,smb0+Q0,x);
        
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
        
        % calculate RMSE of surface speed and elevation change at the
        % final model time
        if t(i)==t_end
            U_RMSE(f) = sqrt(nansum(interp1(x,U,x0)-U_obs(9).U').^2/length(x0));
            dH_tot(f,:) = interp1(x,h,x0)-h0;
            dH_RMSE(f) = sqrt((nansum(dH_tot(f,:)-dH_obs).^2)/length(x0));
        end
        
    end
    
end

% calculate mean misfit
misfit = nanmean(dH_RMSE,2).*nanmean(U_RMSE,2);

% plot results for optimal E
IoptE = find(abs(misfit)==min(abs(misfit)),1,'first');
optE = Efit(IoptE,:);
figure(10); clf
    hold on; grid on; 
    set(gca,'fontsize',14,'linewidth',2); 
    xlabel('Efit'); ylabel('misfit (m*m/s)');
    %title(['E_{best} = ',num2str(Efit(IEbest).fit(2)),'+',num2str(Efit(IEbest).fit(1)),'x']);
    title(['E_{best} = ',num2str(optE(1))]);
    plot(Efit(:,1),misfit,'.','markersize',20);
    plot(Efit(IoptE,1),misfit(IoptE),'*','markersize',15,'linewidth',2);

% save Ebest
if save_optE
    cd([homepath,'inputs-outputs']);
    save('optimalE.mat','optE','x0');
    disp('optimal E saved');
end

%% 2. Tune fresh water depth in crevasses fwd (via c)

close all;

save_optfwd = 1; % = 1 to save fwdbest

% define fwd values to test
fwd0 = 16:1:20; % fresh water depth in crevasses (m)
    
% pre-allocate misfit variables
c_misfit = NaN.*zeros(length(fwd0),length(2009:2018)); % calving front misfit for each model year        
        
% Loop through all fwd values
for f=1:length(fwd0)
    
    % define time stepping (s)
    dt = 0.01*3.1536e7;
    t_start = 0*3.1536e7;
    t_end = 10*3.1536e7;
    t = (t_start:dt:t_end);
    
    % initialize variables
    x=x0; H=H0; U=U0; dUdx=dUdx0; A=A0; h=h0; hb=hb0;
    
    % load optimal parameters 
    E=load('optimalE.mat').optE; % unitless
    
    % define fwd
    fwd = fwd0(f); 
    
    % run flowline model
    for i=1:length(t)
        
        % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
        Rxx = 2*nthroot(dUdx./(E.*A),n); % resistive stress (Pa)
        crev = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*fwd); % crevasse penetration depth (m)
        % calving front located where the inland-most crevasse intersects sea level
        %xcf = interp1(feval(fit(x',(h-crev)','poly1'),x),x,0,'linear','extrap'); % (m along centerline)
        xcf = interp1(h-crev,x,0,'linear','extrap'); % (m along centerline)
        if i==1 % use observed calving front for first iteration
            xcf = x0(c0);
        end
        
        % calculate the thickness required to remain grounded at each grid cell
        Hf = -(rho_sw./rho_i).*hb; % flotation thickness (m)
        % find the location of the grounding line and use a floating
        % geometry from the grounding line to the calving front
        if ~isempty(find(Hf-H>0,1,'first'))
            xgl = x(find(Hf-H>0,1,'first')-1);
            %xgl = interp1(Hf(find(Hf-H>0,1,'first')-1:find(Hf-H>0,1,'first')+1)...
            %    -H(find(Hf-H>0,1,'first')-1:find(Hf-H>0,1,'first')+1),...
            %    x(find(Hf-H>0,1,'first')-1:find(Hf-H>0,1,'first')+1),0,'linear','extrap'); % (m along centerline)
        else
            xgl=xcf;
        end
        if xgl>xcf % grounding line can't be past calving front
            xgl=xcf;
        elseif xgl<0 % grounding line can't be negative
            xgl=xcf;
        end
        
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
        E = interp1(x,E,xn,'linear','extrap');
        beta = interp1(x0,beta0,xn,'linear','extrap');
        x = xn; dx = dxn;
        
        % calculate surface elevation
        h = hb+H; % surface elevation (m a.s.l.)
        h(gl+1:c) = (1-rho_i/rho_sw).*H(gl+1:c); %adjust the surface elevation of ungrounded ice to account for buoyancy
        H(h<0)=0-hb(h<0); h(h<0)=0; % surface cannot go below sea level
        h(h-H<hb) = hb(h-H<hb)+H(h-H<hb); % thickness cannot go beneath bed elevation        
        
        % plot geometry, speed, & calving front position
        if t(i)==t_start
            col = parula(length(t)+20); % color scheme for plots
            figure(1); clf
            set(gcf,'Position',[0 100 1300 400]);
            ax1 = axes('Position',[0.05 0.1 0.28 0.8]); % glacier geometry
            hold on; grid on;
            set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
            title('Glacier Geometry'); legend('Location','northeast');
            xlim([0 65]); ylim([min(hb)-100 max(h)+200]);
            xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)');
            % ice surface
            plot(x(1:c)./10^3,h(1:c),'color',col(i,:),'linewidth',2,'displayname','0');
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
            set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
            title('Ice Speed'); legend('Location','northeast');
            xlim([0 65]); ylim([0 2500]);
            xlabel('Distance Along Centerline (km)'); ylabel('Speed (m yr^{-1})');
            % ice speed
            plot(x(1:c)./10^3,U(1:c).*3.1536e7,'color',col(i,:),'linewidth',2,'displayname','0');
            ax3 = axes('Position',[0.7 0.1 0.28 0.8]); % calving front position
            hold on; grid on;
            set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
            title('Calving Front Position'); legend('Location','best');
            xlim([30 65]); ylim([0 10]);
            xlabel('Distance Along Centerline (km)'); ylabel('Year');
            % calving front position
            plot(x(c)/10^3,t(i)./3.1536e7,'.','markersize',15,'color',col(i,:),'displayname','0');
        elseif mod(i-1,round(length(t)/50))==0 % display every length(t)/50
            figure(1);
            if mod(i-1,round(length(t)/10))==0 % display in legend every length(t)/10
                % ice surface
                plot(ax1,x(1:c)/10^3,h(1:c),'-','color',col(i,:),'linewidth',2,'displayname',num2str(t(i)./3.1536e7));
                % ice speed
                plot(ax2,x(1:c)/10^3,U(1:c).*3.1536e7,'-','Color',col(i,:),'linewidth',2,'DisplayName',num2str(t(i)./3.1536e7));
                % calving front position
                plot(ax3,x(c)./10^3,t(i)./3.1536e7,'.','Color',col(i,:),'markersize',15,'displayname',num2str(t(i)./3.1536e7)); hold on;
            else
                % ice surface
                plot(ax1,x(1:c)/10^3,h(1:c),'-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
                % ice speed
                plot(ax2,x(1:c)/10^3,U(1:c).*3.1536e7,'-','Color',col(i,:),'linewidth',2,'HandleVisibility','off'); hold on;
                % calving front position
                plot(ax3,x(c)./10^3,t(i)./3.1536e7,'.','Color',col(i,:),'markersize',15,'HandleVisibility','off'); hold on;
            end
            % calving front
            plot(ax1,x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
            % floating bed (gl:c)
            plot(ax1,x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off');
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
        [U,dUdx,vm,Un,dUndx,M,T] = U_convergence(x,U,dUdx,H,h,A,E,N,W,dx,c,n,m,beta,rho_i,rho_sw,g,sigma_b,i);
        
        % calculate ice flux
        F = U.*H.*W; % ice flux (m^3 s^-1)
        F(isnan(F))=0;
        
        % implement SMB & SMR
        smr = zeros(1,length(x));
        smr(gl+1:length(x)) = feval(fit([x(gl);x(gl+1);x(c)],[0;smr0;0.95*smr0],'pchip'),x(gl+1:end));
        smb = interp1(x0,smb0+Q0,x);
        
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
   
        % calculate misfit in calving front position at each full model year
        if t(i)==t_end %mod(t(i),3.1536e7)==0 && t(i)~=t_start
            %c_misfit(f,t(i)./3.1536e7) = x(c)-interp1(termDate_obs,termx_obs,t(i)./3.1536e7+2009);
            c_misfit(f) = x(c)-interp1(termDate_obs,termx_obs,t(i)./3.1536e7+2009);
        end
        
    end 

end

% calculate mean misfit
misfit = nanmean(c_misfit,2); %mean(c_misfit,2);

% plot results
Ioptfwd = find(abs(misfit)==min(abs(misfit)),1,'first');
optfwd = fwd0(Ioptfwd); % m.w.e.
figure(10); clf
    set(gca,'fontsize',14,'fontname','Arial','linewidth',2);
    grid on; xlabel('fwd (m)'); ylabel('misfit'); hold on;
    plot(fwd0,misfit,'.','markersize',20);
    plot(optfwd,misfit(Ioptfwd),'*','markersize',15,'linewidth',2);

% save results
if save_optfwd
    cd([homepath,'inputs-outputs']);
    save('optimalfwd.mat','optfwd');
    disp(['Optimal fwd = ',num2str(optfwd),' m saved.']);
end

%% 3. Tune for the back pressure sigma_b (via c) 

close all; 

save_optSigma_b = 1; % = 1 to save sigma_bbest
    
% define sigma_b values to test
sigma_b0 = 0e3:1e3:5e3; % back pressure (Pa)
    
% pre-allocate misfit variables
c_misfit = NaN.*zeros(length(sigma_b0),length(2009:2018)); % calving front misfit for each model year 
    
% loop through all sigma_b values
for f=1:length(sigma_b0)
    
    % define time stepping (s)
    dt = 0.01*3.1536e7;
    t_start = 0*3.1536e7;
    t_end = 10*3.1536e7;
    t = (t_start:dt:t_end);
    
    % initialize variables
    x=x0; H=H0; U=U0; dUdx=dUdx0; A=A0; h=h0; hb=hb0;
    
    % load optimal parameters 
    E=load('optimalE.mat').optE; % unitless
    fwd=load('optimalfwd.mat').optfwd; % m 
    
    % define sigma_b
    sigma_b = sigma_b0(f);
    
    % run flowline model
    for i=1:length(t)
        
        % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
        Rxx = 2*nthroot(dUdx./(E.*A),n); % resistive stress (Pa)
        crev = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*fwd); % crevasse penetration depth (m)
        % calving front located where the inland-most crevasse intersects sea level
        xcf = interp1(h-crev,x,0,'linear','extrap'); % (m along centerline)
        if i==1 % use observed calving front for first iteration
            xcf = x0(c0);
        end
        
        % calculate the thickness required to remain grounded at each grid cell
        Hf = -(rho_sw./rho_i).*hb; % flotation thickness (m)
        % find the location of the grounding line and use a floating
        % geometry from the grounding line to the calving front
        if ~isempty(find(Hf-H>0,1,'first'))
            xgl = x(find(Hf-H>0,1,'first')-1);
            %xgl = interp1(Hf(find(Hf-H>0,1,'first')-1:find(Hf-H>0,1,'first')+1)...
            %    -H(find(Hf-H>0,1,'first')-1:find(Hf-H>0,1,'first')+1),...
            %    x(find(Hf-H>0,1,'first')-1:find(Hf-H>0,1,'first')+1),0,'linear','extrap'); % (m along centerline)
        else
            xgl=xcf;
        end
        if xgl>xcf % grounding line can't be past calving front
            xgl=xcf;
        elseif xgl<0 % grounding line can't be negative
            xgl=xcf;
        end
        
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
        E = interp1(x,E,xn,'linear','extrap');
        beta = interp1(x0,beta0,xn,'linear','extrap');
        x = xn; dx = dxn;
        
        % calculate surface elevation
        h = hb+H; % surface elevation (m a.s.l.)
        h(gl+1:c) = (1-rho_i/rho_sw).*H(gl+1:c); %adjust the surface elevation of ungrounded ice to account for buoyancy
        H(h<0)=0-hb(h<0); h(h<0)=0; % surface cannot go below sea level
        h(h-H<hb) = hb(h-H<hb)+H(h-H<hb); % thickness cannot go beneath bed elevation
        
        % plot geometry, speed, & calving front position
        if t(i)==t_start
            col = parula(length(t)+20); % color scheme for plots
            figure(1); clf
            set(gcf,'Position',[0 100 1300 400]);
            ax1 = axes('Position',[0.05 0.1 0.28 0.8]); % glacier geometry
            hold on; grid on;
            set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
            title('Glacier Geometry'); legend('Location','northeast');
            xlim([0 65]); ylim([min(hb)-100 max(h)+200]);
            xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)');
            % ice surface
            plot(x(1:c)./10^3,h(1:c),'color',col(i,:),'linewidth',2,'displayname','0');
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
            set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
            title('Ice Speed'); legend('Location','northeast');
            xlim([0 65]); ylim([0 2500]);
            xlabel('Distance Along Centerline (km)'); ylabel('Speed (m yr^{-1})');
            % ice speed
            plot(x(1:c)./10^3,U(1:c).*3.1536e7,'color',col(i,:),'linewidth',2,'displayname','0');
            ax3 = axes('Position',[0.7 0.1 0.28 0.8]); % calving front position
            hold on; grid on;
            set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
            title('Calving Front Position'); legend('Location','best');
            xlim([30 65]); ylim([0 10]);
            xlabel('Distance Along Centerline (km)'); ylabel('Year');
            % calving front position
            plot(x(c)/10^3,t(i)./3.1536e7,'.','markersize',15,'color',col(i,:),'displayname','0');
        elseif mod(i-1,round(length(t)/50))==0 % display every length(t)/50
            figure(1);
            if mod(i-1,round(length(t)/10))==0 % display in legend every length(t)/10
                % ice surface
                plot(ax1,x(1:c)/10^3,h(1:c),'-','color',col(i,:),'linewidth',2,'displayname',num2str(t(i)./3.1536e7));
                % ice speed
                plot(ax2,x(1:c)/10^3,U(1:c).*3.1536e7,'-','Color',col(i,:),'linewidth',2,'DisplayName',num2str(t(i)./3.1536e7));
                % calving front position
                plot(ax3,x(c)./10^3,t(i)./3.1536e7,'.','Color',col(i,:),'markersize',15,'displayname',num2str(t(i)./3.1536e7)); hold on;
            else
                % ice surface
                plot(ax1,x(1:c)/10^3,h(1:c),'-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
                % ice speed
                plot(ax2,x(1:c)/10^3,U(1:c).*3.1536e7,'-','Color',col(i,:),'linewidth',2,'HandleVisibility','off'); hold on;
                % calving front position
                plot(ax3,x(c)./10^3,t(i)./3.1536e7,'.','Color',col(i,:),'markersize',15,'HandleVisibility','off'); hold on;
            end
            % calving front
            plot(ax1,x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
            % floating bed (gl:c)
            plot(ax1,x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off');
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
        [U,dUdx,vm,Un,dUndx,M,T] = U_convergence(x,U,dUdx,H,h,A,E,N,W,dx,c,n,m,beta,rho_i,rho_sw,g,sigma_b,i);
        
        % calculate ice flux
        F = U.*H.*W; % ice flux (m^3 s^-1)
        F(isnan(F))=0;
        
        % implement SMB & SMR
        smr = zeros(1,length(x));
        smr(gl+1:length(x)) = feval(fit([x(gl);x(gl+1);x(c)],[0;smr0;0.95*smr0],'pchip'),x(gl+1:end));
        smb = interp1(x0,smb0+Q0,x);
        
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
        
        % calculate RMSE of surface speed and elevation change at the
        % final model time
        if t(i)==t_end
            U_RMSE(f) = sqrt(nansum(interp1(x,U,x0)-U_obs(9).U').^2/length(x0));
            dH_tot(f,:) = interp1(x,h,x0)-h0;
            dH_RMSE(f) = sqrt((nansum(dH_tot(f,:)-dH_obs).^2)/length(x0));
        end        

        % calculate misfit in calving front position at each full model year
        if mod(t(i),3.1536e7)==0 && t(i)~=t_start
            c_misfit(f,t(i)./3.1536e7) = x(c)-interp1(termDate_obs,termx_obs,t(i)./3.1536e7+2009);
        end
        
    end 
    
end

% calculate mean misfit
misfit = mean(c_misfit,2);

% plot results for optimal sigma_b
optSigma_b = sigma_b0(abs(misfit)==min(abs(misfit)));
figure(10); clf; hold on;
    plot(sigma_b0./10^3,misfit,'.','markersize',20);
    plot(optSigma_b./10^3,misfit(abs(misfit)==min(abs(misfit))),'*','markersize',15,'linewidth',2);
    set(gca,'fontsize',14,'fontname','Arial','linewidth',2);
    grid on; xlabel('\sigma_b (kPa)'); ylabel('misfit (m)');
    title('Calving Front Position Misfit');

% save sigma_bbest
if save_optSigma_b
    cd([homepath,'inputs-outputs']);
    save('optimalSigma_b.mat','optSigma_b');
    disp(['Optimal sigma_b = ',num2str(optSigma_b./10^3),' kPa saved.']);
end

%% 4. Conduct sensitivity tests for SMB & SMR

% (1) Run through the the designated number of model years with no change 
%   in SMB or SMR and save the output. 
% (2) Run through simulated 2009-2019, then continue evolving the model  
%   under new scenarios: \Delta(SMB) &/or \Delta(SMR) and compare changes 
%   in geometry and speed with the no change scenario.  

close all; cd([homepath,'inputs-outputs/']);

save_figure = 1;    % = 1 to save resulting figure
save_final = 1;     % = 1 to save final geometry and speed 

% load no change conditions
load('Crane_2100_noChange.mat'); % load no change variables
    
% set up changes in SMB, SMR, & fwd
% note: decrease SMB & SMR in increments of 0.5 m a-1 (1.585e-8) starting
%   from 1 m a-1 until reaching max SMR found at other Antarctic ice shelves: 
%   ~6 m a^-1 (Adusumilli et al., 2020)
delta_smr = 0/3.1536e7; % m/s change in SMR
delta_smb = 0/3.1536e7; % m/s change in SMB
delta_fwd = -0.5; % m change in fwd
    
% define time stepping (s)
dt = 0.01*3.1536e7;
t_start = 0*3.1536e7;
t_end = 91*3.1536e7;
t = (t_start:dt:t_end);

% initialize variables
x=x0; H=H0; U=U0; W=W0; gl=gl0; dUdx=dUdx0; A=A0; h=h0; hb=hb0;

% load optimal parameters 
E=load('optimalE.mat').optE; % unitless
fwd = load('optimalfwd.mat').optfwd; % m 
sigma_b = load('optimalSigma_b.mat').optSigma_b; % Pa
    
% run flowline model
for i=1:length(t)
    
    % implement change to fwd after 10 model years
    if t(i)>10*3.1536e7
        % linearly increase fwd to reach delta_fwd by 2100
        delta_fwdi = delta_fwd/(2100-2019)*(t(i)/3.1536e7-10); % total increase in fwd at this time increment
        fwd = load('optimalfwd.mat').optfwd + delta_fwdi; % m
    end 
    
    % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
    Rxx = 2*nthroot(dUdx./(E.*A),n); % resistive stress (Pa)
    crev = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*fwd); % crevasse penetration depth (m)
    % calving front located where the inland-most crevasse intersects sea level
    xcf = interp1(h-crev,x,0,'linear','extrap'); % (m along centerline)
    if i==1 % use observed calving front for first iteration
        xcf = x0(c0);
    end
        
    % calculate the thickness required to remain grounded at each grid cell
    Hf = -(rho_sw./rho_i).*hb; % flotation thickness (m)
    % find the location of the grounding line and use a floating
    % geometry from the grounding line to the calving front
    if ~isempty(find(Hf-H>0,1,'first'))
        xgl = x(find(Hf-H>0,1,'first')-1);
        %xgl = interp1(Hf(find(Hf-H>0,1,'first')-1:find(Hf-H>0,1,'first')+1)...
        %    -H(find(Hf-H>0,1,'first')-1:find(Hf-H>0,1,'first')+1),...
        %    x(find(Hf-H>0,1,'first')-1:find(Hf-H>0,1,'first')+1),0,'linear','extrap'); % (m along centerline)
    else
        xgl=xcf;
    end
    if xgl>xcf % grounding line can't be past calving front
        xgl=xcf;
    elseif xgl<0 % grounding line can't be negative
        xgl=xcf;
    end

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
    E = interp1(x,E,xn,'linear','extrap');
    beta = interp1(x0,beta0,xn,'linear','extrap');
    x = xn; dx = dxn;

    % calculate surface elevation
    h = hb+H; % surface elevation (m a.s.l.)
    h(gl+1:c) = (1-rho_i/rho_sw).*H(gl+1:c); %adjust the surface elevation of ungrounded ice to account for buoyancy
    H(h<0)=0-hb(h<0); h(h<0)=0; % surface cannot go below sea level
    h(h-H<hb) = hb(h-H<hb)+H(h-H<hb); % thickness cannot go beneath bed elevation

    % plot geometry, speed, & calving front position
    if t(i)==t_start
        col = parula(length(t)+20); % color scheme for plots
        figure(1); clf
        set(gcf,'Position',[0 100 1300 400]);
        ax1 = axes('Position',[0.05 0.1 0.28 0.8]); % glacier geometry
        hold on; grid on;
        set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
        title('Glacier Geometry'); legend('Location','northeast');
        xlim([0 65]); ylim([min(hb)-100 max(h)+200]);
        xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)');
        % ice surface
        plot(x(1:c)./10^3,h(1:c),'color',col(i,:),'linewidth',2,'displayname','0');
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
        set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
        title('Ice Speed'); legend('Location','northeast');
        xlim([0 65]); ylim([0 2500]);
        xlabel('Distance Along Centerline (km)'); ylabel('Speed (m yr^{-1})');
        % ice speed
        plot(x(1:c)./10^3,U(1:c).*3.1536e7,'color',col(i,:),'linewidth',2,'displayname','0');
        ax3 = axes('Position',[0.7 0.1 0.28 0.8]); % calving front position
        hold on; grid on;
        set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
        title('Terminus & Grounding Line Positions'); legend('Location','best');
        xlim([30 65]); ylim([0 t_end./3.1536e7]);
        xlabel('Distance Along Centerline (km)'); ylabel('Year');
        % terminus & grounding line positions
        plot(x(c)/10^3,t(i)./3.1536e7,'.','markersize',15,'color',col(i,:),'displayname','0');
        plot(ax3,x(gl)./10^3,t(i)./3.1536e7,'x','Color',col(i,:),'markersize',10,'linewidth',2,'HandleVisibility','off');        
    elseif mod(i-1,round(length(t)/50))==0 % display every length(t)/50
        figure(1);
        if mod(i-1,round(length(t)/10))==0 % display in legend every length(t)/10
            % ice surface
            plot(ax1,x(1:c)/10^3,h(1:c),'-','color',col(i,:),'linewidth',2,'displayname',num2str(t(i)./3.1536e7));
            % ice speed
            plot(ax2,x(1:c)/10^3,U(1:c).*3.1536e7,'-','Color',col(i,:),'linewidth',2,'DisplayName',num2str(t(i)./3.1536e7));
            % calving front position
            plot(ax3,x(c)./10^3,t(i)./3.1536e7,'.','Color',col(i,:),'markersize',15,'displayname',num2str(t(i)./3.1536e7)); hold on;
        else
            % ice surface
            plot(ax1,x(1:c)/10^3,h(1:c),'-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
            % ice speed
            plot(ax2,x(1:c)/10^3,U(1:c).*3.1536e7,'-','Color',col(i,:),'linewidth',2,'HandleVisibility','off'); hold on;
            % calving front & grounding line positions
            plot(ax3,x(c)./10^3,t(i)./3.1536e7,'.','Color',col(i,:),'markersize',15,'HandleVisibility','off'); hold on;
            plot(ax3,x(gl)./10^3,t(i)./3.1536e7,'x','Color',col(i,:),'markersize',10,'linewidth',2,'HandleVisibility','off');
        end
        % calving front
        plot(ax1,x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
        % floating bed (gl:c)
        plot(ax1,x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off');
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
    [U,dUdx,vm,Un,dUndx,M,T] = U_convergence(x,U,dUdx,H,h,A,E,N,W,dx,c,n,m,beta,rho_i,rho_sw,g,sigma_b,i);

    % calculate ice flux
    F = U.*H.*W; % ice flux (m^3 s^-1)
    F(isnan(F))=0;

    % implement SMB, SMR, delta_SMB, & dela_SMR
    if t(i)/3.1536e7<10 % use original SMB & SMR for first 10 model years
        smr = zeros(1,c);
        smr(gl+1:length(x)) = feval(fit([x(gl);x(gl+1);x(c)],[0;smr0;0.95*smr0],'pchip'),x(gl+1:end));
        smb = interp1(x0,smb0+Q0,x);  
        % plot
        if i==1
            figure(2); clf
            subplot(1,2,1); hold on; 
                set(gca,'fontsize',12,'fontweight','bold','linewidth',2);
                plot(x,smb.*3.1536e7,'color',col(i,:),'linewidth',2);
                xlabel('m along centerline'); ylabel('m a^{-1}'); grid on;
                title('SMB');
            subplot(1,2,2); hold on; 
                set(gca,'fontsize',12,'fontweight','bold','linewidth',2);
                plot(x,smr.*3.1536e7,'color',col(i,:),'linewidth',2);
                xlabel('m along centerline'); ylabel('m a^{-1}'); grid on;
                title('SMR');
        end        
    elseif t(i)/3.1536e7>=10 % implement changes after 10 model years
        smr = zeros(1,c);
        % increase SMR linearly at each time increment to reach delta_smr by 2100
        delta_smri = delta_smr/(2100-2019)*(t(i)/3.1536e7-10); % total increase in smr from 2019 rate
        smr(gl+1:length(x)) = feval(fit([x(gl);x(gl+1);x(c)],[0;smr0+delta_smri;0.95*smr0+delta_smri],'pchip'),x(gl+1:end));
        delta_smbi = delta_smb/(2100-2019)*(t(i)/3.1536e7-10); % total increase in smr from 2019 rate        
        smb = interp1(x0,smb0+Q0+delta_smbi,x);
        % plot
        if mod(t(i)/3.1536e7,10)==0
            figure(2);
            subplot(1,2,1);
                plot(x,smb.*3.1536e7,'color',col(i,:),'linewidth',2);
            subplot(1,2,2);
                plot(x,smr.*3.1536e7,'color',col(i,:),'linewidth',2);
        end
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

% plot results
if t(i)==t_end
    h2=h; H2=H; x2=x; c2=c; gl2=gl; U2=U; fwd2=fwd; % store final geometry & speed
    figure(10); clf % sensitivity test changes
    hold on; grid on;
    set(gcf,'Position',[491 80 886 686]);
    set(gca,'FontSize',14,'linewidth',2,'fontweight','bold');
    xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)');
    legend('Location','east'); xlim([0 70]); ylim([-1200 max(h)+100]);
    title(['SMR = + ',num2str(round(delta_smr.*3.1536e7,1)),'m/a, SMB = + ',...
        num2str(round(delta_smb*3.1536e7,1)),'m/a, fwd = ',num2str(fwd),'m']);
    ax1=get(gca);
        % ice surface
        plot(x1(1:c1)/10^3,h1(1:c1),'-k','linewidth',2,'displayname','no change');
        plot(x2(1:c2)/10^3,h2(1:c2),'color',[0.8 0 0],'linewidth',2,'displayname','change');
        % calving front
        plot(x1(c1)*[1,1]/10^3,[h1(c1)-H1(c1),h1(c1)],'-k','linewidth',2,'HandleVisibility','off');
        plot(x2(c2)*[1,1]/10^3,[h2(c2)-H2(c2),h2(c2)],'color',[0.8 0 0],'linewidth',2,'HandleVisibility','off');
        % floating bed
        plot(x1(gl1:c1)/10^3,h1(gl1:c1)-H1(gl1:c1),'-k','linewidth',2,'HandleVisibility','off');
        plot(x2(gl2:c2)/10^3,h2(gl2:c2)-H2(gl2:c2),'color',[0.8 0 0],'linewidth',2,'HandleVisibility','off');
        % bed elevation
        plot(x1/10^3,hb1,'k','linewidth',2,'HandleVisibility','off');
    % inset plot of terminus
    ax2 = axes('Position',[0.62 0.62 0.28 0.28]);
        hold on; grid on; set(gca,'linewidth',2,'fontweight','bold','fontsize',9);
        delta_L = x2(c2)-x1(c1); % change in length (m)
        delta_H = mean(H2(1:c2))-mean(H1(1:c1)); % change in ice thickness (m)
        delta_U = mean(U2(1:c2))-mean(U1(1:c1)); % change in ice speed (m/s) 
        title(['\Delta L=',num2str(round(delta_L)),'m, ','\DeltaH_{\mu}=',...
            num2str(round(delta_H)),' m, ','\DeltaU_{\mu}=',...
            num2str(round(delta_U.*3.1536e7)),' m/a']);
        xlim([40 60]); ylim([-1000 300]);
        % ice surface
        plot(x1(1:c1)/10^3,h1(1:c1),'-k','linewidth',2,'displayname','no change');
        plot(x2(1:c2)/10^3,h2(1:c2),'color',[0.8 0 0],'linewidth',2,'displayname','no change');
        % calving front
        plot(x1(c1)*[1,1]/10^3,[h1(c1)-H1(c1),h1(c1)],'-k','linewidth',2,'HandleVisibility','off');
        plot(x2(c2)*[1,1]/10^3,[h2(c2)-H2(c2),h2(c2)],'color',[0.8 0 0],'linewidth',2,'HandleVisibility','off');
        % floating bed
        plot(x1(gl1:c1)/10^3,h1(gl1:c1)-H1(gl1:c1),'-k','linewidth',2,'HandleVisibility','off');
        plot(x2(gl2:c2)/10^3,h2(gl2:c2)-H2(gl2:c2),'color',[0.8 0 0],'linewidth',2,'HandleVisibility','off');
        % bed elevation
        plot(x/10^3,hb,'k','linewidth',2,'HandleVisibility','off');
        % mean sea level
        plot([x(1),x(end)]/10^3,[0,0],'k--','HandleVisibility','off');
end

% save figure
if save_figure && ishandle(10)
    cd([homepath,'scripts/modelingWorkflow/results/']);
    % Save figure for test
    if delta_fwd==0
        fileName = ['SMR',num2str(delta_smr),'_SMB',num2str(delta_smb)];
    else
        fileName = ['fwd_',num2str(fwd),'m'];
    end
    saveas(gcf,[fileName,'.png'],'png');
    disp(['Fig. 10 saved in: ',pwd]);
elseif ~ishandle(4)
    disp('figure does not exist.');
else
    disp('no figure saved.');
end

% save geometry & speed
if save_final
    cd([homepath,'scripts/modelingWorkflow/results2']);
    if delta_fwd==0 && delta_smr==0 && delta_smb==0 
        save(['SMR0_SMB0_geom.mat'],'h2','H2','c2','c2','U2','gl2','x2','fwd2');
        save(['fwd_',num2str(fwd),'_geom.mat'],'h2','H2','c2','c2','U2','gl2','x2','fwd2');
    elseif delta_fwd==0 
        fileName = ['SMR',num2str(delta_smr),'_SMB',num2str(delta_smb),'_geom.mat'];
        save(fileName,'h2','H2','c2','c2','U2','gl2','x2','fwd2');
    else
        fileName = ['fwd_',num2str(fwd),'m_geom.mat'];
        save(fileName,'h2','H2','c2','c2','U2','gl2','x2','fwd2');
    end
    disp('geometry saved.');
else
    disp('geometry not saved.');
end

