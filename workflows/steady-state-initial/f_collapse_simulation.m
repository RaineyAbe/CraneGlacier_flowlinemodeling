%% Glacier flowline model: script to simulate glacier response to ice shelf collapse
% Rainey Aberle and Ellyn Enderlin
% Spring 2022
% 
%   0. Load initial model parameters
%   1. Set up future climate scenarios
%   2. Model pre-collapse conditions
%   3. Model ice shelf collapse

% ------------------------------------------
% -----0. Load initial model parameters-----  
% ------------------------------------------

clear all; close all;
warning off; % turn off warnings (velocity coefficient matrix is close to singular)
    
plotTimeSteps = 0; % = 1 to plot geometry, speed, cf/gl positions throughout model time period
plotClimateParams = 0; % = 1 to plot SMB, DFW, TF over time
saveFinal = 1; % = 1 to save pre-collapse and final (2100) conditions

% define home path in directory and add necessary paths
homepath = '/Users/raineyaberle/Research/MS/CraneGlacier_flowlinemodeling/';
addpath([homepath,'workflows/'],...
    [homepath,'workflows/2_steady-state/'],...
    [homepath,'functions/cmocean_v2.0/cmocean'],...
    [homepath,'inputs-outputs/']);

% -----load initialization file
load([homepath,'inputs-outputs/modelInitialization_preCollapse.mat']);
SMR_mean_fit = load('LarsenC_MeanMeltRate.mat').mr_mean_fit;

% -----time stepping [s]-----
dt = 0.001*3.1536e7;

% -----densities and g-----
rho_i = 917; % ice density (kg m^-3)
rho_sw = 1028; % ocean water density (kg m^-3)
rho_fw = 1000; % fresh water density (kg m^-3)
g = 9.81; % acceleration (m s^-2)

% -----stress parameters (unitless)-----
m = 3; % basal sliding exponent
n = 3; % flow law exponent
E = 1; % enhancement factor

% -----calving parameters-----
DFW0 = 0; % fresh water depth in crevasses [m]

% -----instability checks (using thickness and speed)
H_max = 2000; % maximum thickness (m)
H_min = 100;  % minimum thickness (m)
U_min = 100./3.1536e7;  % minimum mean speed (m s^-1)
    
% -----initial conditions-----
dx0 = mean(x0(2:end)-x0(1:end-1));
H0 = h0-b0; % ice thickness (m)
dUdx0 = [(U0(2:end)-U0(1:end-1))./(x0(2:end)-x0(1:end-1)) 0]; % strain rate (1/s) %%EE: flipped this to be consistent w/ other computations 15/09/21
% find the location of the grounding line and the end of the ice-covered domain
Hf = -(rho_sw./rho_i).*b0; % flotation thickness (m)
gl0 = find(Hf-H0>0,1,'first')-1; % grounding line location 
if isempty(gl0)
    % set grounding line to calving front location if all ice is grounded
    gl0=c0; 
else
    % adjust thickness for flotation
    H0(gl0+1:length(x0))=h0(gl0+1:end)*rho_sw/(rho_sw-rho_i); % buoyant thickness using surface
    H0(H0>=(h0-b0))=h0(H0>=(h0-b0))-b0(H0>=(h0-b0)); % thickness can't go beneath bed elevation
    H0(c0+1:end) = 0; % zero ice thickness past calving front 
end

% -----observed calving front positions for comparison
term = load([homepath,'inputs-outputs/terminusPositions_2002-2019.mat']).term;

% --------------------------------------------
% -----1. Set up future climate scenarios-----
% --------------------------------------------

% -----tune TF-SMR relationship to observations
% TF-SMR Eqn from Rignot et al. (2016) and Xu et al. (2013):
% mdot = (C * -b0(gl0) * sum(RO(1:gl))^0.39 + 0.15)*TF0^1.18
% Use initial melt rate estimates from Dryak and Enderlin (2020) to solve
% for C. rearrange Eqn.:
% C= (mdot - 0.15*TF^1.18) / (-b0(gl0) * sum(RO(1:gl))^0.39 * TF^1.18)
% Note: mdot and RO must be in units of m/d.
% initial total runoff = sum(RO0(x) * W(x) * dx(x))
% terminus area-averaged runoff = initial total runoff / (H(c) * W(c))
RO0_sum = sum(RO0(1:c0).*W0(1:c0).*dx0) / (H0(c0) * W0(c0)); % m/s
TF0 = 0.2; % ocean thermal forcing [^oC], estimated from Larsen B icebergs
mdot0 = -SMR0; % m/s
% C = (mdot0*86400 - 0.15*nthroot(TF0,1/1.18)) / (-b0(c0) * nthroot(RO0_sum*86400,1/0.39) * nthroot(TF0,1/1.18)); % <-- results in a negative number
C = (mdot0*86400 - 0.15*nthroot(TF0,1/2.3)) / (-b0(c0) * nthroot(RO0_sum*86400,1/0.39) * nthroot(TF0,1/2.3)); % <-- adjusted TF exponent until C=positive

% -----define potential changes in SMB, DFW, & TF
% - decrease maximum SMB by increments of 0.5 m a-1 down to -10 m a-1
% - increase DFW by increments of 1 m up to 10 m
% - increase TF by increments of 0.1 ^oC up to 1 ^oC
delta_SMB0 = (0:-1:-10)./3.1536e7; % m/s change in SMB at the calving front (used to increase gradient)
delta_DFW0 = 0:10; % m change in DFW 
delta_TF0 = 0:0.1:1; % ^oC change in TF

% -----store mean final SMB for plotting
smb_mean = NaN*zeros(1,length(delta_SMB0));
dSMR_max = 0; 

for j=1:length(delta_SMB0)
    
    % Switch scenarios on and off
    delta_SMB = 0;% delta_SMB0(j);
    delta_DFW = 0; %delta_DFW0(j);
    delta_TF = delta_TF0(j);
    SMB_enhance = 0; % = 1 to increase SMR due to decreased SMB    

    % ------------------------------------------
    % -----2. Model pre-collapse conditions-----
    % ------------------------------------------

    disp('Simulating pre-collapse, steady-state conditions');

    % -----initialize parameters-----
    x=x0; U=U0; W=W0; gl=gl0; dUdx=dUdx0; A=A0; h=h0; b=b0; H=H0; 
    DFW=DFW0; dx=dx0; SMB=SMB0; SMR=SMR0; c=c0;
    beta0 = interp1([0 x0(end)], [1 2], x0);  
    sigma_b = 800e3; 

    col = parula(10e3); % color scheme for plots

    % -----run flowline model-----
    i=1; % counter for iterations
    while i

        % -----establish time
        t(i) = (i-1)*dt; % [s]

        % -----calving front location
        c = c; % constant until steady-state is achieved
        xcf = x(c); % calving front location [m]

        % -----grounding line location
        % calculate the thickness required to remain grounded at each grid cell
        Hf = -(rho_sw./rho_i).*b; % flotation thickness (m)
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

        % -----create coordinate system that hits cf and gl exactly
        % has resolution dxmax near the ice divide
        % has resolution dxmin from gl to c
        % and has smooth variation between
        xl = round(xgl/dx0); % number of ideal grid spaces needed to reach the grounding line
        dx = xgl/xl; % new grid spacing (should be ~dx0)
        xn = 0:dx:xgl; % new distance vector
        if xcf-xgl > 0
            xl = round((xcf-xgl)/dx0);
            dx = (xcf-xgl)/xl;
            xn = [xn xn(end)+dx:dx:xcf];
        end
        clear dx; dxn = [xn(2:end)-xn(1:end-1) xn(end)-xn(end-1)];

        % -----get geometry on new coordinates
        c = length(xn);
        H = interp1(x,H,xn,'linear','extrap');
        % H gradient
        dHdx(1) = (H(2)-H(1))/(xn(2)-xn(1)); % forward difference
        dHdx(2:c-1) = (H(3:c)-H(1:c-2))./(xn(3:c)-xn(1:c-2)); % central difference
        dHdx(c) = (H(c)-H(c-1))./(xn(c)-xn(c-1)); % backward difference            
        % if there is a sudden jump in H (large gradient) past the grounding
        % line, set that to the previous point
        % thickness gradient
        if any(dHdx>50)
            c = find(dHdx>50,1,'first')-1;            
            H = H(1:c);
            xn = xn(1:c);
        end
        gl = dsearchn(xn',xgl); % indices for xcf and xgl
        b = interp1(x0,b0,xn,'linear','extrap');
        W = interp1(x0,W0,xn,'linear','extrap');
        U = interp1(x,U,xn,'linear','extrap');
        A = interp1(x0,A0,xn,'linear','extrap');
        beta = interp1(x0,beta0,xn,'linear','extrap');
        x = xn; dx = dxn; clear xn dxn; %EE: added the clear statement 14/09/21
        XGL(i) = xgl; % save grounding line position over time

        % -----calculate surface elevation
        h = b+H; % surface elevation (m a.s.l.)
        h(gl+1:c) = (1-rho_i/rho_sw).*H(gl+1:c); % adjust the surface elevation of ungrounded ice to account for buoyancy
        H(h<0)=0-b(h<0); h(h<0)=0; % surface cannot go below sea level
        h(h-H<b) = b(h-H<b)+H(h-H<b); % thickness cannot go beneath bed elevation

        % -----plot geometry, speed, & grounding line and calving front positions
        if plotTimeSteps
            if i==1
                drawnow
                figure(1); clf
                set(gcf,'Position',[-1500 100 1300 450]);
                ax1 = axes('Position',[0.05 0.12 0.27 0.78]); % glacier geometry
                hold on; grid on;
                set(gca,'FontSize',14,'linewidth',2);
                title([num2str(t(i)*3.1536e7), 'yrs']);
                xlim([0 95]); ylim([min(b)-100 max(h)+200]);
                xlabel('distance along centerline [km]'); ylabel('elevation [m]');
                % ice surface
                plot(x(1:c)./10^3,h(1:c),'color',col(i,:),'linewidth',2,'displayname','2009');
                % calving front
                plot(x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
                % floating bed
                plot(x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off');
                % bed elevation
                plot(x0./10^3,b0,'k','linewidth',2,'HandleVisibility','off');
                % mean sea level
                plot([x(1),x(end)]/10^3,[0,0],'k--','HandleVisibility','off');
                % ice speed
                ax2 = axes('Position',[0.37 0.1 0.28 0.8]);
                hold on; grid on;
                set(gca,'FontSize',14,'linewidth',2);
                xlim([0 95]); 
    %                     ylim([0 800]);
                xlabel('distance along centerline [km]'); ylabel('speed [m a^{-1}]');
                plot(x(1:c)./10^3,U(1:c).*3.1536e7,'color',col(i,:),'linewidth',2,'displayname','2009');
                % calving front & grounding line positions
                ax3 = axes('Position',[0.7 0.1 0.28 0.8]); 
                hold on; grid on;
                set(gca,'FontSize',14,'linewidth',2);
                xlim([30 95]);
                xlabel('distance along centerline [km]'); 
                ylabel('Year');
                plot(x(c)/10^3,t(i)./3.1536e7,'.','markersize',15,'color',col(i,:),'displayname','2009');
                plot(ax3,x(gl)./10^3,t(i)./3.1536e7,'x','Color',col(i,:),'markersize',10,'linewidth',2,'HandleVisibility','off');
            elseif mod(t(i),dt*200)==0 % display every dt*200
                figure(1);
                % glacier geometry
                plot(ax1,x(1:c)/10^3,h(1:c),'-','color',col(i,:),'linewidth',2,'displayname',num2str(round(t(i)./3.1536e7)+2009));
                plot(ax1,x(gl:c)/10^3,h(gl:c)-H(gl:c),'-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
                plot(ax1,[x(c);x(c)]/10^3,[h(c);h(c)-H(c)],'-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
                title(ax1, [num2str(t(i)/3.1536e7), 'yrs']);
                plot(ax1,x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off'); % calving front
                plot(ax1,x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off'); % floating bed (gl:c)
                % ice speed
                plot(ax2,x(1:c)/10^3,U(1:c).*3.1536e7,'-','Color',col(i,:),'linewidth',2,'DisplayName',num2str(round(t(i)./3.1536e7)+2009));
                % calving front & grounding line positions
                plot(ax3,x(c)/10^3,t(i)/3.1536e7,'.','Color',col(i,:),'markersize',15,'displayname',num2str(round(t(i)./3.1536e7)+2009)); hold on;
                plot(ax3,x(gl)/10^3,t(i)/3.1536e7,'x','Color',col(i,:),'markersize',10,'linewidth',2,'HandleVisibility','off'); hold on;
            end    
        end

        % -----calculate the effective pressure 
        % (ice overburden pressure minus water pressure) assuming an easy & 
        % open connection between the ocean and ice-bed interface
        sl = find(b<=0,1,'first'); % find where the glacier base first drops below sea level
        N_ground = rho_i*g*H(1:sl); % effective pressure where the bed is above sea level (Pa)
        N_marine = rho_i*g*H(sl+1:length(x))+(rho_sw*g*b(sl+1:length(x))); % effective pressure where the bed is below sea level (Pa)
        N = [N_ground N_marine];
        N(N<0)=0; % cannot have negative values

        % -----solve for new velocity
        [U,dUdx] = U_convergence(x,U,H,h,A,E,N,W,dx,c,n,m,beta,rho_i,rho_sw,g,sigma_b);

        % -----calculate ice flux
        F = U.*H.*W; % ice flux (m^3 s^-1)
        F(isnan(F))=0;
        F(1)=F(2);

        % -----save grounding line oce flux (discharge)
        Fgl(i) = F(gl)*917*1e-12*3.1536e7; % Gt/a

        % -----implement SMB, SMR, & RO
        SMB = interp1(x0,SMB0+Q0,x);
        RO = interp1(x0,RO0,x);
        SMR = zeros(1,c);
        % use the Larsen C mean melt rate profile to scale SMR
        % using the max initial SMR
        if gl<c
            SMR(gl+1:c) = SMR0/(SMR_mean_fit.a+1)*feval(SMR_mean_fit,x(gl+1:c)-x(gl+1)); 
        end

        % -----calculate the  change in ice thickness from continuity
        clearvars dHdt
        dHdt(1) = (-1/W(1))*(F(1)-F(2))/(x(1)-x(2)); % forward difference
        dHdt(2:c-1) = (-1./W(2:c-1)).*(F(1:c-2)-F(3:c))./(x(1:c-2)-x(3:c)); % central difference
        dHdt(c:length(x)) = (-1./W(c:length(x))).*(F(c-1:length(x)-1)-F(c:length(x)))./(x(c-1:length(x)-1)-x(c:length(x))); % backward difference
        dH = dHdt.*dt;

        % -----new thickness (change from dynamics, SMB, & SMR)
        Hn = H+dH+(SMB.*dt)+(SMR.*dt)-(interp1(x0,RO0,x)*dt);
        Hn(Hn < 0) = 0; % remove negative values

        % -----stop the model if it behaves unstably (monitored by ice thickness and speed)
        if max(H) > H_max
            disp(['Adjust dt']);
            break;
        end
        if mean(U) < U_min
            disp('Too slow!');
            break;
        end
        if any(~isfinite(H(1:c))) || any(~isfinite(U(1:c))) || any(~isfinite(h(1:c)))
            disp('non finite values');
            break;
        end

        % -----stop model if stead-state conditions reached
        % (change in U at each point is less than set threshold) 
        if all(abs(H-Hn) < 0.00002*abs(H))
            disp('    steady-state conditions achieved'); 
            disp('    continuing...');             
            break;
        else
            H = Hn; 
        end

        % -----continue loop
        i=i+1;
    end

    % -----save
    % reassign variable names
    xi=x; hi=h; bi=b; Hi=H; Ui=U; dUdxi=dUdx; Wi=W; Ai=A; betai=beta; gli=gl; ci=c; 
    sigma_bi=sigma_b; SMBi=SMB; SMRi=SMR; DFWi=DFW; 
    if saveFinal
        save([homepath,'inputs-outputs/steady_state_conditions.mat'],'x', 'h', 'b', 'H',...
            'U', 'dUdx','W', 'A', 'beta', 'gl', 'c', 'sigma_b', 'SMB', 'SMR', 'DFW');
        disp('steady-state conditions saved');
    end

    % -----plot surface crevasses
    % resistive stress (Pa)
    Rxx = 2*nthroot(dUdx./(E.*A),n); 
    % height above buoyancy (m)
    Hab = H+rho_sw/rho_i*(b); 
    % surface crevasse penetration depth (m)
    crev = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*(DFW)); 
    % basal crevasse penetration depth (m)
    crev_b = rho_i/(rho_sw-rho_i).*(Rxx./(rho_i*g)-Hab);     

    figure(10); clf
    subplot(1,2,1); hold on; 
    set(gca,'fontsize',12,'linewidth',1);
    legend('location','southeast');
    grid on;
    plot(x/10^3, h,'-b','linewidth',2,'displayname','h');
    plot(x/10^3, crev,'-m','linewidth',2,'displayname','crev_s');
    plot(x/10^3, crev_b,'-c','linewidth',2,'displayname','crev_b');
    xlabel('distance along centerline [km]');
    ylabel('elevation [m]');
    subplot(1,2,2); hold on; 
    set(gca,'fontsize',12,'linewidth',1);
    legend('location','southwest');
    grid on;
    plot(x/10^3, h-crev_b,'-c','linewidth',2,'displayname','h - crev_b');
    plot(x/10^3, h-crev,'-m','linewidth',2,'displayname','h - crev_s');
    xlabel('distance along centerline [km]');
    ylabel('elevation [m]');
    
    % -------------------------------------
    % -----3. Model ice shelf collapse-----
    % -------------------------------------

    disp('Simulating ice shelf collapse');
    
    col=parula(1000); % color scheme for plotting

    % ----------A. BACKSTRESS REMOVAL----------

    % -----remove backstress
    sigma_b = 0; % Pa 

    % -----run flowline model until Rxx = extensive past the grounding line
    i=1; clear t; % counter for iterations
    while i

        % -----establish time
        t(i) = (i-1)*dt; % [s]

        % -----calving front location 
        c=length(x); % constant until resistive stress criterion met
        % along-flow resistive stress [Pa]
        Rxx = 2*nthroot(dUdx./(E.*A),n);    

        % -----plot geometry, speed, & grounding line and calving front
        % positions, Rxx
        if plotTimeSteps
            if i==1
                drawnow
                figure(1); clf
                ax1 = axes('Position',[0.05 0.12 0.27 0.78]); % glacier geometry
                hold on; grid on;
                set(gca,'FontSize',12,'linewidth',2);
                title([num2str(t(i)*3.1536e7), 'yrs']);
                xlim([0 95]); ylim([min(b)-100 max(h)+200]);
                xlabel('distance along centerline [km]'); ylabel('elevation [m]');
                % ice surface
                plot(x(1:c)./10^3,h(1:c),'color',col(i,:),'linewidth',2,'displayname','2009');
                % calving front
                plot(x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
                % floating bed
                plot(x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off');
                % bed elevation
                plot(x0./10^3,b0,'k','linewidth',2,'HandleVisibility','off');
                % mean sea level
                plot([x(1),x(end)]/10^3,[0,0],'k--','HandleVisibility','off');
                % ice speed
                ax2 = axes('Position',[0.37 0.1 0.28 0.8]);
                hold on; grid on;
                set(gca,'FontSize',12,'linewidth',2);
                xlim([0 95]); 
                xlabel('distance along centerline [km]'); ylabel('speed [m a^{-1}]');
                plot(x(1:c)./10^3,U(1:c).*3.1536e7,'color',col(i,:),'linewidth',2,'displayname','2009');
                % calving front and grounding line position
                ax3 = axes('Position',[0.7 0.1 0.28 0.8]); 
                hold on; grid on;
                set(gca,'FontSize',12,'linewidth',2);
                xlabel('distance along centerline [km]'); 
                ylabel('year');
                plot(ax3,x(c)/10^3,t(i)/3.1536e7,'.','Color',col(i,:),'markersize',10,'linewidth',2,'HandleVisibility','off'); hold on;
                plot(ax3,x(gl)/10^3,t(i)/3.1536e7,'x','Color',col(i,:),'markersize',10,'linewidth',2,'HandleVisibility','off'); hold on;            
                % Rxx
                figure(3); clf; 
                ax4 = gca;
                hold on; grid on;
                set(ax4,'FontSize',12,'linewidth',2);
                xlabel('distance along centerline [km]'); 
                ylabel('R_{xx} [kPa]');
                plot(ax4, x/10^3, Rxx/10^3, 'Color',col(i,:),'markersize',10,'linewidth',2);
            else
                figure(1);
                % glacier geometry
                plot(ax1,x(1:c)/10^3,h(1:c),'-','color',col(i,:),'linewidth',2,'displayname',num2str(round(t(i)./3.1536e7)+2009));
                plot(ax1,x(gl:c)/10^3,h(gl:c)-H(gl:c),'-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
                plot(ax1,[x(c);x(c)]/10^3,[h(c);h(c)-H(c)],'-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
                title(ax1, [num2str(t(i)/3.1536e7), 'yrs']);
                plot(ax1,x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off'); % calving front
                plot(ax1,x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off'); % floating bed (gl:c)
                % ice speed
                plot(ax2,x(1:c)/10^3,U(1:c).*3.1536e7,'-','Color',col(i,:),'linewidth',2,'DisplayName',num2str(round(t(i)./3.1536e7)+2009));
                % calving front & grounding line positions
                plot(ax3,x(c)/10^3,t(i)/3.1536e7,'.','Color',col(i,:),'markersize',10,'linewidth',2,'HandleVisibility','off'); 
                plot(ax3,x(gl)/10^3,t(i)/3.1536e7,'x','Color',col(i,:),'markersize',10,'linewidth',2,'HandleVisibility','off'); 
                % Rxx
                plot(ax4, x/10^3, Rxx/10^3, 'Color',col(i,:),'markersize',10,'linewidth',2);
            end    
        end

        % -----stop model when resistive stress is extensive near calving front (> 0)
        if any(Rxx(gl:c)>=0)
            disp('    extensive conditions reached');
            disp('    continuing...');
            if plotTimeSteps
                % plot
                figure(1);
                % glacier geometry
                plot(ax1,x(1:c)/10^3,h(1:c),'-m','linewidth',2,'displayname',num2str(round(t(i)./3.1536e7)+2009));
                plot(ax1,x(gl:c)/10^3,h(gl:c)-H(gl:c),'-m','linewidth',2,'HandleVisibility','off');
                plot(ax1,[x(c);x(c)]/10^3,[h(c);h(c)-H(c)],'-m','linewidth',2,'HandleVisibility','off');
                title(ax1, [num2str(t(i)/3.1536e7), 'yrs']);
                plot(ax1,x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'-m','linewidth',2,'HandleVisibility','off'); % calving front
                plot(ax1,x(gl:c)/10^3,h(gl:c)-H(gl:c),'-m','linewidth',2,'HandleVisibility','off'); % floating bed (gl:c)
                % ice speed
                plot(ax2,x(1:c)/10^3,U(1:c).*3.1536e7,'-m','linewidth',2,'DisplayName',num2str(round(t(i)./3.1536e7)+2009));
                % calving front & grounding line positions
                plot(ax3,x(c)/10^3,t(i)/3.1536e7,'.m','markersize',10,'linewidth',2,'HandleVisibility','off'); 
                plot(ax3,x(gl)/10^3,t(i)/3.1536e7,'xm','markersize',10,'linewidth',2,'HandleVisibility','off'); 
                % Rxx
                plot(ax4, x/10^3, Rxx/10^3, '-m','markersize',10,'linewidth',2);
            end
            break;
        end

        % -----grounding line location
        % calculate the thickness required to remain grounded at each grid cell
        Hf = -(rho_sw./rho_i).*b; % flotation thickness (m)
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

        % -----create coordinate system that hits gl exactly
        xl = round(xgl/dx0); % number of ideal grid spaces needed to reach the grounding line
        dx = xgl/xl; % new grid spacing (should be ~dx0)
        xn = 0:dx:xgl; % new distance vector
        if xcf-xgl > 0
            xl = round((xcf-xgl)/dx0);
            dx = (xcf-xgl)/xl;
            xn = [xn xn(end)+dx:dx:xcf];
        end
        clear dx; dxn = [xn(2:end)-xn(1:end-1) xn(end)-xn(end-1)];

        % -----get geometry on new coordinates
        c = length(xn); % index for xcf
        H = interp1(x,H,xn,'linear','extrap');
        % H gradient
        dHdx(1) = (H(2)-H(1))/(xn(2)-xn(1)); % forward difference
        dHdx(2:c-1) = (H(3:c)-H(1:c-2))./(xn(3:c)-xn(1:c-2)); % central difference
        dHdx(c) = (H(c)-H(c-1))./(xn(c)-xn(c-1)); % backward difference            
        % if there is a sudden jump in H (large gradient) past the grounding
        % line, set that to the previous point
        % thickness gradient
        if any(dHdx>50)
            c = find(dHdx>50,1,'first')-1;            
            H = H(1:c);
            xn = xn(1:c);
        end
        gl = dsearchn(xn',xgl); % index for xgl
        b = interp1(x0,b0,xn,'linear','extrap');
        W = interp1(x0,W0,xn,'linear','extrap');
        U = interp1(x,U,xn,'linear','extrap');
        A = interp1(x0,A0,xn,'linear','extrap');
        beta = interp1(x0,beta0,xn,'linear','extrap');
        x = xn; dx = dxn; 
        clear xn dxn; % clear to prevent issues in future iterations
        XGL(i) = xgl; % save grounding line position over time

        % -----calculate surface elevation
        h = b+H; % surface elevation (m a.s.l.)
        h(gl+1:c) = (1-rho_i/rho_sw).*H(gl+1:c); % adjust the surface elevation of ungrounded ice to account for buoyancy
        H(h<0)=0-b(h<0); h(h<0)=0; % surface cannot go below sea level
        h(h-H<b) = b(h-H<b)+H(h-H<b); % thickness cannot go beneath bed elevation

        % -----calculate the effective pressure 
        % (ice overburden pressure minus water pressure) assuming an easy & 
        % open connection between the ocean and ice-bed interface
        sl = find(b<=0,1,'first'); % find where the glacier base first drops below sea level
        N_ground = rho_i*g*H(1:sl); % effective pressure where the bed is above sea level (Pa)
        N_marine = rho_i*g*H(sl+1:length(x))+(rho_sw*g*b(sl+1:length(x))); % effective pressure where the bed is below sea level (Pa)
        N = [N_ground N_marine];
        N(N<0)=0; % cannot have negative values

        % -----solve for new velocity
        [U,dUdx] = U_convergence(x,U,H,h,A,E,N,W,dx,c,n,m,beta,rho_i,rho_sw,g,sigma_b);

        % -----calculate ice flux
        F = U.*H.*W; % ice flux (m^3 s^-1)
        F(isnan(F))=0;
        F(1)=F(2);

        % -----save grounding line oce flux (discharge)
        Fgl(i) = F(gl)*917*1e-12*3.1536e7; % Gt/a

        % -----implement SMB, SMR, & RO
        SMB = interp1(x0,SMB0+Q0,x);
        RO = interp1(x0,RO0,x);
        SMR = zeros(1,c);
        % use the Larsen C mean melt rate profile to scale SMR
        % using the max initial SMR
        if gl<c
            SMR(gl+1:c) = SMR0/(SMR_mean_fit.a+1)*feval(SMR_mean_fit,x(gl+1:c)-x(gl+1)); 
        end

        % -----calculate the  change in ice thickness from continuity
        clearvars dHdt
        dHdt(1) = (-1/W(1))*(F(1)-F(2))/(x(1)-x(2)); % forward difference
        dHdt(2:c-1) = (-1./W(2:c-1)).*(F(1:c-2)-F(3:c))./(x(1:c-2)-x(3:c)); % central difference
        dHdt(c:length(x)) = (-1./W(c:length(x))).*(F(c-1:length(x)-1)-F(c:length(x)))./(x(c-1:length(x)-1)-x(c:length(x))); % backward difference
        dH = dHdt.*dt;

        % -----new thickness (change from dynamics, SMB, & SMR)
        Hn = H+dH+(SMB.*dt)+(SMR.*dt)-(interp1(x0,RO0,x)*dt);
        Hn(Hn < 0) = 0; % remove negative values
        H = Hn; % set thickness to new value

        % -----stop the model if it behaves unstably (monitored by ice thickness and speed)
        if max(H) > H_max
            disp(['Adjust dt']);
            break;
        end
        if mean(U) < U_min
            disp('Too slow!');
            %break;
        end
        if any(~isfinite(H(1:c))) || any(~isfinite(U(1:c))) || any(~isfinite(h(1:c)))
            disp('non finite values');
            break;
        end

        i=i+1; % increase counter
    end

    % -----plot surface crevasses
%     % resistive stress (Pa)
%     Rxx = 2*nthroot(dUdx./(E.*A),n); 
%     % height above buoyancy (m)
%     Hab = H+rho_sw/rho_i*(b); 
%     % surface crevasse penetration depth (m)
%     crev = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*(DFW)); 
%     % basal crevasse penetration depth (m)
%     crev_b = rho_i/(rho_sw-rho_i).*(Rxx./(rho_i*g)-Hab);     
%     % plot
%     figure(10); clf
%     subplot(1,2,1); hold on; 
%     set(gca,'fontsize',12,'linewidth',1);
%     legend('location','southeast');
%     grid on;
%     plot(x/10^3, h,'-b','linewidth',2,'displayname','h');
%     plot(x/10^3, crev,'-m','linewidth',2,'displayname','crev_s');
%     plot(x/10^3, crev_b,'-c','linewidth',2,'displayname','crev_b');
%     xlabel('distance along centerline [km]');
%     ylabel('elevation [m]');
%     subplot(1,2,2); hold on; 
%     set(gca,'fontsize',12,'linewidth',1);
%     legend('location','southwest');
%     grid on;
%     plot(x/10^3, h-crev_b,'-c','linewidth',2,'displayname','h - crev_b');
%     plot(x/10^3, h-crev,'-m','linewidth',2,'displayname','h - crev_s');
%     xlabel('distance along centerline [km]');
%     ylabel('elevation [m]');

    % -----solve for calving criteria (DFW)
    % Solve for DFW that satisfies calving criterion at 2002 calving front 
    % position using current steady-state conditions. 
    %   - Surface crevasse penetration depth [m]: 
    %       crev_s = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*(DFW));
    %   - At xcf, crev_s = h(c). Substitute, rearrange to solve for DFW:
    DFW = (rho_i/rho_fw) * (h(c) - (Rxx(c)/(rho_i*g))); % fresh water depth in crevasses [m]

    % ----------B. CALVING FRONT EVOLUTION----------

    disp('Simulating post-ice shelf collapse conditions');
    
    % -----time stepping [s]
    t_start = 0*3.1536e7; % 2002
    t_end = 98*3.1536e7; % 2022
    t = [t_start:dt:t_end];
    
    % -----initialize variables to track throughout model run
    Fgl = zeros(1,length(t)); % grounding line discharge [Gt/a]
    XCF = zeros(1,length(t)); % calving front position [m]
    XGL = zeros(1,length(t)); % grounding line position [m]
    col = parula(length(t)); % color scheme for plotting

    % -----run flowline model
    for i=1:length(t)

        % decrease DFW every 1 year until ~2019
        if i>1 && t(i)/3.1536e7 < 20 %&& t(i)/3.1536e7 > 2
            if DFW>10
                DFW = DFW-0.0017;
            else
                DFW=10;
            end
        end
        % increase DFW linearly at each time increment to reach delta_DFW by 2100
        if t(i)/3.1536e7 > 20
            delta_DFWi = delta_DFW/(2100-2022)*t(i)/3.1536e7; % total increase in DFW from 2022
        else
            delta_DFWi = 0;
        end
        DFW=DFW+delta_DFWi;
        
        % add backstress after year 5 to account for sea ice occurence
        if t(i)/3.1536e7 > 5; sigma_b = 20e3; end
%         if t(i)/3.1536e7 > 4 && t(i)/3.1536e7 < 7
%             % increase linearly until reaching 50 kPa in year 7
%             sigma_b = sigma_b + 20e3/(find(t/3.1536e7 < 7, 1, 'last') - find(t/3.1536e7 > 4, 1, 'first')); % Pa
%         end
        
        % -----calving front location 
        if i==1
            c=c;
            xcf=x(c);
        else
            % resistive stress [Pa]
            Rxx = 2*nthroot(dUdx./(E.*A),n); 
            % surface crevasse penetration depth [m]
            crev = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*DFW); 
            % interpolate where h-crev <=0 
            % make sure sample points are unique for interpolation
            [~,I] = unique(h-crev); sort(I);
            xcf = interp1(h(I)-crev(I), x(I), 0, 'linear', 'extrap'); % [m along centerline]
            % use a linear fit if calving front is not found
            if isempty(xcf) || xcf<=0
                xcf = interp1(polyval(polyfit(x, h-crev, 1),x), x, 0, 'linear', 'extrap'); % [m along centerline]
                disp('poly');
            end
        end
        XCF(i) = xcf; % save calving front position over time

        % -----grounding line location
        % calculate the thickness required to remain grounded at each grid cell
        Hf = -(rho_sw./rho_i).*b; % flotation thickness (m)
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

        % -----create coordinate system that hits gl exactly
        xl = round(xgl/dx0); % number of ideal grid spaces needed to reach the grounding line
        dx = xgl/xl; % new grid spacing (should be ~dx0)
        xn = 0:dx:xgl; % new distance vector
        if xcf-xgl > 0
            xl = round((xcf-xgl)/dx0);
            dx = (xcf-xgl)/xl;
            xn = [xn xn(end)+dx:dx:xcf];
        end
        clear dx; dxn = [xn(2:end)-xn(1:end-1) xn(end)-xn(end-1)];

        % -----get geometry on new coordinates
        c = length(xn); % index for xcf
        H = interp1(x,H,xn,'linear','extrap');
        % H gradient
        dHdx(1) = (H(2)-H(1))/(xn(2)-xn(1)); % forward difference
        dHdx(2:c-1) = (H(3:c)-H(1:c-2))./(xn(3:c)-xn(1:c-2)); % central difference
        dHdx(c) = (H(c)-H(c-1))./(xn(c)-xn(c-1)); % backward difference            
        % if there is a sudden jump in H (large gradient) past the grounding
        % line, set that to the previous point
        % thickness gradient
        if any(dHdx>50)
            c = find(dHdx>50,1,'first')-1;            
            H = H(1:c);
            xn = xn(1:c);
        end
        gl = dsearchn(xn',xgl); % index for xgl
        b = interp1(x0,b0,xn,'linear','extrap');
        W = interp1(x0,W0,xn,'linear','extrap');
        U = interp1(x,U,xn,'linear','extrap');
        A = interp1(x0,A0,xn,'linear','extrap');
        beta = interp1(x0,beta0,xn,'linear','extrap');
        x = xn; dx = dxn; 
        clear xn dxn; % clear to prevent issues in future iterations
        XGL(i) = xgl; % save grounding line position over time

        % -----calculate surface elevation
        h = b+H; % surface elevation (m a.s.l.)
        h(gl+1:c) = (1-rho_i/rho_sw).*H(gl+1:c); % adjust the surface elevation of ungrounded ice to account for buoyancy
        H(h<0)=0-b(h<0); h(h<0)=0; % surface cannot go below sea level
        h(h-H<b) = b(h-H<b)+H(h-H<b); % thickness cannot go beneath bed elevation

        % -----plot geometry, speed, & grounding line and calving front positions
        col = parula(length(t)); % color scheme for plots
        if plotTimeSteps && mod(t(i),dt*1000)==0 % display every dt*1000
            if i==1
                % plot observed terminus positions
                figure(1); 
                plot(ax3, term.x/10^3, term.date-term.date(1), '.k', 'markersize', 20, 'displayname','observed');
            else
                figure(1);
                % glacier geometry
                plot(ax1,x(1:c)/10^3,h(1:c),'-','color',col(i,:),'linewidth',2,'displayname',num2str(round(t(i)./3.1536e7)+2009));
                plot(ax1,x(gl:c)/10^3,h(gl:c)-H(gl:c),'-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
                plot(ax1,[x(c);x(c)]/10^3,[h(c);h(c)-H(c)],'-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
                title(ax1, [num2str(t(i)/3.1536e7), 'yrs, DFW = ',num2str(DFW), 'm']);
                plot(ax1,x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off'); % calving front
                plot(ax1,x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off'); % floating bed (gl:c)
                % ice speed
                plot(ax2,x(1:c)/10^3,U(1:c).*3.1536e7,'-','Color',col(i,:),'linewidth',2,'DisplayName',num2str(round(t(i)./3.1536e7)+2009));
                % calving front & grounding line positions
                plot(ax3,x(c)/10^3,t(i)/3.1536e7,'.','Color',col(i,:),'markersize',15,'displayname',num2str(round(t(i)./3.1536e7)+2009)); hold on;
                plot(ax3,x(gl)/10^3,t(i)/3.1536e7,'x','Color',col(i,:),'markersize',10,'linewidth',2,'HandleVisibility','off'); hold on;
                title(['sigma_b = ',num2str(round(sigma_b/1000)),' kPa']);
            end
        end

        % -----calculate the effective pressure 
        % (ice overburden pressure minus water pressure) assuming an easy & 
        % open connection between the ocean and ice-bed interface
        sl = find(b<=0,1,'first'); % find where the glacier base first drops below sea level
        N_ground = rho_i*g*H(1:sl); % effective pressure where the bed is above sea level (Pa)
        N_marine = rho_i*g*H(sl+1:length(x))+(rho_sw*g*b(sl+1:length(x))); % effective pressure where the bed is below sea level (Pa)
        N = [N_ground N_marine];
        N(N<0)=0; % cannot have negative values

        % -----solve for new velocity
        [U,dUdx] = U_convergence(x,U,H,h,A,E,N,W,dx,c,n,m,beta,rho_i,rho_sw,g,sigma_b);

        % -----calculate ice flux
        F = U.*H.*W; % ice flux (m^3 s^-1)
        F(isnan(F))=0;
        F(1)=F(2);

        % -----save grounding line ice flux (discharge)
        Fgl(i) = F(gl)*917*1e-12*3.1536e7; % Gt/a

        % -----implement SMB, SMR, & RO
        SMB = interp1(x0,SMB0+Q0,x);
        RO = interp1(x0,RO0,x);
        SMR = zeros(1,c);
        if t(i)/3.1536e7 < 20
            % SMR change
            TF = TF0;
            delta_mdot = ((C*-b(gl)*nthroot((sum(RO0(1:gl0))*86400),1/0.39) + 0.15)*(TF^1.18))/86400-mdot0; % m/s
            if gl<c
                SMR(gl+1:c) = (SMR0-delta_mdot)/(SMR_mean_fit.a+1)*feval(SMR_mean_fit,x(gl+1:c)-x(gl));
            end
        else
            % SMB change
            delta_SMBi = delta_SMB/(2100-2022)*(t(i)/3.1536e7-10); % total increase in SMB from 2022 rate 
            for k=1:c
                SMB(k) = SMB(k)+delta_SMBi*(h0(1)-h(k))/(h0(1)-h0(c0)); 
            end
            % TF change
            TFi = delta_TF/(2100-2022)*(t(i)/3.1536e7-10); % total increase in TF from 2022
            TF = TFi + TF0;
            % increase runoff due to surface melting if running the 'enhanced SMB' scenario        
            if SMB_enhance
                RO = (interp1(x0,SMB0,x)-SMB)+interp1(x0,RO0,x);   
            else
                RO = interp1(x0,RO0,x);
            end
            % SMR change
            delta_mdot = ((C*-b(gl)*nthroot((sum(RO0(1:gl0))*86400),1/0.39) + 0.15)*(TF^1.18))/86400-mdot0; % m/s
            if gl<c
                SMR(gl+1:c) = (SMR0-delta_mdot)/(SMR_mean_fit.a+1)*feval(SMR_mean_fit,x(gl+1:c)-x(gl));
            end
        end 
        if plotClimateParams && mod(t(i),dt*1000)==0 % display every dt*1000
            if i==1
                figure(2); clf
                drawnow
                set(gcf,'Position',[50 100 900 600]);
                % SMB
                axA = subplot(2, 2, 1);
                hold on; grid on;
                set(gca,'FontSize',12,'linewidth',2);
                title('SMB');
                ylabel('[m/y]');
                plot(axA, x./10^3, SMB, 'color', col(i,:), 'linewidth', 2, 'displayname', '2022');
                % TF
                axB = subplot(2, 2, 2); 
                hold on; grid on;
                set(gca,'FontSize',12,'linewidth',2);
                title('TF');
                ylabel('[^oC]');
                plot(axB, t(i)./3.1536e7, TF, '.', 'Color', col(i,:), 'markersize', 20, 'displayname', '2022');
                % SMR
                axC = subplot(2, 2, 3);
                hold on; grid on;
                set(gca,'FontSize',12,'linewidth',2);
                title('SMR');
                xlabel('distance along centerline [km]'); ylabel('[m/y]');
                plot(axC, x/10^3, SMR, 'color', col(i,:), 'linewidth', 2, 'displayname', '2022');
                % DFW
                axD = subplot(2, 2, 4); 
                hold on; grid on;
                set(gca,'FontSize',12,'linewidth',2);
                title('DFW');
                xlabel('Year'); ylabel('[m]');
                plot(axD, t(i)./3.1536e7, DFW, '.', 'Color', col(i,:), 'markersize', 20, 'displayname', '2022');
            else
                plot(axA, x./10^3, SMB, 'color', col(i,:), 'linewidth', 2);
                plot(axB, t(i)/3.1536e7, TF, '.', 'Color', col(i,:), 'markersize', 20); 
                plot(axC, x/10^3, SMR, 'color', col(i,:), 'linewidth', 2);
                plot(axD, t(i)./3.1536e7, DFW, '.', 'Color', col(i,:), 'markersize', 20);
            end
        end

        % -----calculate the  change in ice thickness from continuity
        clearvars dHdt
        dHdt(1) = (-1/W(1))*(F(1)-F(2))/(x(1)-x(2)); % forward difference
        dHdt(2:c-1) = (-1./W(2:c-1)).*(F(1:c-2)-F(3:c))./(x(1:c-2)-x(3:c)); % central difference
        dHdt(c:length(x)) = (-1./W(c:length(x))).*(F(c-1:length(x)-1)-F(c:length(x)))./(x(c-1:length(x)-1)-x(c:length(x))); % backward difference
        dH = dHdt.*dt;

        % -----new thickness (change from dynamics, SMB, & SMR)
        Hn = H+dH+(SMB.*dt)+(SMR.*dt)-(interp1(x0,RO0,x)*dt);
        Hn(Hn < 0) = 0; % remove negative values
        H = Hn; % set thickness to new value

        % -----stop the model if it behaves unstably (monitored by ice thickness and speed)
        if max(H) > H_max
            disp(['Adjust dt']);
            break;
        end
        if mean(U) < U_min
            disp('Too slow!');
            break;
        end
        if any(~isfinite(H(1:c))) || any(~isfinite(U(1:c))) || any(~isfinite(h(1:c)))
            disp('non finite values');
            break;
        end

    end

    % save geometry & speed
    if saveFinal
        % store final geometry & speed
        h2=h; H2=H; x2=x; c2=c; gl2=gl; U2=U; Fgl2=Fgl; XCF2=XCF; XGL2=XGL; DFW2 = DFW0+delta_DFW;
        % no change
        if delta_DFW==0 && delta_TF==0 && delta_SMB==0 
            save([homepath,'workflows/steady-state-initial/results/1_SMB_DFW_TF/SMB0_DFW0m_TF0_geom.mat'],'h2','H2','c2','U2','gl2','x2','DFW2','Fgl2','XGL2','XCF2');
            h1=h; H1=H; b1=b; c1=c; U1=U; x1=x; gl1=dsearchn(x1',XGL(end)); DFW1=DFW0; Fgl1=Fgl; XGL1=XGL; XCF1=XCF; 
            save([homepath,'workflows/steady-state-initial/results/1_SMB_DFW_TF/2100_noChange_steady_state_initial.mat'],'h1','H1','b1','c1','U1','gl1','x1','DFW1','Fgl1','XGL1','XCF1');
            disp('2100_noChange_steady_state_initial saved');
        % 2) SMB_enhanced
        elseif SMB_enhance==1 && delta_TF==0
            fileName = ['SMB',num2str(delta_SMB*3.1536e7),'_enh_geom.mat'];
            save([homepath,'workflows/steady-state-initial/results/2_SMB_enh/',fileName],'h2','H2','c2','U2','gl2','x2','DFW2','Fgl2','XGL2','XCF2');
            disp('geometry saved (2)');
        % 3) SMB_enhanced + TF
        elseif SMB_enhance==1 && delta_TF~=0 
            fileName = ['SMB',num2str(delta_SMB*3.1536e7),'_enh_dTF',num2str(delta_TF),'_geom.mat'];
            save([homepath,'workflows/steady-state-initial/results/3_SMB_enh+TF/',fileName],'h2','H2','c2','U2','gl2','x2','DFW2','Fgl2','XGL2','XCF2');
            disp('geometry saved (3)');           
        % 1) SMB, DFW, and TF independent
        else
            fileName = ['SMB',num2str(delta_SMB*3.1536e7),'_DFW',num2str(DFW2),'m_TF',num2str(delta_TF),'_geom.mat'];
            save([homepath,'workflows/steady-state-initial/results/1_SMB_DFW_TF/',fileName],'h2','H2','c2','U2','gl2','x2','DFW2','Fgl2','XGL2','XCF2');
            disp('geometry saved (1)');
        end
    else
        disp('geometry not saved.');
    end

    % plot results
    i=1;
    while i==1
        load([homepath,'inputs-outputs/2100_noChange_steady_state_initial.mat']);
        h2=h; H2=H; x2=x; c2=c; gl2=gl; U2=U; DFW2=DFW0+delta_DFW; Fgl2=Fgl; XCF2=XCF; XGL2=XGL; % store final geometry & speed
        gl1=dsearchn(x2',XGL(end)); 
        figure(11); clf % sensitivity test changes
        hold on; grid on;
        set(gcf,'Position',[250 80 886 686]);
        set(gca,'FontSize',14,'linewidth',2,'fontweight','bold');
        xlabel('distance Along Centerline (km)'); ylabel('elevation (m)');
        legend('Location','east'); xlim([0 100]); ylim([min(b0)-100 max(h)+100]);
        title(['TF = + ',num2str(round(delta_TF,1)),'^oC, SMB = + ',...
            num2str(round(delta_SMB*3.1536e7,1)),'m/a, DFW = ',num2str(DFW2),'m']);
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
            plot(x1/10^3,b1,'k','linewidth',2,'HandleVisibility','off');
        % inset plot of terminus
        ax2 = axes('Position',[0.62 0.62 0.28 0.28]);
            hold on; grid on; set(gca,'linewidth',2,'fontweight','bold','fontsize',9);
            delta_L = x2(c2)-x1(c1); % change in length (m)
            delta_H = mean(H2(1:c2))-mean(H1(1:c1)); % change in ice thickness (m)
            delta_U = mean(U2(1:c2))-mean(U1(1:c1)); % change in ice speed (m/s) 
            title(['\Delta L=',num2str(round(delta_L)),'m, ','\DeltaH_{\mu}=',...
                num2str(round(delta_H)),' m, ','\DeltaU_{\mu}=',...
                num2str(round(delta_U.*3.1536e7)),' m/a']);
            xlim([x(gl)/10^3-5 x(c)/10^3+5]); ylim([min(b0)-100 300]);
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
            plot(x/10^3,b,'k','linewidth',2,'HandleVisibility','off');
            % mean sea level
            plot([x(1),x(end)]/10^3,[0,0],'k--','HandleVisibility','off'); drawnow
            i=i+1;
    end

    disp('Simulation complete');
end
    
    
    