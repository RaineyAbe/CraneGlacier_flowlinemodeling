function [x,U,h,hb,H,gl,c,xcf,dUdx,Fgl,XCF,XGL] = flowlineModel(homepath,plotTimeSteps,plotMisfits,plotClimateParams,dt,t_start,t_end,beta0,DFW0,delta_SMB,delta_DFW,delta_TF)
% Rainey Aberle, 2021
% Adapted from code authored by Enderlin et al. (2013), doi:10.5194/tc-7-1007-2013
% Function to run the flowline model using variables saved in the
% flowlineModelInitialization.m file. 
%
% INPUTS:   homepath = path to "CraneGlacier_flowlinemodeling" in directory 
%                      (string ending with "/")
%           plotTimeSteps = option to plot the model geometry and speed at
%                           several model time steps (0 = no, 1 = yes).
%           plotMisfits = option to plot the misfits between observed and
%                         modeled surface elevation and speed in model 
%                         year 2018 (0 = no, 1 = yes).
%           dt = model time step (s)
%           t_start = model beginning time (s), often set to 0. 
%           t_end = model end time (s). 
%
% OUTPUTS:  x = resulting model grid vector (m)
%           U = modeled centerline speed (m/s)
%           h = modeled centerline surface elevation (m.a.s.l.)
%           hb = modeled centerline bed elevation (m.a.s.l.)
%           H = modeled centerline thickness (m.a.s.l.)
%           gl = modeled grounding line location (index of model grid)
%           c = modeled calving front location (index of model grid)
%           xcf = modeled calving front location (m along centerline)
%           dUdx = modeled centerline strain rates (1/s)
%           beta0
%           d_fw0
%
% NOTE(s):  Requires U_convergence.m 

%----------------------------------
%----------INITIALIZATION----------
%----------------------------------

% Add necessary paths
addpath([homepath,'matlabFunctions/cmocean_v2.0/cmocean/']);
addpath([homepath,'inputs-outputs/']);
addpath([homepath,'scripts/']);

dx0 = 200; % grid spacing (m)

% Load Crane Glacier initialization variables
b=1; % use a while loop to make this section collapsible
while b==1
    
    % load initial conditions
    load('flowlineModelInitialization.mat','x0','h0','hb0','W0','c0',...
        'A0','Q0','SMB0','SMR0','RO0','U0');
    gl0=c0;
    SMR_mean_fit = load('LarsenC_MeanMeltRate.mat').mr_mean_fit;
    % set default values for beta0 and d_fw0 if not provided
%     if isempty(varargin{7})
%         beta0 = load('flowlineModelInitialization.mat','beta0');
%     end
%     if isempty(varargin{8})
%         d_fw0 = load('flowlineModelInitialization.mat','FWD0');
%     end

    % Centerline
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

    % densities and g
    rho_i = 917; % ice density (kg m^-3)
    rho_sw = 1028; % ocean water density (kg m^-3)
    rho_fw = 1000; % fresh water density (kg m^-3)
    g = 9.81; % acceleration (m s^-2)

    % stress parameters (unitless)
    m = 3; % basal sliding exponent
    n = 3; % flow law exponent
    E = 1; % enhancement factor

    % ocean thermal forcing
    TF0 = 0.2; % ^oC - estimated from Larsen B icebergs
    % estimate initial melt rate using Eqn from Slater et al. (2020):
    mdot0 = (3*10^-4*-hb0(gl0)*((sum(RO0(1:gl0)))*86400)^0.39 + 0.15)*TF0^1.18/86400; % m/s
    
    % calving parameters
    Hc = 100; % m -> set the calving front to a default minimum ice thickness value if calving criteria is not met
    sigma_b = 1000; % back pressure (Pa) - similar to that employed at Helheim Glacier (Nick et al., 2009)

    % maximum & minimum thickness & speed cut-off to check for instability
    H_max = 2000; % maximum thickness (m)
    H_min = 100;  % minimum thickness (m)
    U_min = 100./3.1536e7;  % minimum mean speed (m s^-1)

    % initial conditions
    H0 = h0-hb0; % ice thickness (m)
    dUdx0 = [(U0(2:end)-U0(1:end-1))./(x0(2:end)-x0(1:end-1)) 0]; % strain rate (1/s) %%EE: flipped this to be consistent w/ other computations 15/09/21
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
    t = (t_start:dt:t_end);
    
    % initialize variables
    x=x0; U=U0; W=W0; gl=gl0; dUdx=dUdx0; A=A0; h=h0; hb=hb0; H=H0; 
    DFW=DFW0; beta=beta0; dx=dx0;
    Fgl=zeros(1,length(t)); % ice mass flux across the grounding line
    XCF = NaN*ones(1,length(t)); XGL = NaN*ones(1,length(t)); % store xcf and xgl over time

    b=b+1; % exit loop
end
    
%--------------------------------------
%----------RUN FLOWLINE MODEL----------
%--------------------------------------

for i=1:length(t)

    % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
    if i==1 % use observed calving front position for first iteration
        c=c0;
        xcf = x0(c0);
    else
        % to allow for advance, extrapolate U, dUdx, h, hb, and H past the calving front
        U_cf = interp1(x(1:c),U(1:c),x0,'linear','extrap');
        dUdx_cf(1) = (U_cf(2)-U_cf(1))/(x0(2)-x0(1)); % forward difference at upper boundary
        dUdx_cf(2:length(x0)-1) = (U_cf(3:length(x0))-U_cf(1:length(x0)-2))...
            /(x0(3:length(x0))-x0(1:length(x0)-2)); % central difference
        dUdx_cf(length(x0)) = (U_cf(length(x0))-U_cf(length(x0)-1))/(x0(end)-x0(end-1)); % backward difference at end
        h_cf = interp1(x(1:c),h(1:c),x0,'linear','extrap');
            h_cf(c-50:end) = feval(fit(x0(c-50:c)',h(c-50:c)','poly1'),x0(c-50:end))';
        hb_cf = hb0;    
        H_cf = interp1(x(1:c),H(1:c),x0,'linear','extrap'); 
        Rxx = 2*nthroot(dUdx_cf./(E.*A0),n); % resistive stress (Pa)
        % increase DFW linearly at each time increment to reach delta_DFW by 2100
        if t(i) >= 10*3.1536e7
            delta_DFWi = delta_DFW/(2100-2019)*(t(i)/3.1536e7-10); % total increase in smr from 2019 rate
        else
            delta_DFWi = 0;
        end
        DFW=DFW0+delta_DFWi;
        crev_s = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*(DFW)); % surface crevasse penetration depth (m)
        Hab = H_cf+rho_sw/rho_i*(hb_cf); % height above buoyancy (m)
        crev_b = rho_i/(rho_sw-rho_i).*(Rxx./(rho_i*g)-Hab); % basal crevasse depth (m)    % calving front located where the inland-most crevasse intersects sea level        
        xcf_s = interp1(h_cf(find(h_cf-crev_s<0,1,'first')-1:find(h_cf-crev_s<0,1,'first')+1)...
            -crev_s(find(h_cf-crev_s<0,1,'first')-1:find(h_cf-crev_s<0,1,'first')+1),...
            x0(find(h_cf-crev_s<0,1,'first')-1:find(h_cf-crev_s<0,1,'first')+1),0,'linear','extrap'); % (m along centerline)        
%         xcf_b = interp1(h_cf(find(h_cf-crev_b<0,1,'first')-1:find(h_cf-crev_b<0,1,'first')+1)...
%             -crev_b(find(h_cf-crev_b<0,1,'first')-1:find(h_cf-crev_b<0,1,'first')+1),...
%             x0(find(h_cf-crev_b<0,1,'first')-1:find(h_cf-crev_b<0,1,'first')+1),0,'linear','extrap'); % (m along centerline)
%         % calving front located where the inland-most crevasse intersects sea level
%         if xcf_s<xcf_b
%             xcf = xcf_s;
%         else
%             xcf = xcf_b;
%         end
        xcf=xcf_s;
    end
    
    % calculate the thickness required to remain grounded at each grid cell
    Hf = -(rho_sw./rho_i).*hb; % flotation thickness (m)
    % find the location of the grounding line and use a floating
    % geometry from the grounding linU_widthavge to the calving front
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

    % create coordinate system that hits cf and gl exactly
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

    % get geometry on new coordinates
    c = length(xn);
    H = interp1(x,H,xn,'linear','extrap');
    % if there is a sudden jump in H (large gradient) past the grounding
    % line, set that to the previous point
        % thickness gradient
        dHdx = gradient(H);
        if any(dHdx>50)
            c = find(dHdx>50,1,'first')-1;            
            H = H(1:c);
            xn = xn(1:c);
        end
    gl = dsearchn(xn',xgl); % indices for xcf and xgl
    hb = interp1(x0,hb0,xn,'linear','extrap');
    W = interp1(x0,W0,xn,'linear','extrap');
    U = interp1(x,U,xn,'linear','extrap');
    A = interp1(x0,A0,xn,'linear','extrap');
    beta = interp1(x0,beta0,xn,'linear','extrap');
    x = xn; dx = dxn; clear xn dxn; %EE: added the clear statement 14/09/21
    XCF(i) = xcf; % save calving front position over time
    XGL(i) = xgl; % save grounding line position over time

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
    [U,dUdx] = U_convergence(x,U,U0,H,h,A,E,N,W,dx,c,n,m,beta,rho_i,rho_sw,g,sigma_b,i);
    %EE: removed U0 and dUdx from inputs to U_convergence 15/09/21
    
    % calculate ice flux
    F = U.*H.*W; % ice flux (m^3 s^-1)
    F(isnan(F))=0;
    F(1)=F(2);

    % implement SMB, SMR, RO, delta_SMB, & delta_SMR
    if t(i)/3.1536e7<10 % use original SMB & SMR for first 10 model years  
        SMB = interp1(x0,SMB0+Q0,x);
        RO = interp1(x0,RO0,x);
        delta_mdot = 0; % m/s
        SMR = zeros(1,c);
        % use the Larsen C mean melt rate profile to scale SMR
        % using the max initial SMR
        if gl<c
            SMR(gl+1:c) = (SMR0+delta_mdot)/(SMR_mean_fit.a+1)*feval(SMR_mean_fit,x(gl+1:c)-x(gl+1)); 
        end
    elseif t(i)/3.1536e7>=10 % implement changes after 10 model years
        RO = (interp1(x0,SMB0,x)-SMB)+interp1(x0,RO0,x);        
        % calculate additional melt due to the increase in subglacial discharge
        TFi = delta_TF/(2100-2019)*(t(i)/3.1536e7-10); % total increase in TF from 2019 
        delta_mdot = ((3*10^-4*-hb(gl)*((sum(RO(1:gl))*86400)^0.39) + 0.15)*((TFi+TF0)^1.18))/86400-mdot0; % m/s
        SMR = zeros(1,c);
        if gl<c
            SMR(gl+1:c) = (SMR0-delta_mdot)/(SMR_mean_fit.a+1)*feval(SMR_mean_fit,x(gl+1:c)-x(gl));
        end
        delta_SMBi = delta_SMB/(2100-2019)*(t(i)/3.1536e7-10); % total increase in smb from 2019 rate 
        SMB = interp1(x0,SMB0+Q0,x);
        for k=1:c
            SMB(k) = SMB(k)+delta_SMBi*(h0(1)-h(k))/(h0(1)-h0(c0)); 
        end
    end
    if plotClimateParams
        if i==1
            figure(10); clf; set(gcf,'position',[200 500 1000 500]);
            subplot(1,3,1); hold on; grid on; % SMB
            set(gca,'fontsize',14,'fontname','Arial','linewidth',2);
            ylabel('SMB (m a^{-1})');
            plot(x/10^3,SMB/3.1536e7,'linewidth',2,'color',col(i,:));
            subplot(1,3,2); hold on; grid on; % DFW
            set(gca,'fontsize',14,'fontname','Arial','linewidth',2);
            xlabel('Year'); ylabel('DFW (m)');
            plot(t(i)/3.1536e7+2009,DFW,'.','markersize',15,'color',col(i,:));
            subplot(1,3,3); hold on; grid on; % SMR
            set(gca,'fontsize',14,'fontname','Arial','linewidth',2);
            ylabel('SMR (m a^{-1})');
            plot(x/10^3,SMR/3.1536e7,'linewidth',2,'color',col(i,:));
        elseif mod(i-1,round(length(t)/10))==0 % display every length(t)/10
            figure(10);
            subplot(1,3,1); % SMB
            plot(x/10^3,SMB/3.1536e7,'linewidth',2,'color',col(i,:));
            subplot(1,3,2); % DFW
            plot(t(i)/3.1536e7+2009,DFW,'.','markersize',15,'color',col(i,:));
            subplot(1,3,3); % SMR
            plot(x/10^3,SMR/3.1536e7,'linewidth',2,'color',col(i,:));
        end
    end
    
    % calculate the  change in ice thickness from continuity
    clearvars dHdt
    dHdt(1) = (-1/W(1))*(F(1)-F(2))/(x(1)-x(2)); % forward difference
    dHdt(2:c-1) = (-1./W(2:c-1)).*(F(1:c-2)-F(3:c))./(x(1:c-2)-x(3:c)); % central difference
    dHdt(c:length(x)) = (-1./W(c:length(x))).*(F(c-1:length(x)-1)-F(c:length(x)))./(x(c-1:length(x)-1)-x(c:length(x))); % backward difference
    dH = dHdt.*dt;

    % new thickness (change from dynamics, SMB, & SMR)
    Hn = H+dH+(SMB.*dt)+(SMR.*dt)-(interp1(x0,RO0,x)*dt);
    Hn(Hn < 0) = 0; % remove negative values
    H = Hn; % set as the new thickness value

    Fgl(i) = F(gl)*pi*1/4*917*1e-12*3.1536e7; % Gt/a

    % stop the model if it behaves unstably (monitored by ice thickness and speed)
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

% calculate and plot misfits
if plotMisfits
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
        yyaxis right; ylabel('Basal Roughness Factor (s^{1/m} m^{-1/m})'); ylim([min(beta0)-0.2 max(beta)+0.2]);
            plot(x/10^3,movmean(beta,5),'-','linewidth',2,'color',[0 0.4 0.8],'displayname','\beta');
    
    % display grounding line misfit
    disp('Mean misfits at the grounding line:');
    disp(['    Speed: ',num2str((U(gl)-interp1(cl.x,U_obs(10).U,x(gl)))*3.1536e7),' m/yr']);
    disp(['    Elevation: ',num2str((h(gl)-interp1(cl.x,h_obs(36).surface,x(gl)))),' m']);

end

% Plot stresses
% figure(3); clf; 
% set(gcf,'position',[324    83   913   369]);
% col2 = parula(4); 
% hold on; grid on; legend; 
% set(gca,'fontsize',14,'linewidth',2); xlabel('x (km)'); 
% yyaxis left; ylabel('Td (Pa)');plot(x/10^3,Td,'linewidth',2,'color',col2(1,:),'displayname','T_d');
% yyaxis right; ylabel('Pa'); ylim([-1e11 1e11]);
% plot(x/10^3,Tlon,'linewidth',2,'color',col2(2,:),'displayname','T_{lon}');
% plot(x/10^3,Tlatb,'linewidth',2,'color',col2(3,:),'displayname','T_{lat,b}');
% 
% % Plot viscosity, A, U and dUdx
% figure(4); clf;
% set(gcf,'position',[322   331   920   366]);
% subplot(1,2,1);
% hold on; grid on; legend;
% set(gca,'fontsize',14,'linewidth',2); xlabel('x (km)'); 
% yyaxis left; ylabel('vm (Pa s)'); plot(x/10^3,vm,'-m','linewidth',2,'displayname','vm');
% yyaxis right; ylabel('A (Pa^{-3} a^{-1})'); plot(x/10^3,A*3.1536e7,'-b','linewidth',2,'displayname','A');
% plot(x/10^3,5.6e-17*ones(1,length(x)),'--b','displayname','A (Nick et al., 2010)','linewidth',2);
% subplot(1,2,2);
% hold on; grid on; legend;
% set(gca,'fontsize',14,'linewidth',2); xlabel('x (km)'); 
% yyaxis left; ylabel('m/y'); plot(x/10^3,U*3.1536e7,'-k','linewidth',2,'displayname','U');
% yyaxis right; ylabel('1/y'); plot(x/10^3,dUdx*3.1536e7,'--k','linewidth',2,'displayname','dUdx');
% 
% % plot SMB, SMR, H
% figure(5); clf; 
% set(gca,'fontsize',14,'linewidth',2);
% hold on; grid on; legend; xlabel('x (km)');
% yyaxis left; ylabel('(m/yr)');
% plot(x/10^3,smb.*3.1536e7,'-k','displayname','SMB','linewidth',2);
% plot(x/10^3,smr.*3.1536e7,'--k','displayname','SMR','linewidth',2);
% yyaxis right; ylabel('H (m)');
% plot(x/10^3,H,'-b','displayname','H','linewidth',2);
% plot(x/10^3,dH,'--b','displayname','dH','linewidth',2);
