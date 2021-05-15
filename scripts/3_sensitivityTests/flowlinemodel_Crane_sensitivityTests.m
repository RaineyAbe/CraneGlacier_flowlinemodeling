%% Glacier Flowline Model: Sensitivity Tests for Crane Glacier, Antarctic Peninsula
% Last edited: Spring 2021
% Adapted by Rainey Aberle from Ellyn Enderlin's flowline model package 
% (Enderlin et al., 2013)
% 
%   0. Define time and space independent variables by loading the flowline 
%       initialization file and regridding to the desired grid spacing. 
%   1. Run sensitivity tests for surface mass balance (SMB), submarine
%       melting rate (SMR), and fresh water depth in crevasses (fwd). 

%% 0. define time and space independent variables

clear all; close all;
warning off; % turn off warnings (velocity coefficient matrix is close to singular)
    
% define home path in directory and add necessary paths
homepath = '/Users/raineyaberle/Desktop/Research/CraneGlacier_flowlinemodeling/';
addpath([homepath,'scripts/sensitivityTests/']); % add path to U_convergence
addpath([homepath,'../matlabFunctions/cmocean_v2.0/cmocean']);
cd([homepath,'inputs-outputs/']);

dx0 = 200; % grid spacing (m)

% Load Crane Glacier initialization variables
load('Crane_flowlineModelInitialization.mat');
F0=0; % inner boundary flux (m^3/s)

% Load observed conditions
% dH
dH_obs = load('Crane_dHdt_2009-2018.mat').dHdt.dH_total; % (m) total change in thickness 2009-2018
h_obs = load('Crane_surfaceElevationObs.mat').h;
h_obs_2019 = load('Crane_surfaceElevationObs.mat').h(36).surface;
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
H0(gl0:length(x0))=h0(gl0:end)*rho_sw/(rho_sw-rho_i); % buoyant thickness using surface
H0(H0>=(h0-hb0))=h0(H0>=(h0-hb0))-hb0(H0>=(h0-hb0)); % thickness can't go beneath bed elevation
H0(c0+1:end) = 0; % zero ice thickness past calving front

% extend variables to model length if necessary (set extended points to last value)
A0(find(isnan(A0),1,'first'):end) = A0(find(isnan(A0),1,'first')-1);
U0(find(isnan(U0),1,'first'):end) = U0(find(isnan(U0),1,'first')-1);
dUdx0(find(isnan(dUdx0),1,'first'):end) = dUdx0(find(isnan(dUdx0),1,'first')-1);
beta0(find(isnan(beta0),1,'first'):end) = beta0(find(isnan(beta0),1,'first')-1);
hb0(find(isnan(hb0),1,'first'):end) = hb0(find(isnan(hb0),1,'first')-1);
W0(find(isnan(W0),1,'first'):end) = W0(find(isnan(W0),1,'first')-1);    

%% 1. conduct sensitivity tests for SMB, SMR, & fwd

% (1) Run through the the designated number of model years with no change 
%   in SMB or SMR and save the output. 
% (2) Run through simulated 2009-2019, then continue evolving the model  
%   under new scenarios: \Delta(SMB) &/or \Delta(SMR) and compare changes 
%   in geometry and speed with the no change scenario.  

close all; 

save_figure = 0;    % = 1 to save resulting figure
save_final = 1;     % = 1 to save final geometry and speed 
plot_timeSteps = 0; % = 1 to plot geometry, speed, cf/gl positions every decade
timeseries_save = 0; % = 1 to save figures for time series

% load no change conditions
cd([homepath,'inputs-outputs/']);
load('Crane_2100_noChange.mat'); % load no change variables
    
% set up changes in SMB, SMR, & fwd
% note: decrease SMB & SMR in increments of 0.5 m a-1 (1.585e-8) starting
%   from 1 m a-1 until reaching max SMR found at other Antarctic ice shelves: 
%   ~6 m a^-1 (Adusumilli et al., 2020)
delta_smr0 = (0:-1:-10)./3.1536e7; % m/s change in maxmimum SMR
delta_smb0 = (0:-1:-10)./3.1536e7; % m/s change in SMB at the calving front (used to increase gradient)
delta_fwd0 = 0:0.5:5; % m change in fwd

% define time stepping (s)
dt = 0.01*3.1536e7;
t_start = 0*3.1536e7;
t_end = 91*3.1536e7;
t = (t_start:dt:t_end);

% initialize ice flux
Fgl=zeros(1,length(t)); % ice mass flux across the grounding line

for j=1:length(delta_smb0)

    % initialize variables
    x=x0; H=H0; U=U0; W=W0; gl=gl0; dUdx=dUdx0; A=A0; h=h0; hb=hb0; fwd=fwd0; beta=beta0;
    XCF = NaN*ones(1,length(t)); XGL = NaN*ones(1,length(t)); % store xcf and xgl over time
    
    delta_smr = 0;%delta_smr0(j);
    delta_smb = 0;%delta_smb0(j);
    delta_fwd = 0;%delta_fwd0(j);
    
    %try
        % run flowline model
        for i=1:length(t)

            % implement change to fwd after 10 model years
            if t(i)>10*3.1536e7
                % linearly increase fwd to reach delta_fwd by 2100
                delta_fwdi = delta_fwd/(2100-2019)*(t(i)/3.1536e7-10); % total increase in fwd at this time increment
                fwd = fwd0 + delta_fwdi; % m
            end 
            
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
            beta = interp1(x0,beta0,xn,'linear','extrap'); beta(gl+1:end)=0;
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
                    xlim([0 65]); ylim([0 2500]);
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
            [U,dUdx,~,~,~,~,~] = U_convergence(x,U,dUdx,H,h,A,E,N,W,dx,c,n,m,beta,rho_i,rho_sw,g,sigma_b,i);

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
                smb = interp1(x0,smb0+Q0,x);
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
%                         title('FWD');
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

        % plot results
        if t(i)==t_end
            h2=h; H2=H; x2=x; c2=c; gl2=gl; U2=U; fwd2=fwd; Fgl2=Fgl; XCF2=XCF; XGL2=XGL; % store final geometry & speed
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
                plot([x(1),x(end)]/10^3,[0,0],'k--','HandleVisibility','off'); drawnow
        end

        % save figure
        if save_figure && ishandle(10)
            cd([homepath,'scripts/sensitivityTests/results/']);
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
            cd([homepath,'scripts/sensitivityTests/results/']);
            if delta_fwd==0 && delta_smr==0 && delta_smb==0 
                save(['SMR0_SMB0_geom.mat'],'h2','H2','c2','U2','gl2','x2','fwd2','Fgl2','XGL2','XCF2');
            elseif delta_fwd==0 
                fileName = ['SMR',num2str(delta_smr),'_SMB',num2str(delta_smb),'_geom.mat'];
                save(fileName,'h2','H2','c2','U2','gl2','x2','fwd2','Fgl2','XGL2','XCF2');
            else
                fileName = ['fwd_',num2str(fwd),'m_geom.mat'];
                save(fileName,'h2','H2','c2','U2','gl2','x2','fwd2','Fgl2','XGL2','XCF2');
            end
            disp('geometry saved.');
        else
            disp('geometry not saved.');
        end
        
    %catch
    %    disp(['iteration ',num2str(j),' failed']);
    %end
    
end


