%% Glacier Flowline Model: Sensitivity Tests for Crane Glacier, Antarctic Peninsula
% Rainey Aberle 
% Fall 2021
% 
%   0. Define time and space independent variables by loading the flowline 
%       initialization file and regridding to the desired grid spacing. 
%   1. Run sensitivity tests for surface mass balance (SMB), submarine
%       melting rate (SMR), and fresh water depth in crevasses (FWD) independently. 

%% 0. define necessary paths in directory

clear all; close all;
warning off; % turn off warnings (velocity coefficient matrix is close to singular)
    
% define home path in directory and add necessary paths
homepath = '/Users/raineyaberle/Desktop/Research/CraneModeling/CraneGlacier_flowlinemodeling/';
addpath([homepath,'scripts/']); % add path to U_convergence
addpath([homepath,'../matlabFunctions/cmocean_v2.0/cmocean']);
addpath([homepath,'inputs-outputs/']);

%% 1. conduct sensitivity tests for SMB & DFW independently

% (1) Run through the the designated number of model years with no change 
% (2) Run through simulated 2009-2019, then continue evolving the model  
%   under new scenarios: Delta_SMB & Delta_FWD and compare changes 
%   in geometry and speed with the no change scenario.  

close all; 

saveFinal = 1;         % = 1 to save final conditions in homepath/3_sensitivityTests/results/
plotTimeSteps = 1;     % = 1 to plot geometry, speed, cf/gl positions every decade
plotMisfits = 0;       % = 1 to plot misfit with 2018 conditions
plotClimateParams = 1; % = 1 to plot climate parameters

% load no change conditions
load('2100_noChange.mat'); % load no change variables
    
% set up changes in SMB & DFW
% - decrease maximum SMB by increments of 0.5 m a-1 down to -10 m a-1
% - increase DFW by increments of 1 m up to 10 m
delta_SMB0 = (0:-1:-10)./3.1536e7; % m/s change in SMB at the calving front (used to increase gradient)
delta_DFW0 = 0:10; % m change in DFW 

% define time stepping (s)
dt = 0.01*3.1536e7;
t_start = 0*3.1536e7;
t_end = 91*3.1536e7;

% loop through scenarios
for j=1%:length(delta_SMB0)

    load('flowlineModelInitialization.mat');
    delta_SMB = 0;%delta_SMB0(j);
    delta_DFW = delta_DFW0(j);
    delta_TF = 0;

    %try
        % run flowline model
        [x,U,h,hb,H,gl,c,xcf,dUdx,Fgl,XCF,XGL] = flowlineModel(homepath,plotTimeSteps,plotMisfits,plotClimateParams,dt,t_start,t_end,beta0,DFW0,delta_SMB,delta_DFW,delta_TF);
        
        % save geometry & speed
        if saveFinal
            cd([homepath,'scripts/3_sensitivityTests/results/1_SMB_DFW/']);
            if delta_DFW==0 && delta_SMR==0 && delta_SMB==0 
                h2=h; H2=H; x2=x; c2=c; gl2=gl; U2=U; DFW2=DFW0; Fgl2=Fgl; XCF2=XCF; XGL2=XGL; % store final geometry & speed
                save('SMB0_DFW0_geom.mat','h2','H2','c2','U2','gl2','x2','DFW2','Fgl2','XGL2','XCF2');
                cd([homepath,'inputs-outputs/']);
                h1=h; H1=H; hb1=hb; c1=c; U1=U; x1=x; gl1=dsearchn(x1',XGL(end)); DFW1=DFW0; Fgl1=Fgl; XGL1=XGL; XCF1=XCF; 
                save('2100_noChange.mat','h1','H1','hb1','c1','U1','gl1','x1','DFW1','Fgl1','XGL1','XCF1');
                disp('2100_noChange saved');
            else
                fileName = ['SMB',num2str(delta_SMB*3.1536e7),'_DFW',num2str(delta_DFW),'m_geom.mat'];
                DFW2 = DFW0+delta_DFW;
                h2=h; H2=H; x2=x; c2=c; gl2=gl; U2=U; DFW2=DFW0; Fgl2=Fgl; XCF2=XCF; XGL2=XGL; % store final geometry & speed
                save(fileName,'h2','H2','c2','U2','gl2','x2','DFW2','Fgl2','XGL2','XCF2');
            end
            disp('geometry saved.');
        else
            disp('geometry not saved.');
        end

        % plot results
        b=1;
        while b==1
            cd([homepath,'inputs-outputs/']);
            load('2100_noChange.mat');
            h2=h; H2=H; x2=x; c2=c; gl2=gl; U2=U; DFW2=DFW0; Fgl2=Fgl; XCF2=XCF; XGL2=XGL; % store final geometry & speed
            gl1=dsearchn(x2',XGL(end)); 
            figure(10); clf % sensitivity test changes
            hold on; grid on;
            set(gcf,'Position',[491 80 886 686]);
            set(gca,'FontSize',14,'linewidth',2,'fontweight','bold');
            xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)');
            legend('Location','east'); xlim([0 85]); ylim([min(hb0)-100 max(h)+100]);
            title(['SMR = + ',num2str(round(delta_SMR.*3.1536e7,1)),'m/a, SMB = + ',...
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
                xlim([x(gl)/10^3-5 x(c)/10^3+5]); ylim([min(hb0)-100 300]);
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
            b=b+1;
        end

    %catch
    %    disp(['iteration ',num2str(j),' failed']);
    %end
    
end

%% 2. conduct sensitivity tests for ocean thermal forcing TF independently
%   Simulate increase in submarine melting due to increased ocean thermal
%   forcing up to 1^o C
% -------------------------------------------------------------------------
%   Submarine melt rate relationship from Slater et al. (2020):
%   m' = (3*10−4 * h * q^0.39 + 0.15) * TF^1.18     (1)
%   where   h = grounding line depth [m]
%           q = mean annual subglacial runoff normalized by calving front
%           area [m/d]
%           TF = ocean thermal forcing [^oC]
% -------------------------------------------------------------------------

close all; 

saveFinal = 1;         % = 1 to save final conditions in homepath/3_sensitivityTests/results/
plotTimeSteps = 1;     % = 1 to plot geometry, speed, cf/gl positions every decade
plotMisfits = 0;       % = 1 to plot misfit with 2018 conditions
plotClimateParams = 1; % = 1 to plot climate parameters

% load no change conditions
load('2100_noChange.mat'); % load no change variables
    
% set up changes in TF
delta_TF0 = (0:0.1:1); % ^oC change in TF 

% define time stepping (s)
dt = 0.01*3.1536e7;
t_start = 0*3.1536e7;
t_end = 91*3.1536e7;

% define thermal forcing
TF0 = 0.2; % ^oC - estimated from Larsen B icebergs

% initialize ice flux
Fgl=zeros(1,length(t)); % ice mass flux across the grounding line

% estimate initial melt rate using Eqn from Slater et al. (2020):
mdot0 = (3*10^-4*-hb0(gl0)*((sum(RO0(1:gl0)))*86400)^0.39 + 0.15)*TF0^1.18/86400; % m/s

% loop through scenarios
for j=1:length(delta_TF0)

    delta_TF = delta_TF0(j);

    %try
        % run flowline model
        for i=1:length(t)
           
            % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
            Rxx = 2*nthroot(dUdx./(E.*A),n); % resistive stress (Pa)
            crev_s = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*FWD); % surface crevasse penetration depth (m)
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
                if xcf_s<0; xcf_s=NaN; end
                xcf=xcf_s;
                if xcf<20e3 || xcf > 100e3 || isnan(xcf)
                    xcf = x(dsearchn(x',x(c)));
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
            beta = interp1(x0,beta0,xn,'linear','extrap'); beta(gl+1:end)=0;
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
            [U,dUdx,Td,Tlatb,Tlon,vm] = U_convergence(x,U,U0,dUdx,H,h,A,E,N,W,dx,c,n,m,beta,rho_i,rho_sw,g,sigma_b,i);
            
            % calculate ice flux
            F = U.*H.*W; % ice flux (m^3 s^-1)
            F(isnan(F))=0;
            F(1)=F(2)+F0;
            
            % implement SMB, SMR, delta_SMB, & delta_TF
            if t(i)/3.1536e7<10 % use original SMB & SMR for first 10 model years
                SMR = zeros(1,c);
                % Define SMB and estimated runoff
                SMB = interp1(x0,SMB0,x);
                RO = interp1(x0,RO0,x);
                % no additional melt
                delta_mdot = 0;                
                % plot
                if i==1 && plotClimateParams
                    figure(2); clf
                    set(gcf,'position',[200 300 1000 500]);
                    subplot(2,2,1); hold on; 
                        set(gca,'fontsize',12,'fontweight','bold','linewidth',2);
                        plot(x/10^3,SMB.*3.1536e7,'color',col(i,:),'linewidth',2);
                        xlabel('km along centerline'); ylabel('m a^{-1}'); grid on;
                        title('SMB');
                    subplot(2,2,2); hold on; 
                        set(gca,'fontsize',12,'fontweight','bold','linewidth',2);
                        plot(x/10^3,SMR.*3.1536e7,'color',col(i,:),'linewidth',2);
                        xlabel('km along centerline'); ylabel('m a^{-1}'); grid on;
                        title('SMR');
                    subplot(2,2,3); hold on;
                        set(gca,'fontsize',12,'fontweight','bold','linewidth',2);
                        plot(x/10^3,RO*3.1536e7,'color',col(i,:),'linewidth',2);
                        xlabel('km along centerline'); ylabel('m a^{-1}'); grid on;
                        title('Runoff');
                    subplot(2,2,4); hold on;
                        set(gca,'fontsize',12,'fontweight','bold','linewidth',2);
                        plot(t(i)/3.1536e7+2009,TF,'.','markersize',10,'color',col(i,:));
                        xlabel('Year'); ylabel('^oC'); grid on;
                        title('TF');
                end                            
            elseif t(i)/3.1536e7>=10 % implement changes after 10 model years
                SMB = interp1(x0,SMB0,x);
                RO = (interp1(x0,SMB0,x)-SMB)+interp1(x0,RO0,x);
                % calculate additional melt due to the increase in subglacial discharge
                TFi = delta_TF/(2100-2019)*(t(i)/3.1536e7-10); % total increase in TF from 2019 
                delta_mdot = ((3*10^-4*-hb(gl)*((sum(RO(1:gl))*86400)^0.39) + 0.15)*((TFi+TF0)^1.18))/86400-mdot0; % m/s
                SMR = zeros(1,c);                
                SMR(gl+1:c) = (SMR0-delta_mdot)/(SMR_mean_fit.a+1)*feval(SMR_mean_fit,x(gl+1:c)-x(gl));
                % plot
                if mod(t(i)/3.1536e7,10)==0 && plotClimateParams
                    figure(2);
                    subplot(2,2,1);
                        plot(x/10^3,SMB.*3.1536e7,'color',col(i,:),'linewidth',2);
                    subplot(2,2,2);
                        plot(x/10^3,SMR.*3.1536e7,'color',col(i,:),'linewidth',2);
                    subplot(2,2,3);
                        plot(x/10^3,RO*3.1536e7,'color',col(i,:),'linewidth',2);
                    subplot(2,2,4);
                        plot(t(i)/3.1536e7+2009,TFi+TF0,'.','markersize',10,'color',col(i,:));
                end
            end
            
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

        % plot results
        if t(i)==t_end
            h2=h; H2=H; x2=x; c2=c; gl2=gl; U2=U; DFW2=FWD; Fgl2=Fgl; XCF2=XCF; XGL2=XGL; % store final geometry & speed
            figure(10); clf % sensitivity test changes
            hold on; grid on;
            set(gcf,'Position',[491 80 886 686]);
            set(gca,'FontSize',14,'linewidth',2,'fontweight','bold');
            xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)');
            legend('Location','east'); xlim([0 85]); ylim([min(hb0)-100 max(h)+100]);
            title(['TF = + ',num2str(TFi),'^oC']);
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
                xlim([x(gl)/10^3-5 x(c)/10^3+5]); ylim([min(hb0)-100 300]);
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

        % save geometry & speed
        if saveFinal
            cd([homepath,'scripts/3_sensitivityTests/results/2_TF/']);
            if delta_TF==0  
                save('TF0_geom.mat','h2','H2','c2','U2','gl2','x2','FWD2','Fgl2','XGL2','XCF2');
            else
                fileName = ['TF',num2str(delta_TF),'_geom.mat'];
                save(fileName,'h2','H2','c2','U2','gl2','x2','FWD2','Fgl2','XGL2','XCF2');
            end
            disp('geometry saved.');
        else
            disp('geometry not saved.');
        end
        
    %catch
    %    disp(['iteration ',num2str(j),' failed']);
    %end
    
end

%% 3. conduct sensitivity tests for SMB + enhanced SMR:
%   Simulate enhanced runoff due to increased surface melt -> 
%   increased meltwater discharge / convection due to plume
% -------------------------------------------------------------------------
%   Submarine melt rate relationship from Slater et al. (2020):
%   m' = (3*10−4 * h * q^0.39 + 0.15) * TF^1.18     (1)
%   where   h = grounding line depth [m]
%           q = mean annual subglacial runoff normalized by calving front
%           area [m/d]
%           TF = ocean thermal forcing [^oC]
% -------------------------------------------------------------------------

close all; 

saveFinal = 1;     % = 1 to save final geometry and speed 
plotTimeSteps = 0; % = 1 to plot geometry, speed, cf/gl positions every decade
plotClimateParams = 1; % = 1 to plot climate parameters

% load no change conditions
cd([homepath,'inputs-outputs/']);
load('2100_noChange.mat'); % load no change variables
    
% set up changes in SMB
delta_SMB0 = (0:-1:-10)./3.1536e7; % m/s change in SMB at the calving front (used to increase gradient)

% define time stepping (s)
dt = 0.01*3.1536e7;
t_start = 0*3.1536e7;
t_end = 91*3.1536e7;
t = (t_start:dt:t_end);

% define thermal forcing
TF0 = 0.2; % ^oC - estimated from Larsen B icebergs

% initialize ice flux
Fgl=zeros(1,length(t)); % ice mass flux across the grounding line

% estimate initial melt rate using Eqn from Slater et al. (2020):
mdot0 = (3*10^-4*-hb0(gl0)*((sum(RO0(1:gl0)))*86400)^0.39 + 0.15)*TF0^1.18/86400; % m/s

% loop through scenarios
for j=1:length(delta_SMB0)

    % initialize variables
    x=x0; U=U0; W=W0; gl=gl0; dUdx=dUdx0; A=A0; h=h0; hb=hb0; H=H0; beta=beta0; FWD=FWD0; 
    
    XCF = NaN*ones(1,length(t)); XGL = NaN*ones(1,length(t)); % store xcf and xgl over time

    delta_SMB = delta_SMB0(j);

    %try
        % run flowline model
        for i=1:length(t)
           
            % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
            Rxx = 2*nthroot(dUdx./(E.*A),n); % resistive stress (Pa)
            crev_s = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*FWD); % surface crevasse penetration depth (m)
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
                if xcf_s<0; xcf_s=NaN; end
                xcf=xcf_s;
                if xcf<20e3 || xcf > 100e3 || isnan(xcf)
                    xcf = x(dsearchn(x',x(c)));
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
            beta = interp1(x0,beta0,xn,'linear','extrap'); beta(gl+1:end)=0;
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
            [U,dUdx,Td,Tlatb,Tlon,vm] = U_convergence(x,U,U0,dUdx,H,h,A,E,N,W,dx,c,n,m,beta,rho_i,rho_sw,g,sigma_b,i);
            
            % calculate ice flux
            F = U.*H.*W; % ice flux (m^3 s^-1)
            F(isnan(F))=0;
            F(1)=F(2)+F0;
            
            % implement SMB, SMR, delta_SMB, & delta_SMR
            if t(i)/3.1536e7<10 % use original SMB & SMR for first 10 model years
                SMR = zeros(1,c);
                % Define SMB and estimated runoff
                SMB = interp1(x0,SMB0,x);
                RO = interp1(x0,RO0,x);
                % no additional melt
                delta_mdot = 0;                
                % use the Larsen C mean melt rate profile to scale SMR
                % using the max initial SMR
                SMR(gl+1:c) = (SMR0+delta_mdot)/(SMR_mean_fit.a+1)*feval(SMR_mean_fit,x(gl+1:c)-x(gl+1));
                % plot
                if i==1 && plotClimateParams
                    figure(2); clf
                    set(gcf,'position',[200 300 1000 500]);
                    subplot(2,2,1); hold on; 
                        set(gca,'fontsize',12,'fontweight','bold','linewidth',2);
                        plot(x/10^3,SMB.*3.1536e7,'color',col(i,:),'linewidth',2);
                        xlabel('km along centerline'); ylabel('m a^{-1}'); grid on;
                        title('SMB');
                    subplot(2,2,2); hold on; 
                        set(gca,'fontsize',12,'fontweight','bold','linewidth',2);
                        plot(x/10^3,SMR.*3.1536e7,'color',col(i,:),'linewidth',2);
                        xlabel('km along centerline'); ylabel('m a^{-1}'); grid on;
                        title('SMR');
                    subplot(2,2,3); hold on;
                        set(gca,'fontsize',12,'fontweight','bold','linewidth',2);
                        plot(x/10^3,RO*3.1536e7,'color',col(i,:),'linewidth',2);
                        xlabel('km along centerline'); ylabel('m a^{-1}'); grid on;
                        title('Runoff');
                    subplot(2,2,4); hold on;
                        set(gca,'fontsize',12,'fontweight','bold','linewidth',2);
                        plot(t(i)/3.1536e7+2009,FWD,'.','markersize',10,'color',col(i,:));
                        xlabel('Year'); ylabel('m'); grid on;
                        title('FWD');
                end                            
            elseif t(i)/3.1536e7>=10 % implement changes after 10 model years
                delta_SMBi = delta_SMB/(2100-2019)*(t(i)/3.1536e7-10); % total increase in smb from 2019 rate 
                SMB = interp1(x0,SMB0,x);
                for k=1:c
                    SMB(k) = SMB(k)+delta_SMBi*(h0(1)-h(k))/(h0(1)-h0(c0)); 
                end
                RO = (interp1(x0,SMB0,x)-SMB)+interp1(x0,RO0,x);
                % calculate additional melt due to the increase in subglacial discharge
                delta_mdot = ((3*10^-4*-hb(gl)*((sum(RO(1:gl))*86400)^0.39) + 0.15)*((TF0)^1.18))/86400-mdot0; % m/s
                SMR = zeros(1,c);                
                SMR(gl+1:c) = (SMR0-delta_mdot)/(SMR_mean_fit.a+1)*feval(SMR_mean_fit,x(gl+1:c)-x(gl));
                % plot
                if mod(t(i)/3.1536e7,10)==0 && plotClimateParams
                    figure(2);
                    subplot(2,2,1);
                        plot(x/10^3,SMB.*3.1536e7,'color',col(i,:),'linewidth',2);
                    subplot(2,2,2);
                        plot(x/10^3,SMR.*3.1536e7,'color',col(i,:),'linewidth',2);
                    subplot(2,2,3);
                        plot(x/10^3,RO*3.1536e7,'color',col(i,:),'linewidth',2);
                    subplot(2,2,4);
                        plot(t(i)/3.1536e7+2009,FWD,'.','markersize',10,'color',col(i,:));
                end
            end
            % save mean final SMB
            if t(i)==t_end
                SMB_mean(j) = mean(SMB,'omitnan');
            end
            if t(i)==t_end && j==length(delta_SMB0)
                cd([homepath,'inputs-outputs/']);
                save('2100SMB_mean.mat','SMB_mean');
                disp('smb_mean saved');
            end
            
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

        % plot results
        if t(i)==t_end
            h2=h; H2=H; x2=x; c2=c; gl2=gl; U2=U; DFW2=FWD; Fgl2=Fgl; XCF2=XCF; XGL2=XGL; % store final geometry & speed
            figure(10); clf % sensitivity test changes
            hold on; grid on;
            set(gcf,'Position',[491 80 886 686]);
            set(gca,'FontSize',14,'linewidth',2,'fontweight','bold');
            xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)');
            legend('Location','east'); xlim([0 85]); ylim([min(hb0)-100 max(h)+100]);
            title(['SMR = + ',num2str(round(delta_SMR.*3.1536e7,1)),'m/a, SMB = + ',...
                num2str(round(delta_SMB*3.1536e7,1)),'m/a, FWD = ',num2str(FWD),'m']);
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
                xlim([x(gl)/10^3-5 x(c)/10^3+5]); ylim([min(hb0)-100 300]);
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

        % save geometry & speed
        if saveFinal
            cd([homepath,'scripts/3_sensitivityTests/results/2_SMB+enhancedSMR/']);
            if delta_SMR==0 && delta_SMB==0 
                save(['TF0_SMB0_geom.mat'],'h2','H2','c2','U2','gl2','x2','FWD2','Fgl2','XGL2','XCF2');
            else
                fileName = ['SMB',num2str(delta_SMB),'_geom.mat'];
                save(fileName,'h2','H2','c2','U2','gl2','x2','FWD2','Fgl2','XGL2','XCF2');
            end
            disp('geometry saved.');
        else
            disp('geometry not saved.');
        end
        
    %catch
    %    disp(['iteration ',num2str(j),' failed']);
    %end
    
end

%% 4. conduct sensitivity tests for SMB + enhanced SMR + increased ocean thermal forcing: 
% Simulate an increase in surface melt, increased submarine melt due to
% increased discharge, and increased thermal forcing

close all; 

saveFinal = 1;     % = 1 to save final geometry and speed 
plotTimeSteps = 0; % = 1 to plot geometry, speed, cf/gl positions every decade
plotClimateParams = 0; % = 1 to plot climate parameters

% load no change conditions
cd([homepath,'inputs-outputs/']);
load('2100_noChange.mat'); % load no change variables
    
% set up changes in SMB and TF
delta_SMB0 = (0:-1:-10)./3.1536e7; % m/s change in SMB at the calving front (used to increase gradient)
delta_TF0 = 0:0.1:1; % ^oC change in thermal forcing

% define time stepping (s)
dt = 0.01*3.1536e7;
t_start = 0*3.1536e7;
t_end = 91*3.1536e7;
t = (t_start:dt:t_end);

% define thermal forcing
TF0 = 0.2; % ^oC - estimated from Larsen B icebergs

% initialize ice flux
Fgl=zeros(1,length(t)); % ice mass flux across the grounding line

% estimate initial melt rate using Eqn from Slater et al. (2020):
mdot0 = (3*10^-4*-hb0(gl0)*((sum(RO0(1:gl0)))*86400)^0.39 + 0.15)*TF0^1.18/86400; % m/s

% loop through scenarios
for j=1:length(delta_SMB0)

    % initialize variables
    x=x0; U=U0; W=W0; gl=gl0; dUdx=dUdx0; A=A0; h=h0; hb=hb0; H=H0; beta=beta0; FWD=FWD0; 
    
    XCF = NaN*ones(1,length(t)); XGL = NaN*ones(1,length(t)); % store xcf and xgl over time

    % implement changes to climate variables
    delta_SMB = delta_SMB0(j);
    delta_TF = delta_TF0(j);

    %try
        % run flowline model
        for i=1:length(t)

            % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
            Rxx = 2*nthroot(dUdx./(E.*A),n); % resistive stress (Pa)
            crev_s = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*FWD); % surface crevasse penetration depth (m)
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
                if xcf_s<0; xcf_s=NaN; end
                xcf=xcf_s;
                if xcf<20e3 || xcf > 100e3 || isnan(xcf)
                    xcf = x(dsearchn(x',x(c)));
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
            beta = interp1(x0,beta0,xn,'linear','extrap'); beta(gl+1:end)=0;
            x = xn; dx = dxn;

            % calculate surface elevation
            h = hb+H; % surface elevation (m a.s.l.)
            h(gl+1:c) = (1-rho_i/rho_sw).*H(gl+1:c); %adjust the surface elevation of ungrounded ice to account for buoyancy
            H(h<0)=0-hb(h<0); h(h<0)=0; % surface cannot go below sea level
            h(h-H<hb) = hb(h-H<hb)+H(h-H<hb); % thickness cannot go beneath bed elevation

            % plot geometry, speed, & grounding line and calving front positions
            col = parula(length(t)+20); % color scheme for plots
            if plotTimeSteps && t(i)==t_start
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
            elseif plotTimeSteps && mod(i-1,round(length(t)/10))==0 % display every length(t)/10
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
            
            % calculate ice flux
            F = U.*H.*W; % ice flux (m^3 s^-1)
            F(isnan(F))=0;
            F(1)=F(2)+F0;
            
            % implement SMB, SMR, delta_SMB, & delta_SMR
            if t(i)/3.1536e7<10 % use original SMB & SMR for first 10 model years
                SMR = zeros(1,c);
                % Define SMB and estimated runoff
                SMB = interp1(x0,SMB0,x);
                RO = interp1(x0,RO0,x);
                % no additional melt
                delta_mdot = 0;                
                % use the Larsen C mean melt rate profile to scale SMR
                % using the max initial SMR
                SMR(gl+1:c) = (SMR0-delta_mdot)/(SMR_mean_fit.a+1)*feval(SMR_mean_fit,x(gl+1:c)-x(gl));
                % plot
                if i==1 && plotClimateParams
                    figure(2); clf
                    set(gcf,'position',[200 300 1000 500]);
                    subplot(2,2,1); hold on; 
                        set(gca,'fontsize',12,'fontweight','bold','linewidth',2);
                        plot(x/10^3,SMB.*3.1536e7,'color',col(i,:),'linewidth',2);
                        xlabel('km along centerline'); ylabel('m a^{-1}'); grid on;
                        title('SMB');
                    subplot(2,2,2); hold on; 
                        set(gca,'fontsize',12,'fontweight','bold','linewidth',2);
                        plot(x/10^3,SMR.*3.1536e7,'color',col(i,:),'linewidth',2);
                        xlabel('km along centerline'); ylabel('m a^{-1}'); grid on;
                        title('SMR');
                    subplot(2,2,3); hold on;
                        set(gca,'fontsize',12,'fontweight','bold','linewidth',2);
                        plot(x/10^3,RO*3.1536e7,'color',col(i,:),'linewidth',2);
                        xlabel('km along centerline'); ylabel('m a^{-1}'); grid on;
                        title('Runoff');
                    subplot(2,2,4); hold on;
                        set(gca,'fontsize',12,'fontweight','bold','linewidth',2);
                        plot(t(i)/3.1536e7+2009,TF0,'.','markersize',10,'color',col(i,:));
                        xlabel('Year'); ylabel('m'); grid on;
                        title('TF');
                end                            
            elseif t(i)/3.1536e7>=10 % implement changes after 10 model years
                delta_SMBi = delta_SMB/(2100-2019)*(t(i)/3.1536e7-10); % total increase in smb from 2019 rate 
                SMB = interp1(x0,SMB0,x);
                for k=1:c
                    SMB(k) = SMB(k)+delta_SMBi*(h0(1)-h(k))/(h0(1)-h0(c0)); 
                end
                RO = (interp1(x0,SMB0,x)-SMB)+interp1(x0,RO0,x);
                % implement change in thermal forcing linearly over time
                delta_TFi = delta_TF/(2100-2019)*(t(i)/3.1536e7-10);
                % calculate additional melt due to the increase in subglacial discharge
                delta_mdot = ((3*10^-4*-hb(gl)*((sum(RO(1:gl))*86400)^0.39) + 0.15)*((TF0+delta_TFi)^1.18)/86400)-mdot0; % m/s
                % increase SMR linearly at each time increment to reach delta_smr by 2100
                SMR = zeros(1,c);
                SMR(gl+1:c) = (SMR0-delta_mdot)/(SMR_mean_fit.a+1)*feval(SMR_mean_fit,x(gl+1:c)-x(gl));
                % plot
                if mod(t(i)/3.1536e7,10)==0 && plotClimateParams
                    figure(2);
                    subplot(2,2,1);
                        plot(x/10^3,SMB.*3.1536e7,'color',col(i,:),'linewidth',2);
                    subplot(2,2,2);
                        plot(x/10^3,SMR.*3.1536e7,'color',col(i,:),'linewidth',2);
                    subplot(2,2,3);
                        plot(x/10^3,RO*3.1536e7,'color',col(i,:),'linewidth',2);
                    subplot(2,2,4);
                        plot(t(i)/3.1536e7+2009,TF0,'.','markersize',10,'color',col(i,:));
                end
            end
            % save mean final SMB
            if t(i)==t_end
                SMB_mean(j) = mean(SMB,'omitnan');
            end
            if t(i)==t_end && j==length(delta_SMB0)
                cd([homepath,'inputs-outputs/']);
                save('2100SMB_mean.mat','SMB_mean');
                disp('smb_mean saved');
            end
            
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

        % plot results
        if t(i)==t_end
            h2=h; H2=H; x2=x; c2=c; gl2=gl; U2=U; DFW2=FWD; Fgl2=Fgl; XCF2=XCF; XGL2=XGL; % store final geometry & speed
            figure(10); clf % sensitivity test changes
            hold on; grid on;
            set(gcf,'Position',[491 80 886 686]);
            set(gca,'FontSize',14,'linewidth',2,'fontweight','bold');
            xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)');
            legend('Location','east'); xlim([0 85]); ylim([min(hb0)-100 max(h)+100]);
            title(['SMR = + ',num2str(round(delta_SMR.*3.1536e7,1)),'m/a, SMB = + ',...
                num2str(round(delta_SMB*3.1536e7,1)),'m/a, FWD = ',num2str(FWD),'m']);
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
                xlim([x(gl)/10^3-5 x(c)/10^3+5]); ylim([min(hb0)-100 300]);
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

        % save geometry & speed
        if saveFinal
            cd([homepath,'scripts/3_sensitivityTests/results/3_SMB+enhancedSMR+SMR/']);
            if delta_DFW==0 && delta_SMR==0 && delta_SMB==0 
                save(['TF0_SMB0_geom.mat'],'h2','H2','c2','U2','gl2','x2','FWD2','Fgl2','XGL2','XCF2');
            elseif delta_DFW==0 
                fileName = ['TF',num2str(delta_TF),'_SMB',num2str(delta_SMB),'_geom.mat'];
                save(fileName,'h2','H2','c2','U2','gl2','x2','FWD2','Fgl2','XGL2','XCF2');
            else
                fileName = ['FWD_',num2str(FWD),'m_geom.mat'];
                save(fileName,'h2','H2','c2','U2','gl2','x2','FWD2','Fgl2','XGL2','XCF2');
            end
            disp('geometry saved.');
        else
            disp('geometry not saved.');
        end
        
    %catch
    %    disp(['iteration ',num2str(j),' failed']);
    %end
    
end
