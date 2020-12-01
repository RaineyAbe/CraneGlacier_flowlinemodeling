%% Glacier Flowline Model
% Script to tune the upper boundary and SMB to nearly match observations of
% ice surface elevation (and dHDt)
%
% Rainey Aberle
% Fall 2020
% Adapted from Ellyn Enderlin's flowline model demo code

clear all; close all;
warning off; % turn off warnings (velocity coefficient matrix is close to singular)

%% define time and space independent variables
    
dx0 = 150; % desired grid spacing (m)
dx=dx0;
            
use_binavg = 1;     % = 1 to use average within bins proportional to dx0
plot_stresses = 0;  % = 1 to plot stress balance variables

% define home path in directory
    homepath = '/Users/raineyaberle/Desktop/Research/CraneGlacier_modeling/';
    cd([homepath,'inputs-outputs']);
    addpath([homepath,'scripts/tuningCalving']); % add path to U_convergence

% Load Crane Glacier initialization variables
    load('Crane_flowline_initialization.mat');
    n = find(isnan(W0),1,'first');
    W0(n:end) = W0(n-1).*ones(1,length(W0(n:end))); clear n
    A0(end+1:length(x0)) = A0(end).*ones(1,length(x0)-length(A0)); 
    
% Load observations of dH to help tune SMB
    %dH_obs = load('dHdt.mat').dHdt.dH_total; % (m) total change in thickness 2009-2018
% Load observations of terminus position to help tune calving parameter
    term = load('Crane_TerminusPosition_2002-2019.mat').term;
    for i=1:length(term)
        termx_obs(i) = term(i).x;
        termDate_obs(i) = term(i).decidate;
    end
    % start in 2009
    termx_obs(1:4)=[]; termDate_obs(1:4)=[];
    % fit a quadratic to the terminus positions
    termx_obs = feval(fit(termDate_obs',termx_obs','poly2'),termDate_obs');
    clear term
    
% densities and g
    rho_i = 917; % ice density (kg m^-3)
    rho_sw = 1028; % ocean water density (kg m^-3)
    rho_fw = 1000; % fresh water density (kg m^-3)
    g = 9.81; % acceleration (m s^-2)

% time stepping (s)
    dt = 0.1*3.1536e7; 
    t_start = 0*3.1536e7; 
    t_end = 10*3.1536e7;    
    t = (t_start:dt:t_end);

% stress parameters (unitless)
    m = 1; % basal sliding exponent
    n = 3; % flow law exponent
    E = 1; % enhancement factor
    
% maximum thickness cut-off to check for instability
    H_max = 2000; %maximum thickness (m)

% regrid the initialization data to match the desired grid spacing, 
    xi = x0; % rename initialization distance vector
    L = 70e3; % length of glacier (m)
    xi = 0:dx0:L; % desired distance vector (m from ice divide)  
    
    % calving front location
    c = dsearchn(transpose(xi),x0(term_2009.x)); % 2009 terminus location (index)
         
    % If the desired grid spacing is smaller than the original, use the
    % interp1 function to determine each spatial vector.
    % Otherwise, take the average within each bin at every point for each
    % spatial variable. 
    if length(xi)>length(x0)
        hi = interp1(x0,h0,xi);
        hbi = interp1(x0,hb0,xi);
        Wi = interp1(x0,W0,xi);
        Ui = interp1(x0,U0,xi);
        Ai = interp1(x0,A0,xi);
        betai = interp1(x0,beta0,xi);            
    else
        % Use a staggered grid for bin averages
        xm = (xi(1:end-1)+xi(2:end))./2;
        for k=1:length(xi)
            if k==1
              hi(k) = mean(h0(1:dsearchn(x0',xm(k))));
              hbi(k) = mean(hb0(1:dsearchn(x0',xm(k))));
              Wi(k) = mean(W0(1:dsearchn(x0',xm(k))));
              Ui(k) = mean(U0(1:dsearchn(x0',xm(k))));
              Ai(k) = mean(A0(1:dsearchn(x0',xm(k))));
              betai(k) = mean(beta0(1:dsearchn(x0',xm(k))));               
            elseif k==length(xi)
              hi(k) = mean(h0(dsearchn(x0',xm(k-1)):c));
              hbi(k) = mean(hb0(dsearchn(x0',xm(k-1)):c));
              Wi(k) = mean(W0(dsearchn(x0',xm(k-1)):c));
              Ui(k) = mean(U0(dsearchn(x0',xm(k-1)):c));
              Ai(k) = mean(A0(dsearchn(x0',xm(k-1)):c));
              betai(k) = mean(beta0(dsearchn(x0',xm(k-1)):c));              
            else
              hi(k) = mean(h0(dsearchn(x0',xm(k-1)):dsearchn(x0',xm(k))));
              hbi(k) = mean(hb0(dsearchn(x0',xm(k-1)):dsearchn(x0',xm(k))));
              Wi(k) = mean(W0(dsearchn(x0',xm(k-1)):dsearchn(x0',xm(k))));
              Ui(k) = mean(U0(dsearchn(x0',xm(k-1)):dsearchn(x0',xm(k))));
              Ai(k) = mean(A0(dsearchn(x0',xm(k-1)):dsearchn(x0',xm(k))));
              betai(k) = mean(beta0(dsearchn(x0',xm(k-1)):dsearchn(x0',xm(k))));                            
            end
        end
    end 
    
    dUdxi = [(Ui(2:end)-Ui(1:end-1))./(xi(2:end)-xi(1:end-1)) 0]; % strain rate
    Hi = hi-hbi; % thickness (m)    

    % add the calving front conditions for each spatial variable
    hi(c+1:length(xi)) = zeros(1,length(hi(c+1:length(xi))));
    Ui(c+1:length(xi)) = zeros(1,length(Ui(c+1:length(xi))));
    Hi(c+1:length(xi)) = zeros(1,length(Hi(c+1:length(xi))));
    dUdxi(c+1:length(xi)) = zeros(1,length(dUdxi(c+1:length(xi))));
    betai(c+1:length(xi)) = zeros(1,length(betai(c+1:length(xi))));
    
    % find the location of the grounding line and the end of the ice-covered domain
    Hf = -(rho_sw./rho_i).*hbi; % flotation thickness (m)
    gl = find(Hf-Hi<0,1,'last')-1; %grounding line location 
    Hi(gl:end)=hi(gl:end)*rho_sw/(rho_sw-rho_i); % buoyant thickness using surface
    for i=1:length(xi) % thickness can't go beneath bed elevation
        if Hi(i) >= hi(i)-hbi(i)
            Hi(i) = hi(i)-hbi(i);
        end
    end    
    
    % add a dummy ice end (hi & Hi)
    for i=c+1:length(xi)
        %hi(i) = hi(i-1)-5; % decrease by 5m until at 0m  
        Hi(i) = Hi(i-1)-30; % decrease by 20m until 0 
        if Hi(i)>=hi(i)-hbi(i)
            Hi(i)=hi(i)-hbi(i); % can't go beneath bed elevation
        end
    end  
    hi(hi<0)=0; % surface can't go below sea level
    Hi(Hi<0)=0; % no negative thicknesses
    
    % find the end of the ice-covered domain (m along centerline)
    ice_end = (find(Hi<=0,1,'first')); 
    
    % extend other variables from c+1:ice_end (Ui,Ai)
    Ui(c+1:ice_end) = Ui(c).*ones(1,length(xi(c+1:ice_end)));
    Ai(c+1:ice_end) = Ai(c).*ones(1,length(xi(c+1:ice_end)));
    betai(gl:end) = zeros(1,length(betai(gl:end)));
    
    % use 90% the observed velocity at the upper bounds
    Ui(1:round(0.2*length(Ui))) = 0.9.*Ui(1:round(0.2*length(Ui)));
    
    % rename initial vectors
    x=xi; h=hi; hb=hbi; W=Wi; H=Hi; A=Ai; beta=betai; U=Ui; dUdx=dUdxi;
    
%% Run the flowline model   
    
% calving parameters
    Hc = 400; % m -> set the calving front to a default minimum ice thickness value
    fwd = 30;%5:5:100; % fresh water depth in crevasses (m)
    term_rmse = zeros(length(Hc),length(fwd),length(t));
    
for p=1:length(Hc)    
    for q=1:length(fwd)  
        
        % re-initialize variables
        x=xi; h=hi; hb=hbi; W=Wi; H=Hi; A=Ai; beta=betai; U=Ui; dUdx=dUdxi;
        % find the end of the ice-covered domain (m along centerline)
        ice_end = (find(Hi<=0,1,'first')); 
        % calving front location
        c = dsearchn(transpose(xi),x0(term_2009.x)); % 2009 terminus location (index)
        
        % close figures
        close all; 
        
        % continue to the next loop if there is an error
        try
            for i=1:length(t)

                % set up figures, plot geometries at t==0, then every t/10 iterations
                if t(i)==t_start
                    col = parula(length(t)+10); %Color scheme for plots
                    figure; % glacier geometry
                        hold on; grid on;
                        set(gcf,'Position',[0 50 500 400]);
                        set(gca,'FontSize',14,'linewidth',2,'fontweight','bold'); 
                        legend('Location','east'); xlim([0 65]); ylim([min(hb)-100 max(h)+100]);
                        title('a) Glacier Geometry'); 
                        xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)'); 
                        % ice surface
                        plot(x(1:c)./10^3,h(1:c),'color',col(i,:),'linewidth',2,'displayname','2009');
                        % calving front
                        plot(x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
                        % ice end
                        plot(x(c+1:ice_end)./10^3,h(c+1:ice_end),'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');            
                        plot(x(c+1:ice_end)./10^3,h(c+1:ice_end)-H(c+1:ice_end),'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');            
                        % floating bed
                        plot(x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off');
                        % bed elevation
                        plot(x/10^3,hb,'k','linewidth',2,'HandleVisibility','off'); 
                        % mean sea level
                        plot([x(1),x(end)]/10^3,[0,0],'k--','HandleVisibility','off'); 
                        % ice end            
                    figure; % ice speed
                        hold on; grid on; 
                        set(gcf,'Position',[500 50 500 400]);
                        set(gca,'FontSize',14,'linewidth',2,'fontweight','bold'); 
                        title('b) Ice Speed Profile');  
                        xlim([0 65]); ylim([0 2200]);
                        xlabel('Distance Along Centerline (km)'); ylabel('Speed (m yr^{-1})'); 
                        legend('Location','east'); 
                        % 1:c
                        plot(x(1:c)./10^3,U(1:c).*3.1536e7,'color',col(i,:),'linewidth',2,'displayname','2009');
                        % c:ice_end
                        plot(x(c:ice_end)./10^3,U(c:ice_end).*3.1536e7,'--','color',col(i,:),'linewidth',2,'HandleVisibility','off');            
                    figure; % terminus position
                        hold on; grid on; 
                        set(gcf,'Position',[1000 50 500 400]);
                        set(gca,'FontSize',14,'linewidth',2,'fontweight','bold');
                        title('c) Terminus Position'); ylabel('Year');
                        xlabel('Distance Along Centerline (km)'); 
                        xlim([42 52]); ylim([2008 2020]); legend('Location','east');
                        plot(termx_obs./10^3,termDate_obs,'.-k','markersize',10,'displayname','Observed');
                        % 2009 terminus position
                        plot(x(c)./10^3,2009,'.','markersize',10,'color',col(i,:),...
                            'linewidth',1.5,'displayname','Modeled'); 
                    if plot_stresses==1
                        figure; % fig4 = basal resistance and driving stress
                            set(gcf,'Position',[350 600 700 600]);
                    end
                else
                    figure(1); hold on; % Plot geometries every 10 time iterations
                    if mod(i-1,round(length(t)/10))==0
                        % ice surface
                        plot(x(1:c)/10^3,h(1:c),'-','color',col(i,:),'linewidth',2,'displayname',num2str(t(i)./3.1536e7+2009)); 
                        % calving front
                        plot(x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
                        % floating bed
                        plot(x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off');
                    else 
                        % ice surface
                        plot(x(1:c)/10^3,h(1:c),'-','color',col(i,:),'linewidth',2,'HandleVisibility','off'); 
                        % calving front
                        plot([x(c)/10^3 x(c)/10^3],[h(c) h(c)-H(c)],'-','color',col(i,:),'linewidth',2,'HandleVisibility','off'); 
                        % ice end
                        plot(x(c+1:ice_end)./10^3,h(c+1:ice_end),'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');            
                        plot(x(c+1:ice_end)./10^3,h(c+1:ice_end)-H(c+1:ice_end),'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');            
                        % floating bed
                        plot(x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off');
                    end 

                    figure(2); hold on; % Plot velocity every 10 time iterations
                    if mod(i-1,round(length(t)/10))==0
                        % 1:c
                        plot(x(1:c)/10^3,U(1:c).*3.1536e7,'-','Color',col(i,:),'linewidth',2,'DisplayName',num2str(t(i)./3.1536e7+2009)); hold on;
                        % c:ice_end
                        plot(x(c:ice_end)./10^3,U(c:ice_end).*3.1536e7,'--','color',col(i,:),'linewidth',2,'HandleVisibility','off');            
                    else 
                        % 1:c
                        plot(x(1:c)/10^3,U(1:c).*3.1536e7,'-','Color',col(i,:),'linewidth',2,'HandleVisibility','off'); hold on;           
                        % c:ice_end
                        plot(x(c:ice_end)./10^3,U(c:ice_end).*3.1536e7,'--','color',col(i,:),'linewidth',2,'HandleVisibility','off');                    
                    end

                    figure(3); hold on; % Plot terminus position every 10 time iterations
                    if mod(i-1,round(length(t)/10))==0
                        plot(x(c)./10^3,t(i)./3.1536e7+2009,'.','Color',col(i,:),'markersize',15,'linewidth',1.5,'displayname',num2str(t(i)./3.1536e7+2009)); hold on;
                    else 
                        plot(x(c)./10^3,t(i)./3.1536e7+2009,'.','Color',col(i,:),'markersize',10,'linewidth',1.5,'HandleVisibility','off'); hold on;           
                    end          
                end 

                % calculate RMSE
                term_rmse(p,q,i) = termRMSE(x,c,t(i),termx_obs,termDate_obs);

                % calculate the thickness required to remain grounded at each grid cell
                    Hf = -(rho_sw./rho_i).*hb; %flotation thickness (m)
                    % find the location of the grounding line and use a floating
                    % geometry from the grounding line to the calving front
                    gl = find(Hf-H<0,1,'last')-1; %grounding line location 
                    ice_end = (find(H<=0,1,'first')); %end of ice-covered domain
                    if isempty(ice_end) || ice_end>length(x)
                        ice_end = length(x);
                        disp('ice end criteria not met.');
                    end

                % calculate the glacier's surface elevation and slope
                    h = hb+H; %h = surface elevation (m a.s.l.)
                    h(gl:length(x)) = (1-rho_i/rho_sw).*H(gl:length(x)); %adjust the surface elevation of ungrounded ice to account for buoyancy
                    dhdx = [(h(2:end)-h(1:end-1))./(x(2:end)-x(1:end-1)) 0]; %surface slope (unitless)

                % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
                    Rxx = 2*nthroot((dUdx./(E*A(1:length(dUdx)))),n); %resistive stress (Pa)
                    crev = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*fwd(q)); %crevasse penetration depth (m)
                    c = find(h(1:ice_end-1)-crev(1:ice_end-1)<=0,1,'first'); %calving front located where the inland-most crevasse intersects sea level
                    %if the crevasses never intersect sea level
                    if isempty(c) == 1
                        c = find(H<Hc(p),1,'first'); %set the calving front to a default minimum ice thickness value
                    end
                    if isempty(c)==1
                        c=length(x);
                        disp('calving criteria not met');
                    end
                    % if the crevasses first intersect sea level inland of the grounding line 
                    if c <= gl
                        c = gl; %set the grounding line as the calving front
                    end
                    % use observed terminus position for first time increment
                    if i==1
                        c=dsearchn(transpose(x),x0(term_2009.x));            
                    end 

                %calculate the effective pressure (ice overburden pressure minus water
                %pressure) assuming an easy & open connection between the ocean and
                %ice-bed interface
                sl = find(hb<=0,1,'first'); %find where the glacier base first drops below sea level
                N_ground = rho_i*g*H(1:sl); %effective pressure where the bed is above sea level (Pa)
                N_marine = rho_i*g*H(sl+1:length(x))+(rho_sw*g*hb(sl+1:length(x))); %effective pressure where the bed is below sea level (Pa)
                N = [N_ground N_marine];
                N(N<0)=1; %cannot have negative values

                % Solve for new velocity
                [U,dUdx,vm,T] = U_convergence(x,U,U0,dUdx,dhdx,H,A,E,N,W,dx,c,ice_end,n,m,beta,rho_i,rho_sw,g); 

                % calculate ice flux
                F = U.*H.*W; % ice flux (m^3 s^-1)
                F(isnan(F))=0;

                % calculate the  change in ice thickness from continuity
                dHdt = -(1./W).*gradient(F,x);
                dH = dHdt.*dt;

                % surface mass balance
                yr = round(t(i)/3.1536e7)+1;
                clear smb sigma_smb % clear to avoid changing size with changing x

                % interpolate smb0 to centerline, add tributary flux Q0 to smb
                smb = [interp1(x0,smb0(yr).smb_linear'+Q0,x)]./3.1536e7; % m/s
                    smb(c+1:length(x)) = zeros(1,length(smb(c+1:length(x)))); % zeros past c
                sigma_smb = [interp1(x0,smb0(yr).sigma_smb+Q0_err,x)]./3.1536e7; % m/s
                    sigma_smb(c+1:length(x)) = zeros(1,length(smb(c+1:length(x)))); % zeros past c

                % adjust smb to minimize misfit of surface observations 
                smb = smb-0.1e-5; % m/s
                smb(1:30) = smb(1:30)-0.12e-5; 
                smb(50:70) = smb(50:70)+0.05e-5;
                smb(50:100) = smb(50:100)+0.08e-5; 
                smb(115:290) = smb(115:290)-0.06e-5; 

                % new thickness (change from dynamics, SMB, & submarine melting)
                Hn = H+dH+(smb.*dt); 
                Hn(Hn < 0) = 0; % remove negative values 
                H = Hn; %set as the new thickness value

                % stop the model if it behaves unstably (monitored by ice thickness)
                if max(H) > H_max
                    disp(['Adjust dt']);
                    break
                end

                % find the precise location of the grounding line (where H=Hf)
                %xf = interp1((H(1:ice_end)-Hf(1:ice_end)),x(1:ice_end),0,'spline','extrap'); 
                xf = find(Hf-H>0,1,'first')-1;

                %adjust the grid spacing so the grounding line is continuously tracked
                xl = round(xf/dx0); %number of ideal grid spaces needed to reach the grounding line
                dx = xf/xl; %new grid spacing (should be ~dx0)
                xn = 0:dx:L; %new distance vector    

                %adjust the space-dependent variables to the new distance vector
                hb = interp1(x0,hb0,xn);
                W = interp1(x0,W0,xn);
                H = interp1(x,H,xn,'linear','extrap'); %ice thickness (m)
                Hf = interp1(x,Hf,xn,'linear','extrap');
                U = interp1(x(1:c),U(1:c),xn,'linear','extrap'); %speed (m s^-1)
                A = interp1(x,A,xn,'linear','extrap'); %rate factor (Pa^-n s^-1)
                beta = interp1(x,beta,xn,'linear','extrap'); % basal roughness factor

                %find the location of the grounding line and end of the ice-covered domain for the adjusted data
                gl = find(Hf-H<0,1,'last');
                ice_end = (find(H<=0,1,'first'));
                if isempty(ice_end) || ice_end>length(xn)
                    ice_end = length(xn);
                    disp('ice end criteria not met.')
                end

                %rename the distance vector
                x = xn; %distance from the divide (m)

                %calculate the new surface elevation and slope
                h = hb+H; %grounded ice surface elevation (m a.s.l.)
                h(gl:length(x)) = (1-rho_i/rho_sw).*H(gl:length(x)); %floating ice surface elevation (m a.s.l.)

                % calculate new strain rate
                dUdx = [(U(2:end)-U(1:end-1))./(x(2:end)-x(1:end-1)) 0]; % strain rate

                % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
                    Rxx = 2*nthroot((dUdx./(E*A(1:length(dUdx)))),n); %resistive stress (Pa)
                    crev = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*fwd(q)); %crevasse penetration depth (m)
                    c = find(h(1:ice_end-1)-crev(1:ice_end-1)<=0,1,'first'); %calving front located where the inland-most crevasse intersects sea level
                    %if the crevasses never intersect sea level
                    if isempty(c) == 1
                        c = find(H<Hc(p),1,'first'); %set the calving front to a default minimum ice thickness value
                    end
                    if isempty(c)==1
                        c=length(x);
                        disp('calving criteria not met');
                    end
                    % if the crevasses first intersect sea level inland of the grounding line 
                    if c <= gl
                        c = gl; %set the grounding line as the calving front
                    end
                    % use observed terminus position for first time increment
                    if i==1
                        c=dsearchn(transpose(x),x0(term_2009.x));            
                    end 

                % Calculate resistive and driving stresses
                    % Eqn. (2) in Flowline Model User Guide (Enderlin, 2013)
                    if plot_stresses==1 && t(i)==t(end)

                        % Calculate basal resistance
                            Rb = -beta(1:c).*N(1:c).*(U(1:c).^(1/m)); 

                        % Calculate longitudinal and lateral stresses (Pa)
                        Rlon_piece = zeros(1,c); R_lat = zeros(1,c); 
                        for k=1:c
                            Rlon_piece(k) = H(k)*vm(k)*dUdx(k); % Calculate inner matrix before taking gradient
                            Rlat(k) = -2*H(k)/W(k)*((5*U(k)/(E*A(k)*W(k)))^(1/n));
                        end 
                            % Use a forward difference to calculate d/dx(R_lon_piece)
                            Rlon = zeros(1,c);
                            for k=1:length(x)
                                if k==length(x)
                                    Rlon(k) = NaN;
                                else
                                    Rlon(k) = (Rlon_piece(k+1)-Rlon_piece(k))/(x(k+1)-x(k));
                                end 
                            end

                        % Calculate total resistive stress
                        R_tot = Rb+Rlon+Rlat; 

                        % Plot results onto figure 5
                        figure(4); clf
                        subplot(3,2,1); 
                        hold on;
                        plot(x(1:c)/10^3,Rb(1:c),'-b','linewidth',2,'displayname','BR');
                        grid on; set(gca,'FontSize',12,'fontname','Arial');        
                        title('R_{basal}'); ylabel('(Pa)');

                        subplot(3,2,2); 
                        hold on; 
                        plot(x(1:c-1)/10^3,T(1:c-1),'-b','linewidth',2);
                        title('T_d');
                        grid on; set(gca,'FontSize',12,'fontname','Arial');

                        subplot(3,2,3); 
                        hold on; 
                        plot(x(1:c-1)/10^3,Rlon(1:c-1),'-b','linewidth',2);
                        title('R_{longitudinal}'); ylabel('(Pa)');
                        grid on; set(gca,'FontSize',12,'fontname','Arial');

                        subplot(3,2,4); 
                        hold on; 
                        plot(x(1:c-1)/10^3,Rlat(1:c-1),'-b','linewidth',2);
                        title('R_{lateral}');
                        grid on; set(gca,'FontSize',12,'fontname','Arial');

                        subplot(3,2,[5 6]);
                        hold on; 
                        plot(x(1:c-1)/10^3,R_tot(1:c-1),'-m','linewidth',2,'displayname','R_tot');
                        title('R_{total}'); 
                        grid on; set(gca,'FontSize',12,'fontname','Arial');
                        xlabel('Distance Along Centerline (km)'); ylabel('(Pa)');

                    end 

            end 
        catch
            disp(['iteration p=',num2str(p),' q=',num2str(q),' broke the loop.']);
            term_rmse(p,q,i:end) = NaN; % insert NaNs for remaining times
            continue; % continue to the next loop
        end
        
    end 
end

term_rmse(term_rmse==0)=NaN;
term_Cts = zeros(length(Hc),length(fwd)); % # of successful iterations per combo
term_rmse_mean = zeros(length(Hc),length(fwd)); % 

for p=1:length(Hc)
    for q=1:length(fwd)
        
        % Count the number of time iterations in each combo before crashing
        term_Cts(p,q) = length(find(~isnan(term_rmse(p,q,:))));
        
        % Calculate mean RMSE for each combo of Hc and fwd
        term_rmse_mean(p,q) = nanmean(term_rmse(p,q,:));
        
    end
end

% Display results
figure(6); clf; hold on; 
plot(fwd,term_rmse_mean,'.-k','markersize',10,'linewidth',2);
set(gca,'fontsize',14,'linewidth',2);
xlabel('fwd (m)'); ylabel('Mean error');
yyaxis right
plot(fwd,term_Cts,'.-r','markersize',5,'linewidth',1);
ylabel('# of points');

% Save image
cd([homepath,'scripts/tuningCalving/']);
saveas(gcf,'fwdError.png','png');



