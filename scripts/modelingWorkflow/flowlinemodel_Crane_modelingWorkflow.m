%% Glacier Flowline Model: Modeling workflow for Crane Glacier, Antarctic Peninsula
% Rainey Aberle
% Last edited: Spring 2021
% Adapted from Ellyn Enderlin's flowline model package (Enderlin et al.,
% 2013)
%
% Workflow: Run step 0 before any other step to initialize variables.
%   0. Define time and space independent variables by loading the flowline 
%       initialization file and regridding to the desired grid spacing. 
%   1. SKIP IF OUTPUT ALREADY SAVED. Run flowline model for 100 years 
%       to reach near steady-state conditions using 2009 observations. 
%   2. Tune the inner boundary flux F0 using observations of surface speed 
%       and elevation. 
%   3. Tune the backstress sigma_b using the resistive stress term to 
%       minimize misfit between modeled and observed terminus position 
%       to account for sea ice or sikkusak (Nick et al., 2010)
%   4. Tune the enhancement factor E along the centerline to minimize
%       misfit between modeled and observed ice surface elevation.
%   5. Tune the fresh water depth in crevasses fwd to minimize the misfit
%       between modeled and observed terminus positions 2009-2019.
%   6. Run sensitivity tests for surface mass balance (SMB) & submarine
%       melting rate (SMR). 

%% 0. define time and space independent variables
  
clear all; close all;
warning off; % turn off warnings (velocity coefficient matrix is close to singular)

dx0 = 250; % desired grid spacing (m)
dx=dx0;
            
% define home path in directory and add necessary paths
    homepath = '/Users/raineyaberle/Desktop/Research/CraneGlacier_modeling/';
    cd([homepath,'scripts/modelingWorkflow']);
    addpath([homepath,'inputs-outputs']);
    addpath([homepath,'scripts/modelingWorkflow']); % add path to U_convergence

% Load Crane Glacier initialization variables
    load('Crane_flowline_initialization.mat');

% Load observed conditions
    % dH
    dH_obs = load('dHdt_2009-2018.mat').dHdt.dH_total; % (m) total change in thickness 2009-2018
    % terminus position 
    term = load('Crane_TerminusPosition_2002-2019.mat').term;
    for i=1:length(term)
        termx_obs(i) = term(i).x;
        termDate_obs(i) = term(i).decidate;
    end
    % fit a quadratic to the terminus positions
    termx_obs = feval(fit(termDate_obs',termx_obs','poly2'),termDate_obs');
    term_obs = interp1(termDate_obs',termx_obs,2009:2017);
    clear term 
    % ice speed
    U_obsi = load('Crane_CenterlineSpeeds_2007-2017.mat').U;
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
    m = 1; % basal sliding exponent
    n = 3; % flow law exponent

% calving parameters
    Hc = 400; % m -> set the calving front to a default minimum ice thickness value
    fwd = 35; % fresh water depth in crevasses (m)
    sigma_b = 10e3; % back pressure (kPa) - similar to that employed at Helheim Glacier (Nick et al., 2009)
     
% maximum & minimum thickness cut-off to check for instability
    H_max = 2000; % maximum thickness (m)
    H_min = 100;  % minimum thickness (m)

% regrid the initialization data to match the desired grid spacing 
    L = 70e3; % length of glacier (m)
    xi = 0:dx0:L; % desired distance vector (m from ice divide)  
         
    % If the desired grid spacing is smaller than the original, use the
    % interp1 function to determine each spatial vector.
    % Otherwise, use the average value within each bin for spatial variables. 
    ci = dsearchn(xi',x0(c0)); c=ci; % regrid the calving front location
    if (xi(2)-xi(1))>(x0(2)-x0(1))
        hi = interp1(x0,h0,xi);
        hbi = interp1(x0,hb0,xi);
        Wi = interp1(x0,W0,xi);
        Ui = interp1(x0,U0,xi);
        Ai = interp1(x0,A0_adj(1,:),xi);
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
    
    dUdxi = [(Ui(2:end)-Ui(1:end-1))./(xi(2:end)-xi(1:end-1)) 0]; % strain rate (1/s)
    dH_obs = interp1(x0,dH_obs,xi); % regrid observed dH (m)
    Hi = hi-hbi; % thickness (m)    
    
    % find the location of the grounding line and the end of the ice-covered domain
    Hf = -(rho_sw./rho_i).*hbi; % flotation thickness (m)
    gl = find(Hf-Hi>0,1,'first')-1; % grounding line location 
    Hi(gl:length(xi))=hi(gl:end)*rho_sw/(rho_sw-rho_i); % buoyant thickness using surface
    Hi(Hi>=(hi-hbi))=hi(Hi>=(hi-hbi))-hbi(Hi>=(hi-hbi)); % thickness can't go beneath bed elevation
    
    % add a dummy ice end (hi & Hi)
    for i=c:length(xi)
        hi(i) = hi(i-1)-5; % decrease until reaching 0m  
        Hi(i) = Hi(i-1)-50; % decrease until reaching 0m 
        if Hi(i)>=hi(i)-hbi(i)
            Hi(i)=hi(i)-hbi(i); % can't go beneath bed elevation
        end
    end  
    hi(hi<0)=0; % surface can't go below sea level
    Hi(Hi<0)=0; % no negative thicknesses
    
    % find the end of the ice-covered domain (m along centerline)
    ice_endi = (find(Hi<=0,1,'first')); 
    if isempty(ice_endi) || ice_endi>length(xi)
        ice_endi=length(xi);
    end
    
    % extend variables past where there are NaNs (set to last value)
    Ai(find(isnan(Ai),1,'first'):end) = Ai(find(isnan(Ai),1,'first')-1);
    Ui(find(isnan(Ui),1,'first'):end) = Ui(find(isnan(Ui),1,'first')-1);
    betai(find(isnan(betai),1,'first'):end) = betai(find(isnan(betai),1,'first')-1);
    hbi(find(isnan(hbi),1,'first'):end) = hbi(find(isnan(hbi),1,'first')-1);
    Wi(find(isnan(Wi),1,'first'):end) = Wi(find(isnan(Wi),1,'first')-1);    
    
    % use 90% the observed velocity at the upper bounds
    Ui(1:round(0.2*length(Ui))) = 0.9.*Ui(1:round(0.2*length(Ui)));
    
    % rename initial vectors
    x=xi; h=hi; hb=hbi; W=Wi; H=Hi; A=Ai; beta=betai; U=Ui; dUdx=dUdxi;
    ice_end=ice_endi;
    E = ones(1,length(x)); % enhancement factor
    
%% 1. Run the flowline model for 100 years to reach near steady-state 
%   conditions using 2009 observations
  
save_final = 1; % = 1 to save final glacier geometry

% define time stepping (s)
    dt = 0.001*3.1536e7; 
    t_start = 0*3.1536e7; 
    t_end = 100*3.1536e7;    
    t = (t_start:dt:t_end);
    
% reinitialize vectors
    x=xi; h=hi; hb=hbi; W=Wi; H=Hi; A=Ai; beta=betai; U=Ui; dUdx=dUdxi;
    ice_end=ice_endi; c=ci;    
    
% run flowline model
for i=1:length(t)

    % plot geometry, speed, & calving front position        
    if t(i)==t_start
        col = parula(length(t)+20); %Color scheme for plots
        figure(1); clf
        set(gcf,'Position',[0 100 1400 500]); 
        ax1 = axes('Position',[0.05 0.1 0.28 0.8]); % glacier geometry
            hold on; grid on;
            set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
            title('Glacier Geometry'); legend('Location','northeast'); 
            xlim([0 70]); ylim([min(hb)-100 max(h)+200]);
            xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)');
            % ice surface (1:c)
            plot(x(1:c)./10^3,h(1:c),'color',col(i,:),'linewidth',2,'displayname','0');
            % calving front
            plot(x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
            % ice surface (c:ice_end)
            plot(x(c+1:ice_end)./10^3,h(c+1:ice_end),'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');
            plot(x(c+1:ice_end)./10^3,h(c+1:ice_end)-H(c+1:ice_end),'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');
            % floating bed
            plot(x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off');
            % bed elevation
            plot(x/10^3,hb,'k','linewidth',2,'HandleVisibility','off');
            % mean sea level
            plot([x(1),x(end)]/10^3,[0,0],'k--','HandleVisibility','off');
        ax2 = axes('Position',[0.38 0.1 0.28 0.8]); % ice speed
            hold on; grid on;
            set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
            title('Ice Speed'); legend('Location','northeast');
            xlim([0 70]); ylim([0 4000]); 
            xlabel('Distance Along Centerline (km)'); ylabel('Speed (m yr^{-1})');
            % ice speed (1:c)
            plot(x(1:c)./10^3,U(1:c).*3.1536e7,'color',col(i,:),'linewidth',2,'displayname','0');
            % ice speed (c:ice_end)
            plot(x(c:ice_end)./10^3,U(c:ice_end).*3.1536e7,'--','color',col(i,:),'linewidth',2,'HandleVisibility','off');
        ax3 = axes('Position',[0.7 0.1 0.28 0.8]); % calving front position
            hold on; grid on;
            set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
            title('Calving Front Position'); legend('Location','best');
            xlim([35 70]); ylim([0 10]);
            xlabel('Distance Along Centerline (km)'); ylabel('Year');
            % calving front position
            plot(x(c)/10^3,t(i)./3.1536e7,'.','markersize',15,'color',col(i,:),'displayname','0');
    elseif mod(i-1,round(length(t)/50))==0 % display every length(t)/50  
        figure(1); 
        if mod(i-1,round(length(t)/10))==0 % display in legend every length(t)/10 
            % ice surface (1:c)
            plot(ax1,x(1:c)/10^3,h(1:c),'-','color',col(i,:),'linewidth',2,'displayname',num2str(t(i)./3.1536e7));
            % ice speed (1:c)
            plot(ax2,x(1:c)/10^3,U(1:c).*3.1536e7,'-','Color',col(i,:),'linewidth',2,'DisplayName',num2str(t(i)./3.1536e7));  
            % calving front position
            plot(ax3,x(c)./10^3,t(i)./3.1536e7,'.','Color',col(i,:),'markersize',15,'displayname',num2str(t(i)./3.1536e7)); hold on;                
        else
            % ice surface (1:c)
            plot(ax1,x(1:c)/10^3,h(1:c),'-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
            % ice speed (1:c)
            plot(ax2,x(1:c)/10^3,U(1:c).*3.1536e7,'-','Color',col(i,:),'linewidth',2,'HandleVisibility','off'); hold on;  
            % calving front position
            plot(ax3,x(c)./10^3,t(i)./3.1536e7,'.','Color',col(i,:),'markersize',15,'HandleVisibility','off'); hold on;                
        end
        % ice surface (c+1:ice_end)
        plot(ax1,x(c+1:ice_end)./10^3,h(c+1:ice_end),'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');            
        % calving front
        plot(ax1,x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
        % floating bed (1:c)
        plot(ax1,x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off');
        % floating bed (c+1:ice_end)
        plot(ax1,x(c+1:ice_end)./10^3,h(c+1:ice_end)-H(c+1:ice_end),'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');
        % ice speed (c:ice_end)
        plot(ax2,x(c:ice_end)./10^3,U(c:ice_end).*3.1536e7,'--','color',col(i,:),'linewidth',2,'HandleVisibility','off');
    end

    % calculate the thickness required to remain grounded at each grid cell
    Hf = -(rho_sw./rho_i).*hb; % flotation thickness (m)
    % find the location of the grounding line and use a floating
    % geometry from the grounding line to the calving front
    gl = find(Hf-H>0,1,'first')-1; % grounding line location
    ice_end = find(H<=0,1,'first'); % end of ice-covered domain
    if isempty(ice_end) || ice_end>length(x) || ice_end<c
        ice_end = length(x);
        disp('ice end criteria not met.');
    end

    %calculate the glacier's surface elevation and slope
    h = hb+H; % surface elevation (m a.s.l.)
    h(gl+1:length(x)) = (1-rho_i/rho_sw).*H(gl+1:length(x)); % adjust the surface elevation of ungrounded ice to account for buoyancy
    dhdx = [(h(2:end)-h(1:end-1))./(x(2:end)-x(1:end-1)) 0]; % surface slope (unitless)
    if any(h<0) % surface cannot go below sea level
        H(c:end) = 0;
        h(c:end)=0;
    end

    % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
    Rxx = 2*nthroot((dUdx./(E.*A(1:length(dUdx)))),n)-sigma_b; % resistive stress (Pa)
    crev = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*fwd); % crevasse penetration depth (m)
    c = find(h(1:ice_end-1)-crev(1:ice_end-1)<=0,1,'first'); % calving front located where the inland-most crevasse intersects sea level
    %if the crevasses never intersect sea level
    if isempty(c) == 1
        c = find(H<Hc,1,'first'); %set the calving front to a default minimum ice thickness value
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
        c=ci;
    end

    % calculate the effective pressure (ice overburden pressure minus water
    % pressure) assuming an easy & open connection between the ocean and
    % ice-bed interface
    sl = find(hb<=0,1,'first'); % find where the glacier base first drops below sea level
    N_ground = rho_i*g*H(1:sl); % effective pressure where the bed is above sea level (Pa)
    N_marine = rho_i*g*H(sl+1:length(x))+(rho_sw*g*hb(sl+1:length(x))); % effective pressure where the bed is below sea level (Pa)
    N = [N_ground N_marine];
    N(N<0)=1; % cannot have negative values

    % Solve for new velocity
    [U,~,vm,T] = U_convergence(x,U,U0,dUdx,dhdx,H,A,E,N,W,dx,c,ice_end,n,m,beta,rho_i,rho_sw,g);

    % calculate ice flux
    F = U.*H.*W; % ice flux (m^3 s^-1)
    F(isnan(F))=0;

    % calculate the  change in ice thickness from continuity
    dHdt = -(1./W).*gradient(F,x);
    dH = dHdt.*dt;

    % implement SMB & SMR
    % calculate SMR where ice is ungrounded using a fitting function 
    % past the grounding line to avoid large spatial gradients
    smr = zeros(1,length(x));
        smr(gl+1:ice_end) = feval(fit([x(gl); x(c); x(ice_end)],...
            [0;smr0.*[1;1]],'poly1'),x(gl+1:ice_end));
        smr(ice_end+1:length(x)) = smr(ice_end); 
        smr = movmean(smr,10);
    smb = interp1(x0,smb0+Q0,x); 
        smb = movmean(smb,10);

    % new thickness (change from dynamics, SMB, & SMR)
    Hn = H+dH+(smb.*dt)+(smr.*dt);
    Hn(Hn < 0) = 0; % remove negative values
    H = Hn; %set as the new thickness value
    % smooth out the points near the ice divide
    H(1:10) = ones(1,length(H(1:10))).*nanmean(H(1:10),'all');

    % thickness & surface past calving front
    for j=c+1:length(xi)
        h(j) = h(j-1)-5; % decrease until reaching 0m
        H(j) = H(j-1)-30; % decrease until reaching 0m
        if H(j)>=h(j)-hb(j)
            H(j)=h(j)-hb(j); % can't go beneath bed elevation
        end
    end
    h(h<0)=0; % surface can't go below sea level
    H(H<0)=0; % no negative thicknesses

    % stop the model if it behaves unstably (monitored by ice thickness and speed)
    if max(H) > H_max
        disp(['Adjust dt']);
        break;
    end
    if mean(U) < 200/3.1536e7
        disp('Too slow!');
        break;
    end

    % find the precise location of the grounding line (where H=Hf)
    xf = x(find(Hf-H>0,1,'first')-1);

    %adjust the grid spacing so the grounding line is continuously tracked
    xl = round(xf/dx0); % number of ideal grid spaces needed to reach the grounding line
    dx = xf/xl; % new grid spacing (should be ~dx0)
    xn = 0:dx:L; % new distance vector

    %adjust the space-dependent variables to the new distance vector
    hb = interp1(x,hb,xn); hb(isnan(hb)) = hb(find(~isnan(hb),1,'last')); % glacier bed elevation (m)
    W = interp1(x,W,xn); W(isnan(W)) = W(find(~isnan(W),1,'last')); % glacier width (m)
    H = interp1(x,H,xn,'linear','extrap'); H(isnan(H)) = H(find(~isnan(H),1,'last')); % ice thickness (m)
    Hf = interp1(x,Hf,xn,'linear','extrap'); Hf(isnan(Hf)) = Hf(find(~isnan(Hf),1,'last'));% flotation thickness (m)
    U = interp1(x,U,xn,'linear','extrap'); U(isnan(U)) = U(find(~isnan(U),1,'last')); % speed (m s^-1)
    A = interp1(x,A,x,'linear','extrap'); A(isnan(A)) = A(find(~isnan(A),1,'last'));% rate factor (Pa^-n s^-1)
    beta = interp1(x,beta,xn,'linear','extrap'); beta(isnan(beta)) = beta(find(~isnan(beta),1,'last')); % basal roughness factor

    %find the location of the grounding line and end of the ice-covered domain for the adjusted data
    gl = find(Hf-H<0,1,'last');
    ice_end = find(H<=0,1,'first');
    if isempty(ice_end) || ice_end>length(xn) || ice_end<c
        ice_end = length(xn);
        disp('ice end criteria not met.')
    end

    %rename the distance vector
    x = xn; %distance from the divide (m)

    %calculate the new surface elevation and slope
    h = hb+H; % grounded ice surface elevation (m a.s.l.)
    h(gl+1:length(x)) = (1-rho_i/rho_sw).*H(gl+1:length(x)); % floating ice surface elevation (m a.s.l.)
    dhdx = [(h(2:end)-h(1:end-1))./(x(2:end)-x(1:end-1)) 0]; % surface slope (unitless)
    if any(h<0) % surface cannot go below sea level
        c = find(h<0,1,'first');
        H(c:end) = 0;
        h(c:end)=0;
    end
    H(H>=(h-hb))=h(H>=(h-hb))-hb(H>=(h-hb)); % thickness can't go beneath bed elevation

    % calculate new strain rate
    dUdx = [(U(2:end)-U(1:end-1))./(x(2:end)-x(1:end-1)) 0]; % strain rate

    % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
    Rxx = 2*nthroot((dUdx./(E.*A(1:length(dUdx)))),n)-sigma_b; % resistive stress (Pa)
    crev = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*fwd); % crevasse penetration depth (m)
    c = find(h(1:ice_end-1)-crev(1:ice_end-1)<=0,1,'first'); % calving front located where the inland-most crevasse intersects sea level
    % if the crevasses never intersect sea level
    if isempty(c) == 1
        c = find(H<Hc,1,'first'); %set the calving front to a default minimum ice thickness value
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
        c=ci;
    end
end 

% save conditions
if save_final
    cd([homepath,'inputs-outputs']);
    xj=x; hj=h; hbj=hb; Wj=W; Hj=H; Aj=A; betaj=beta; Uj=U; dUdxj=dUdx;
    ice_endj=ice_end; cj=c;  
    save('Crane_flowline_100yr_output.mat','xj','hj','Wj','Hj','Aj',...
        'betaj','Uj','dUdxj','ice_endj','cj','hbj');
    disp('saved year 100 conditions.');
end
    
%% 2. Tune the inner boundary flux F0 (via U & h)

close all;

save_F0best = 1; % = 1 to save best F0

% define time stepping (s)
    dt = 0.001*3.1536e7; 
    t_start = 0*3.1536e7; 
    t_end = 10*3.1536e7;    
    t = (t_start:dt:t_end);
    
% load 100 yr output conditions
    load('Crane_flowline_100yr_output.mat');

% define inner boundary flux values to test
    F0 = 0:10:50; % (m^3 s^-1)
    
% pre-allocate misfit variables
    U_misfit = NaN.*zeros(length(F0),length(xi));
    dH_tot = NaN.*zeros(length(F0),length(xi)); 
    dH_misfit = NaN.*zeros(length(F0),length(xi));
    
% loop through all inner boundary flux values
for f=1:length(F0)
    
    % initialize variables
    x=xi; h=hi; hb=hbi; W=Wi; H=Hi; A=Ai; beta=betai; U=Ui; dUdx=dUdxi; 
    ice_end=ice_endi; c=ci;
        
    % run flowline model
    for i=1:length(t)

        % plot geometry, speed, & calving front position
        if t(i)==t_start
            col = parula(length(t)+20); %Color scheme for plots
            figure(1); clf
            set(gcf,'Position',[0 100 1400 500]); 
            ax1 = axes('Position',[0.05 0.1 0.28 0.8]); % glacier geometry
                hold on; grid on;
                set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
                title('Glacier Geometry'); legend('Location','northeast'); 
                xlim([0 70]); ylim([min(hb)-100 max(h)+200]);
                xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)');
                % ice surface (1:c)
                plot(x(1:c)./10^3,h(1:c),'color',col(i,:),'linewidth',2,'displayname','0');
                % calving front
                plot(x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
                % ice surface (c:ice_end)
                plot(x(c+1:ice_end)./10^3,h(c+1:ice_end),'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');
                plot(x(c+1:ice_end)./10^3,h(c+1:ice_end)-H(c+1:ice_end),'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');
                % floating bed
                plot(x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off');
                % bed elevation
                plot(x/10^3,hb,'k','linewidth',2,'HandleVisibility','off');
                % mean sea level
                plot([x(1),x(end)]/10^3,[0,0],'k--','HandleVisibility','off');
            ax2 = axes('Position',[0.38 0.1 0.28 0.8]); % ice speed
                hold on; grid on;
                set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
                title('Ice Speed'); legend('Location','northeast');
                xlim([0 70]); ylim([0 4000]); 
                xlabel('Distance Along Centerline (km)'); ylabel('Speed (m yr^{-1})');
                % ice speed (1:c)
                plot(x(1:c)./10^3,U(1:c).*3.1536e7,'color',col(i,:),'linewidth',2,'displayname','0');
                % ice speed (c:ice_end)
                plot(x(c:ice_end)./10^3,U(c:ice_end).*3.1536e7,'--','color',col(i,:),'linewidth',2,'HandleVisibility','off');
            ax3 = axes('Position',[0.7 0.1 0.28 0.8]); % calving front position
                hold on; grid on;
                set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
                title('Calving Front Position'); legend('Location','best');
                xlim([35 70]); ylim([0 10]);
                xlabel('Distance Along Centerline (km)'); ylabel('Year');
                % calving front position
                plot(x(c)/10^3,t(i)./3.1536e7,'.','markersize',15,'color',col(i,:),'displayname','0');
        elseif mod(i-1,round(length(t)/50))==0 % display every length(t)/50  
            figure(1); 
            if mod(i-1,round(length(t)/10))==0 % display in legend every length(t)/10 
                % ice surface (1:c)
                plot(ax1,x(1:c)/10^3,h(1:c),'-','color',col(i,:),'linewidth',2,'displayname',num2str(t(i)./3.1536e7));
                % ice speed (1:c)
                plot(ax2,x(1:c)/10^3,U(1:c).*3.1536e7,'-','Color',col(i,:),'linewidth',2,'DisplayName',num2str(t(i)./3.1536e7));  
                % calving front position
                plot(ax3,x(c)./10^3,t(i)./3.1536e7,'.','Color',col(i,:),'markersize',15,'displayname',num2str(t(i)./3.1536e7)); hold on;                
            else
                % ice surface (1:c)
                plot(ax1,x(1:c)/10^3,h(1:c),'-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
                % ice speed (1:c)
                plot(ax2,x(1:c)/10^3,U(1:c).*3.1536e7,'-','Color',col(i,:),'linewidth',2,'HandleVisibility','off'); hold on;  
                % calving front position
                plot(ax3,x(c)./10^3,t(i)./3.1536e7,'.','Color',col(i,:),'markersize',15,'HandleVisibility','off'); hold on;                
            end
            % ice surface (c+1:ice_end)
            plot(ax1,x(c+1:ice_end)./10^3,h(c+1:ice_end),'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');            
            % calving front
            plot(ax1,x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
            % floating bed (1:c)
            plot(ax1,x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off');
            % floating bed (c+1:ice_end)
            plot(ax1,x(c+1:ice_end)./10^3,h(c+1:ice_end)-H(c+1:ice_end),'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');
            % ice speed (c:ice_end)
            plot(ax2,x(c:ice_end)./10^3,U(c:ice_end).*3.1536e7,'--','color',col(i,:),'linewidth',2,'HandleVisibility','off');
        end

        % calculate the thickness required to remain grounded at each grid cell
        Hf = -(rho_sw./rho_i).*hb; % flotation thickness (m)
        % find the location of the grounding line and use a floating
        % geometry from the grounding line to the calving front
        gl = find(Hf-H>0,1,'first')-1; % grounding line location
        ice_end = find(H<=50,1,'first'); % end of ice-covered domain
        if isempty(ice_end) || ice_end>length(x) || ice_end<c
            ice_end = length(x);
            disp('ice end criteria not met.');
        end

        %calculate the glacier's surface elevation and slope
        h = hb+H; % surface elevation (m a.s.l.)
        h(gl+1:length(x)) = (1-rho_i/rho_sw).*H(gl+1:length(x)); % adjust the surface elevation of ungrounded ice to account for buoyancy
        dhdx = [(h(2:end)-h(1:end-1))./(x(2:end)-x(1:end-1)) 0]; % surface slope (unitless)
        if any(h<0) % surface cannot go below sea level
            H(c:end) = 0;
            h(c:end)=0;
        end

        % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
        Rxx = 2*nthroot((dUdx./(E.*A(1:length(dUdx)))),n)-sigma_b; % resistive stress (Pa)
        crev = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*fwd); % crevasse penetration depth (m)
        c = find(h(1:ice_end-1)-crev(1:ice_end-1)<=0,1,'first'); % calving front located where the inland-most crevasse intersects sea level
        %if the crevasses never intersect sea level
        if isempty(c) == 1
            c = find(H<Hc,1,'first'); %set the calving front to a default minimum ice thickness value
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
            c=ci;
        end

        % calculate the effective pressure (ice overburden pressure minus water
        % pressure) assuming an easy & open connection between the ocean and
        % ice-bed interface
        sl = find(hb<=0,1,'first'); % find where the glacier base first drops below sea level
        N_ground = rho_i*g*H(1:sl); % effective pressure where the bed is above sea level (Pa)
        N_marine = rho_i*g*H(sl+1:length(x))+(rho_sw*g*hb(sl+1:length(x))); % effective pressure where the bed is below sea level (Pa)
        N = [N_ground N_marine];
        N(N<0)=1; % cannot have negative values

        % Solve for new velocity
        [U,~,vm,T] = U_convergence(x,U,U0,dUdx,dhdx,H,A,E,N,W,dx,c,ice_end,n,m,beta,rho_i,rho_sw,g);

        % calculate ice flux
        F = U.*H.*W; % ice flux (m^3 s^-1)
        F(isnan(F))=0;
        F(1) = F0(f); 

        % calculate the  change in ice thickness from continuity
        dHdt = -(1./W).*gradient(F,x);
        dH = dHdt.*dt;

        % implement SMB & SMR
        % calculate SMR where ice is ungrounded using a fitting function 
        % past the grounding line to avoid large spatial gradients
        smr = zeros(1,length(x));
            smr(gl+1:ice_end) = feval(fit([x(gl); x(c); x(ice_end)],...
                [0;smr0.*[1;1]],'poly1'),x(gl+1:ice_end));
            smr(ice_end+1:length(x)) = smr(ice_end); 
            smr = movmean(smr,10);
        smb = interp1(x0,smb0+Q0,x); 
            smb = movmean(smb,10);

        % new thickness (change from dynamics, SMB, & SMR)
        Hn = H+dH+(smb.*dt)+(smr.*dt);
        Hn(Hn < 0) = 0; % remove negative values
        H = Hn; %set as the new thickness value
        % smooth out the points near the ice divide
        H(1:10) = ones(1,length(H(1:10))).*nanmean(H(1:10),'all');

        % thickness & surface near ice_end
        for j=ice_end-10:length(xi)
            h(j) = h(j-1)-5; % decrease until reaching 0m
            H(j) = H(j-1)-30; % decrease until reaching 0m
            if H(j)>=h(j)-hb(j)
                H(j)=h(j)-hb(j); % can't go beneath bed elevation
            end
        end
        h(h<0)=0; % surface can't go below sea level
        H(H<0)=0; % no negative thicknesses

        % stop the model if it behaves unstably (monitored by ice thickness and speed)
        if max(H) > H_max
            disp(['Adjust dt']);
            break;
        end
        if mean(U) < 200/3.1536e7
            disp('Too slow!');
            break;
        end

        % find the precise location of the grounding line (where H=Hf)
        xf = x(find(Hf-H>0,1,'first')-1);

        %adjust the grid spacing so the grounding line is continuously tracked
        xl = round(xf/dx0); % number of ideal grid spaces needed to reach the grounding line
        dx = xf/xl; % new grid spacing (should be ~dx0)
        xn = 0:dx:L; % new distance vector

        %adjust the space-dependent variables to the new distance vector
        hb = interp1(x,hb,xn); hb(isnan(hb)) = hb(find(~isnan(hb),1,'last')); % glacier bed elevation (m)
        W = interp1(x,W,xn); W(isnan(W)) = W(find(~isnan(W),1,'last')); % glacier width (m)
        H = interp1(x,H,xn,'linear','extrap'); H(isnan(H)) = H(find(~isnan(H),1,'last')); % ice thickness (m)
        Hf = interp1(x,Hf,xn,'linear','extrap'); Hf(isnan(Hf)) = Hf(find(~isnan(Hf),1,'last'));% flotation thickness (m)
        U = interp1(x,U,xn,'linear','extrap'); U(isnan(U)) = U(find(~isnan(U),1,'last')); % speed (m s^-1)
        A = interp1(x,A,x,'linear','extrap'); A(isnan(A)) = A(find(~isnan(A),1,'last'));% rate factor (Pa^-n s^-1)
        beta = interp1(x,beta,xn,'linear','extrap'); beta(isnan(beta)) = beta(find(~isnan(beta),1,'last')); % basal roughness factor

        %find the location of the grounding line and end of the ice-covered domain for the adjusted data
        gl = find(Hf-H<0,1,'last');
        ice_end = find(H<=50,1,'first');
        if isempty(ice_end) || ice_end>length(xn) || ice_end<c
            ice_end = length(xn);
            disp('ice end criteria not met.')
        end

        %rename the distance vector
        x = xn; %distance from the divide (m)

        %calculate the new surface elevation and slope
        h = hb+H; % grounded ice surface elevation (m a.s.l.)
        h(gl+1:length(x)) = (1-rho_i/rho_sw).*H(gl+1:length(x)); % floating ice surface elevation (m a.s.l.)
        dhdx = [(h(2:end)-h(1:end-1))./(x(2:end)-x(1:end-1)) 0]; % surface slope (unitless)
        if any(h<0) % surface cannot go below sea level
            c = find(h<0,1,'first');
            H(c:end) = 0;
            h(c:end)=0;
        end
        H(H>=(h-hb))=h(H>=(h-hb))-hb(H>=(h-hb)); % thickness can't go beneath bed elevation

        % calculate new strain rate
        dUdx = [(U(2:end)-U(1:end-1))./(x(2:end)-x(1:end-1)) 0]; % strain rate

        % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
        Rxx = 2*nthroot((dUdx./(E.*A(1:length(dUdx)))),n)-sigma_b; % resistive stress (Pa)
        crev = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*fwd); % crevasse penetration depth (m)
        c = find(h(1:ice_end-1)-crev(1:ice_end-1)<=0,1,'first'); % calving front located where the inland-most crevasse intersects sea level
        % if the crevasses never intersect sea level
        if isempty(c) == 1
            c = find(H<Hc,1,'first'); %set the calving front to a default minimum ice thickness value
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
            c=ci;
        end

        % calculate the misfit in surface speed and elevation change at the
        % final model time
        if t(i)==t_end
            U_misfit(f,:) = interp1(x0,interp1(x,U,x0)-U_obs(9).U',xi);
            dH_tot(f,:) = interp1(x,h,xi)-hi; 
            dH_misfit(f,:) = dH_tot(f,:)-dH_obs;
        end
    end
    
end

% calculate mean misfit
misfit = nanmean(U_misfit,2).*nanmean(dH_misfit,2);

% plot results
figure(10); clf
set(gca,'fontsize',14,'linewidth',2); grid on; hold on;
xlabel('F0 (m^3/s)'); ylabel('misfit'); 
for i=1:length(F0)
    if any(isnan(U_misfit(i,:)))
        plot(F0(i),misfit(i),'*r','markersize',15,'linewidth',2);
    else
        plot(F0(i),misfit(i),'*b','markersize',15,'linewidth',2);
    end
end

% plot results for optimal F0
F0best = F0(abs(misfit)==min(abs(misfit))); % m^3/s
disp(['Best F0 = ',num2str(F0best),' m^3/s']);

% save F0best
if save_F0best
    cd([homepath,'inputs-outputs']);
    save('F0best.mat','F0best');
    disp('Best F0 saved.');
end

%% 3. Tune for the back pressure sigma_b (via c) 

% Rxx = 2(dUdx/(EA))^(1/n) - sigma_b
%   where sigma_b = back pressure 

close all; 

save_sigma_bbest = 1; % = 1 to save sigma_bbest

% define time stepping (s)
    dt = 0.001*3.1536e7; 
    t_start = 0*3.1536e7; 
    t_end = 10*3.1536e7;    
    t = (t_start:dt:t_end);
    
% load 100 yr output conditions and optimal parameters
    load('Crane_flowline_100yr_output.mat');
    F0=load('F0best.mat').F0best;
    
% define sigma_b values to test
    sigma_b = 0:20e3:60e3; % back pressure (Pa)

% pre-allocate misfit variables
    c_misfit = NaN.*zeros(length(sigma_b),length(2009:2018)); 
    
% loop through all sigma_b values
for f=1:length(sigma_b)
    
    % initialize variables
    x=xi; h=hi; hb=hbi; W=Wi; H=Hi; A=Ai; beta=betai; U=Ui; dUdx=dUdxi; 
    ice_end=ice_endi; c=ci;
        
    % run flowline model
    for i=1:length(t)

        % plot geometry, speed, & calving front position        
        if t(i)==t_start
            col = parula(length(t)+20); %Color scheme for plots
            figure(1); clf
            set(gcf,'Position',[0 100 1400 500]); 
            ax1 = axes('Position',[0.05 0.1 0.28 0.8]); % glacier geometry
                hold on; grid on;
                set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
                title('Glacier Geometry'); legend('Location','northeast'); 
                xlim([0 70]); ylim([min(hb)-100 max(h)+200]);
                xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)');
                % ice surface (1:c)
                plot(x(1:c)./10^3,h(1:c),'color',col(i,:),'linewidth',2,'displayname','0');
                % calving front
                plot(x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
                % ice surface (c:ice_end)
                plot(x(c+1:ice_end)./10^3,h(c+1:ice_end),'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');
                plot(x(c+1:ice_end)./10^3,h(c+1:ice_end)-H(c+1:ice_end),'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');
                % floating bed
                plot(x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off');
                % bed elevation
                plot(x/10^3,hb,'k','linewidth',2,'HandleVisibility','off');
                % mean sea level
                plot([x(1),x(end)]/10^3,[0,0],'k--','HandleVisibility','off');
            ax2 = axes('Position',[0.38 0.1 0.28 0.8]); % ice speed
                hold on; grid on;
                set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
                title('Ice Speed'); legend('Location','northeast');
                xlim([0 70]); ylim([0 4000]); 
                xlabel('Distance Along Centerline (km)'); ylabel('Speed (m yr^{-1})');
                % ice speed (1:c)
                plot(x(1:c)./10^3,U(1:c).*3.1536e7,'color',col(i,:),'linewidth',2,'displayname','0');
                % ice speed (c:ice_end)
                plot(x(c:ice_end)./10^3,U(c:ice_end).*3.1536e7,'--','color',col(i,:),'linewidth',2,'HandleVisibility','off');
            ax3 = axes('Position',[0.7 0.1 0.28 0.8]); % calving front position
                hold on; grid on;
                set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
                title('Calving Front Position'); legend('Location','best');
                xlim([35 70]); ylim([0 10]);
                xlabel('Distance Along Centerline (km)'); ylabel('Year');
                % calving front position
                plot(x(c)/10^3,t(i)./3.1536e7,'.','markersize',15,'color',col(i,:),'displayname','0');
        elseif mod(i-1,round(length(t)/50))==0 % display every length(t)/50  
            figure(1); 
            if mod(i-1,round(length(t)/10))==0 % display in legend every length(t)/10 
                % ice surface (1:c)
                plot(ax1,x(1:c)/10^3,h(1:c),'-','color',col(i,:),'linewidth',2,'displayname',num2str(t(i)./3.1536e7));
                % ice speed (1:c)
                plot(ax2,x(1:c)/10^3,U(1:c).*3.1536e7,'-','Color',col(i,:),'linewidth',2,'DisplayName',num2str(t(i)./3.1536e7));  
                % calving front position
                plot(ax3,x(c)./10^3,t(i)./3.1536e7,'.','Color',col(i,:),'markersize',15,'displayname',num2str(t(i)./3.1536e7)); hold on;                
            else
                % ice surface (1:c)
                plot(ax1,x(1:c)/10^3,h(1:c),'-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
                % ice speed (1:c)
                plot(ax2,x(1:c)/10^3,U(1:c).*3.1536e7,'-','Color',col(i,:),'linewidth',2,'HandleVisibility','off'); hold on;  
                % calving front position
                plot(ax3,x(c)./10^3,t(i)./3.1536e7,'.','Color',col(i,:),'markersize',15,'HandleVisibility','off'); hold on;                
            end
            % ice surface (c+1:ice_end)
            plot(ax1,x(c+1:ice_end)./10^3,h(c+1:ice_end),'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');            
            % calving front
            plot(ax1,x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
            % floating bed (1:c)
            plot(ax1,x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off');
            % floating bed (c+1:ice_end)
            plot(ax1,x(c+1:ice_end)./10^3,h(c+1:ice_end)-H(c+1:ice_end),'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');
            % ice speed (c:ice_end)
            plot(ax2,x(c:ice_end)./10^3,U(c:ice_end).*3.1536e7,'--','color',col(i,:),'linewidth',2,'HandleVisibility','off');
        end

        % calculate the thickness required to remain grounded at each grid cell
        Hf = -(rho_sw./rho_i).*hb; %flotation thickness (m)
        % find the location of the grounding line and use a floating
        % geometry from the grounding line to the calving front
        gl = find(Hf-H>0,1,'first')-1; %grounding line location
        ice_end = find(H<=0,1,'first'); %end of ice-covered domain
        if isempty(ice_end) || ice_end>length(x) || ice_end<c
            ice_end = length(x);
            disp('ice end criteria not met.');
        end

        %calculate the glacier's surface elevation and slope
        h = hb+H; %h = surface elevation (m a.s.l.)
        h(gl+1:length(x)) = (1-rho_i/rho_sw).*H(gl+1:length(x)); %adjust the surface elevation of ungrounded ice to account for buoyancy
        dhdx = [(h(2:end)-h(1:end-1))./(x(2:end)-x(1:end-1)) 0]; % surface slope (unitless)
        if any(h<0) % surface cannot go below sea level
            H(c:end) = 0;
            h(c:end)=0;
        end

        % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
        Rxx = 2*nthroot((dUdx./(E.*A(1:length(dUdx)))),n)-sigma_b(f); % resistive stress (Pa)
        crev = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*fwd); % crevasse penetration depth (m)
        c = find(h(1:ice_end-1)-crev(1:ice_end-1)<=0,1,'first'); % calving front located where the inland-most crevasse intersects sea level
        %if the crevasses never intersect sea level
        if isempty(c) == 1
            c = find(H<Hc,1,'first'); %set the calving front to a default minimum ice thickness value
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
            c=ci;
        end

        % calculate the effective pressure (ice overburden pressure minus water
        % pressure) assuming an easy & open connection between the ocean and
        % ice-bed interface
        sl = find(hb<=0,1,'first'); % find where the glacier base first drops below sea level
        N_ground = rho_i*g*H(1:sl); % effective pressure where the bed is above sea level (Pa)
        N_marine = rho_i*g*H(sl+1:length(x))+(rho_sw*g*hb(sl+1:length(x))); % effective pressure where the bed is below sea level (Pa)
        N = [N_ground N_marine];
        N(N<0)=1; % cannot have negative values

        % Solve for new velocity
        [U,~,vm,T] = U_convergence(x,U,U0,dUdx,dhdx,H,A,E,N,W,dx,c,ice_end,n,m,beta,rho_i,rho_sw,g);

        % calculate ice flux
        F = U.*H.*W; % ice flux (m^3 s^-1)
        F(isnan(F))=0;
        F(1) = F0; 

        % calculate the  change in ice thickness from continuity
        dHdt = -(1./W).*gradient(F,x);
        dH = dHdt.*dt;

        % implement SMB & SMR
        % calculate SMR where ice is ungrounded and fit a smoothing spline 
        % past the grounding line to avoid large spatial gradients
        smr = zeros(1,length(x));
            smr(gl+1:ice_end) = feval(fit([x(gl); x(c); x(ice_end)],...
                [0;smr0.*[1;1]],'poly1'),x(gl+1:ice_end));
            smr(ice_end+1:length(x)) = smr(ice_end); 
            smr = movmean(smr,10);
        smb = interp1(x0,smb0+Q0,x); 
            smb = movmean(smb,10);

        % new thickness (change from dynamics, SMB, & SMR)
        Hn = H+dH+(smb.*dt)+(smr.*dt);
        Hn(Hn < 0) = 0; % remove negative values
        H = Hn; % set as the new thickness value
        % smooth out the points near the ice divide
        H(1:10) = ones(1,length(H(1:10))).*nanmean(H(1:10),'all');

        % thickness & surface past calving front
        for j=c+1:length(xi)
            h(j) = h(j-1)-5; % decrease until reaching 0m
            H(j) = H(j-1)-30; % decrease until reaching 0m
            if H(j)>=h(j)-hb(j)
                H(j)=h(j)-hb(j); % can't go beneath bed elevation
            end
        end
        h(h<0)=0; % surface can't go below sea level
        H(H<0)=0; % no negative thicknesses

        % stop the model if it behaves unstably (monitored by ice thickness and speed)
        if max(H) > H_max
            disp(['Adjust dt']);
            break;
        end
        if mean(U) < 200/3.1536e7
            disp('Too slow!');
            break;
        end

        % find the precise location of the grounding line (where H=Hf)
        xf = x(find(Hf-H>0,1,'first')-1);

        %adjust the grid spacing so the grounding line is continuously tracked
        xl = round(xf/dx0); % number of ideal grid spaces needed to reach the grounding line
        dx = xf/xl; % new grid spacing (should be ~dx0)
        xn = 0:dx:L; % new distance vector

        %adjust the space-dependent variables to the new distance vector
        hb = interp1(x,hb,xn); hb(isnan(hb)) = hb(find(~isnan(hb),1,'last')); % glacier bed elevation (m)
        W = interp1(x,W,xn); W(isnan(W)) = W(find(~isnan(W),1,'last')); % glacier width (m)
        H = interp1(x,H,xn,'linear','extrap'); H(isnan(H)) = H(find(~isnan(H),1,'last')); % ice thickness (m)
        Hf = interp1(x,Hf,xn,'linear','extrap'); Hf(isnan(Hf)) = Hf(find(~isnan(Hf),1,'last'));% flotation thickness (m)
        U = interp1(x,U,xn,'linear','extrap'); U(isnan(U)) = U(find(~isnan(U),1,'last')); % speed (m s^-1)
        A = interp1(x,A,x,'linear','extrap'); A(isnan(A)) = A(find(~isnan(A),1,'last'));% rate factor (Pa^-n s^-1)
        beta = interp1(x,beta,xn,'linear','extrap'); beta(isnan(beta)) = beta(find(~isnan(beta),1,'last')); % basal roughness factor

        %find the location of the grounding line and end of the ice-covered domain for the adjusted data
        gl = find(Hf-H<0,1,'last');
        ice_end = find(H<=0,1,'first');
        if isempty(ice_end) || ice_end>length(xn) || ice_end<c
            ice_end = length(xn);
            disp('ice end criteria not met.')
        end

        %rename the distance vector
        x = xn; %distance from the divide (m)

        %calculate the new surface elevation and slope
        h = hb+H; % grounded ice surface elevation (m a.s.l.)
        h(gl+1:length(x)) = (1-rho_i/rho_sw).*H(gl+1:length(x)); % floating ice surface elevation (m a.s.l.)
        dhdx = [(h(2:end)-h(1:end-1))./(x(2:end)-x(1:end-1)) 0]; % surface slope (unitless)
        if any(h<0) % surface cannot go below sea level
            c = find(h<0,1,'first');
            H(c:end) = 0;
            h(c:end)=0;
        end
        H(H>=(h-hb))=h(H>=(h-hb))-hb(H>=(h-hb)); % thickness can't go beneath bed elevation

        % calculate new strain rate
        dUdx = [(U(2:end)-U(1:end-1))./(x(2:end)-x(1:end-1)) 0]; % strain rate

        % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
        Rxx = 2*nthroot((dUdx./(E.*A(1:length(dUdx)))),n)-sigma_b(f); % resistive stress (Pa)
        crev = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*fwd); % crevasse penetration depth (m)
        c = find(h(1:ice_end-1)-crev(1:ice_end-1)<=0,1,'first'); % calving front located where the inland-most crevasse intersects sea level
        % if the crevasses never intersect sea level
        if isempty(c) == 1
            c = find(H<Hc,1,'first'); %set the calving front to a default minimum ice thickness value
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
            c=ci;
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
figure(10); clf
    plot(sigma_b./10^3,misfit,'*b','markersize',15,'linewidth',2);
    set(gca,'fontsize',14,'fontname','Arial','linewidth',2);
    grid on; xlabel('\sigma_b (kPa)'); ylabel('misfit');

% save sigma_bbest
if save_sigma_bbest
    cd([homepath,'inputs-outputs']);
    sigma_bbest = sigma_b(abs(misfit)==min(abs(misfit)));
    save('sigma_bbest.mat','sigma_bbest');
    disp(['Best sigma_b = ',num2str(sigma_bbest./10^3),' kPa saved.']);
end

%% 4. Tune enhancement factor E (via U & h)

close all; 

save_Ebest = 1; % = 1 to save Ebest

% define time stepping (s)
    dt = 0.001*3.1536e7; 
    t_start = 0*3.1536e7; 
    t_end = 10*3.1536e7;    
    t = (t_start:dt:t_end);
    
% load 100 yr output conditions
    load('Crane_flowline_100yr_output.mat');
    F0 = load('F0best.mat').F0best;
    sigma_b = load('sigma_bbest.mat').sigma_bbest;
    
% define Efit values to tests 
    clear Efit
    ub_E = [0.5 2]; % upper bound (minimum E near ice divide)
    lb_E = [0.5 2]; % lower bound (maximum E at terminus)
    deg_E = 0:2; % polynomial degrees
    col_E = winter(length(deg_E)*length(lb_E)*length(ub_E)); % color scheme for plotting
    % set up figure
    figure(11); clf
    set(gca,'fontsize',14,'linewidth',2);
    xlabel('distance along centerline (km)'); ylabel('enhancement factor');
    grid on; hold on; title('Enhancement Factor Fits');
    count = 0; % counter
    % loop through upper bounds
    for i=1:length(deg_E)
        % loop through lower bounds
        for j=1:length(ub_E)
            % loop through polynomial degrees
            for k=1:length(lb_E)
                count=count+1;
                Efit(count).degree = deg_E(i);
                Efit(count).fit = polyfit([x(1) x(end)],[ub_E(j) lb_E(k)],deg_E(i));
                Efit(count).E = polyval(Efit(count).fit,x);
                figure(10); hold on;
                plot(x./10^3,Efit(count).E,'color',col_E(count,:),'linewidth',2);
            end
        end
    end
    clear ub_E lb_E deg_E

% pre-allocate dH variables
    dH_tot = NaN.*zeros(length(Efit),length(x)); % total dH on each iteration
    dH_misfit = NaN.*zeros(length(Efit),length(x)); % misfit of total dH 
        
% loop through all Efit values
for f=1:length(Efit)
    
    % define E
    E = Efit(f).E; 
    
    % initialize variables
    x=xi; h=hi; hb=hbi; W=Wi; H=Hi; A=Ai; beta=betai; U=Ui; dUdx=dUdxi; 
    ice_end=ice_endi; c=ci;
        
    % run flowline model
    for i=1:length(t)
    
        % plot geometry, speed, & calving front position
        if t(i)==t_start
            col = parula(length(t)+20); %Color scheme for plots
            figure(1); clf
            set(gcf,'Position',[0 100 1400 400]); 
            ax1 = axes('Position',[0.05 0.1 0.28 0.8]); % glacier geometry
                hold on; grid on;
                set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
                title('Glacier Geometry'); legend('Location','northeast'); 
                xlim([0 70]); ylim([min(hb)-100 max(h)+200]);
                xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)');
                % ice surface (1:c)
                plot(x(1:c)./10^3,h(1:c),'color',col(i,:),'linewidth',2,'displayname','0');
                % calving front
                plot(x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
                % ice surface (c:ice_end)
                plot(x(c+1:ice_end)./10^3,h(c+1:ice_end),'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');
                plot(x(c+1:ice_end)./10^3,h(c+1:ice_end)-H(c+1:ice_end),'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');
                % floating bed
                plot(x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off');
                % bed elevation
                plot(x/10^3,hb,'k','linewidth',2,'HandleVisibility','off');
                % mean sea level
                plot([x(1),x(end)]/10^3,[0,0],'k--','HandleVisibility','off');
            ax2 = axes('Position',[0.38 0.1 0.28 0.8]); % ice speed
                hold on; grid on;
                set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
                title('Ice Speed'); legend('Location','northeast');
                xlim([0 70]); ylim([0 4000]); 
                xlabel('Distance Along Centerline (km)'); ylabel('Speed (m yr^{-1})');
                % ice speed (1:c)
                plot(x(1:c)./10^3,U(1:c).*3.1536e7,'color',col(i,:),'linewidth',2,'displayname','0');
                % ice speed (c:ice_end)
                plot(x(c:ice_end)./10^3,U(c:ice_end).*3.1536e7,'--','color',col(i,:),'linewidth',2,'HandleVisibility','off');
            ax3 = axes('Position',[0.7 0.1 0.28 0.8]); % calving front position
                hold on; grid on;
                set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
                title('Calving Front Position'); legend('Location','best');
                xlim([35 70]); ylim([0 10]);
                xlabel('Distance Along Centerline (km)'); ylabel('Year');
                % calving front position
                plot(x(c)/10^3,t(i)./3.1536e7,'.','markersize',15,'color',col(i,:),'displayname','0');
        elseif mod(i-1,round(length(t)/50))==0 % display every length(t)/50  
            figure(1); 
            if mod(i-1,round(length(t)/10))==0 % display in legend every length(t)/10 
                % ice surface (1:c)
                plot(ax1,x(1:c)/10^3,h(1:c),'-','color',col(i,:),'linewidth',2,'displayname',num2str(t(i)./3.1536e7));
                % ice speed (1:c)
                plot(ax2,x(1:c)/10^3,U(1:c).*3.1536e7,'-','Color',col(i,:),'linewidth',2,'DisplayName',num2str(t(i)./3.1536e7));  
                % calving front position
                plot(ax3,x(c)./10^3,t(i)./3.1536e7,'.','Color',col(i,:),'markersize',15,'displayname',num2str(t(i)./3.1536e7)); hold on;                
            else
                % ice surface (1:c)
                plot(ax1,x(1:c)/10^3,h(1:c),'-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
                % ice speed (1:c)
                plot(ax2,x(1:c)/10^3,U(1:c).*3.1536e7,'-','Color',col(i,:),'linewidth',2,'HandleVisibility','off'); hold on;  
                % calving front position
                plot(ax3,x(c)./10^3,t(i)./3.1536e7,'.','Color',col(i,:),'markersize',15,'HandleVisibility','off'); hold on;                
            end
            % ice surface (c+1:ice_end)
            plot(ax1,x(c+1:ice_end)./10^3,h(c+1:ice_end),'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');            
            % calving front
            plot(ax1,x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
            % floating bed (1:c)
            plot(ax1,x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off');
            % floating bed (c+1:ice_end)
            plot(ax1,x(c+1:ice_end)./10^3,h(c+1:ice_end)-H(c+1:ice_end),'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');
            % ice speed (c:ice_end)
            plot(ax2,x(c:ice_end)./10^3,U(c:ice_end).*3.1536e7,'--','color',col(i,:),'linewidth',2,'HandleVisibility','off');
        end

        % calculate the thickness required to remain grounded at each grid cell
        Hf = -(rho_sw./rho_i).*hb; % flotation thickness (m)
        % find the location of the grounding line and use a floating
        % geometry from the grounding line to the calving front
        gl = find(Hf-H>0,1,'first')-1; % grounding line location
        ice_end = find(H<=50,1,'first'); % end of ice-covered domain
        if isempty(ice_end) || ice_end>length(x) || ice_end<c
            ice_end = length(x);
            disp('ice end criteria not met.');
        end

        %calculate the glacier's surface elevation and slope
        h = hb+H; % surface elevation (m a.s.l.)
        h(gl+1:length(x)) = (1-rho_i/rho_sw).*H(gl+1:length(x)); % adjust the surface elevation of ungrounded ice to account for buoyancy
        dhdx = [(h(2:end)-h(1:end-1))./(x(2:end)-x(1:end-1)) 0]; % surface slope (unitless)
        if any(h<0) % surface cannot go below sea level
            H(c:end) = 0;
            h(c:end)=0;
        end

        % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
        Rxx = 2*nthroot((dUdx./(E.*A(1:length(dUdx)))),n)-sigma_b; % resistive stress (Pa)
        crev = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*fwd); % crevasse penetration depth (m)
        c = find(h(1:ice_end-1)-crev(1:ice_end-1)<=0,1,'first'); % calving front located where the inland-most crevasse intersects sea level
        %if the crevasses never intersect sea level
        if isempty(c) == 1
            c = find(H<Hc,1,'first'); %set the calving front to a default minimum ice thickness value
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
            c=ci;
        end

        % calculate the effective pressure (ice overburden pressure minus water
        % pressure) assuming an easy & open connection between the ocean and
        % ice-bed interface
        sl = find(hb<=0,1,'first'); % find where the glacier base first drops below sea level
        N_ground = rho_i*g*H(1:sl); % effective pressure where the bed is above sea level (Pa)
        N_marine = rho_i*g*H(sl+1:length(x))+(rho_sw*g*hb(sl+1:length(x))); % effective pressure where the bed is below sea level (Pa)
        N = [N_ground N_marine];
        N(N<0)=1; % cannot have negative values

        % Solve for new velocity
        [U,~,vm,T] = U_convergence(x,U,U0,dUdx,dhdx,H,A,E,N,W,dx,c,ice_end,n,m,beta,rho_i,rho_sw,g);

        % calculate ice flux
        F = U.*H.*W; % ice flux (m^3 s^-1)
        F(isnan(F))=0;
        F(1) = F0; 

        % calculate the  change in ice thickness from continuity
        dHdt = -(1./W).*gradient(F,x);
        dH = dHdt.*dt;

        % implement SMB & SMR
        % calculate SMR where ice is ungrounded using a fitting function 
        % past the grounding line to avoid large spatial gradients
        smr = zeros(1,length(x));
            smr(gl+1:ice_end) = feval(fit([x(gl); x(c); x(ice_end)],...
                [0;smr0.*[1;1]],'poly1'),x(gl+1:ice_end));
            smr(ice_end+1:length(x)) = smr(ice_end); 
            smr = movmean(smr,10);
        smb = interp1(x0,smb0+Q0,x); 
            smb = movmean(smb,10);

        % new thickness (change from dynamics, SMB, & SMR)
        Hn = H+dH+(smb.*dt)+(smr.*dt);
        Hn(Hn < 0) = 0; % remove negative values
        H = Hn; %set as the new thickness value
        % smooth out the points near the ice divide
        H(1:10) = ones(1,length(H(1:10))).*nanmean(H(1:10),'all');

        % thickness & surface past calving front
        for j=c+1:length(xi)
            h(j) = h(j-1)-5; % decrease until reaching 0m
            H(j) = H(j-1)-50; % decrease until reaching 0m
            if H(j)>=h(j)-hb(j)
                H(j)=h(j)-hb(j); % can't go beneath bed elevation
            end
        end
        h(h<0)=0; % surface can't go below sea level
        H(H<0)=0; % no negative thicknesses

        % stop the model if it behaves unstably (monitored by ice thickness and speed)
        if max(H) > H_max
            disp(['Adjust dt']);
            break;
        end
        if mean(U) < 200/3.1536e7
            disp('Too slow!');
            break;
        end

        % find the precise location of the grounding line (where H=Hf)
        xf = x(find(Hf-H>0,1,'first')-1);

        %adjust the grid spacing so the grounding line is continuously tracked
        xl = round(xf/dx0); % number of ideal grid spaces needed to reach the grounding line
        dx = xf/xl; % new grid spacing (should be ~dx0)
        xn = 0:dx:L; % new distance vector

        %adjust the space-dependent variables to the new distance vector
        hb = interp1(x,hb,xn); hb(isnan(hb)) = hb(find(~isnan(hb),1,'last')); % glacier bed elevation (m)
        W = interp1(x,W,xn); W(isnan(W)) = W(find(~isnan(W),1,'last')); % glacier width (m)
        H = interp1(x,H,xn,'linear','extrap'); H(isnan(H)) = H(find(~isnan(H),1,'last')); % ice thickness (m)
        Hf = interp1(x,Hf,xn,'linear','extrap'); Hf(isnan(Hf)) = Hf(find(~isnan(Hf),1,'last'));% flotation thickness (m)
        U = interp1(x,U,xn,'linear','extrap'); U(isnan(U)) = U(find(~isnan(U),1,'last')); % speed (m s^-1)
        A = interp1(x,A,x,'linear','extrap'); A(isnan(A)) = A(find(~isnan(A),1,'last'));% rate factor (Pa^-n s^-1)
        beta = interp1(x,beta,xn,'linear','extrap'); beta(isnan(beta)) = beta(find(~isnan(beta),1,'last')); % basal roughness factor

        %find the location of the grounding line and end of the ice-covered domain for the adjusted data
        gl = find(Hf-H<0,1,'last');
        ice_end = find(H<=50,1,'first');
        if isempty(ice_end) || ice_end>length(xn) || ice_end<c
            ice_end = length(xn);
            disp('ice end criteria not met.')
        end

        %rename the distance vector
        x = xn; %distance from the divide (m)

        %calculate the new surface elevation and slope
        h = hb+H; % grounded ice surface elevation (m a.s.l.)
        h(gl+1:length(x)) = (1-rho_i/rho_sw).*H(gl+1:length(x)); % floating ice surface elevation (m a.s.l.)
        dhdx = [(h(2:end)-h(1:end-1))./(x(2:end)-x(1:end-1)) 0]; % surface slope (unitless)
        if any(h<0) % surface cannot go below sea level
            c = find(h<0,1,'first');
            H(c:end) = 0;
            h(c:end)=0;
        end
        H(H>=(h-hb))=h(H>=(h-hb))-hb(H>=(h-hb)); % thickness can't go beneath bed elevation

        % calculate new strain rate
        dUdx = [(U(2:end)-U(1:end-1))./(x(2:end)-x(1:end-1)) 0]; % strain rate

        % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
        Rxx = 2*nthroot((dUdx./(E.*A(1:length(dUdx)))),n)-sigma_b; % resistive stress (Pa)
        crev = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*fwd); % crevasse penetration depth (m)
        c = find(h(1:ice_end-1)-crev(1:ice_end-1)<=0,1,'first'); % calving front located where the inland-most crevasse intersects sea level
        % if the crevasses never intersect sea level
        if isempty(c) == 1
            c = find(H<Hc,1,'first'); %set the calving front to a default minimum ice thickness value
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
            c=ci;
        end

        % calculate the misfit in surface speed and elevation change at the
        % final model time
        if t(i)==t_end
            U_misfit(f,:) = interp1(x0,interp1(x,U,x0)-U_obs(9).U',xi);
            dH_tot(f,:) = interp1(x,h,xi)-hi; 
            dH_misfit(f,:) = dH_tot(f,:)-dH_obs;
        end
    end
    
end 

% calculate mean misfit
misfit = nanmean(dH_misfit,2);

% plot results for optimal E
figure(10); clf
    IEbest = find(abs(misfit)==min(abs(misfit)),1,'first');
    Ebest = Efit(IEbest).E; 
    hold on; grid on; legend;
    set(gca,'fontsize',14,'linewidth',2); 
    xlabel('distance along centerline (km)'); ylabel('dH (m)');
    eqn = ['E = ',num2str(Efit(IEbest).fit(1)),'x + ',num2str(Efit(IEbest).fit(2))];
    plot(x./10^3,dH_tot(IEbest,:),'color',col_E(IEbest,:),'linewidth',2,...
        'displayname',eqn);
    plot(xi./10^3,dH_obs,'-k','linewidth',3,'displayname','observed dH');

% save Ebest
if save_Ebest
    cd([homepath,'inputs-outputs']);
    save('Ebest.mat','Ebest','x');
    disp('Ebest saved');
end
    
%% 5. Tune fresh water depth in crevasses fwd (via c)

close all;

save_fwdbest = 1; % = 1 to save fwdbest

% define time stepping (s)
    dt = 0.001*3.1536e7; 
    t_start = 0*3.1536e7; 
    t_end = 10*3.1536e7;    
    t = (t_start:dt:t_end);
    
% load optimal parameters and 100 yr output variables
    F0 = load('F0best.mat').F0best; % inner boundary flux (m s^-1)
    sigma_b = load('sigma_bbest.mat').sigma_bbest; % back pressure (Pa)
    E = load('Ebest.mat').Ebest; % enhancement factor
    load('Crane_flowline_100yr_output.mat');

% define fwd values to test
    fwd = 35:5:50; % fresh water depth in crevasses (m)

% pre-allocate misfit variables
    c_misfit = NaN.*zeros(length(fwd),length(2009:2018));        
        
% Loop through all fwd values
for f=1:length(fwd)
    
    % initialize variables
    x=xi; h=hi; hb=hbi; W=Wi; H=Hi; A=Ai; beta=betai; U=Ui; dUdx=dUdxi; 
    ice_end=ice_endi; c=ci;
        
    % run flowline model
    for i=1:length(t)

        % plot geometry, speed, & calving front position        
        if t(i)==t_start
            col = parula(length(t)+20); %Color scheme for plots
            figure(1); clf
            set(gcf,'Position',[0 100 1400 500]); 
            ax1 = axes('Position',[0.05 0.1 0.28 0.8]); % glacier geometry
                hold on; grid on;
                set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
                title('Glacier Geometry'); legend('Location','northeast'); 
                xlim([0 70]); ylim([min(hb)-100 max(h)+200]);
                xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)');
                % ice surface (1:c)
                plot(x(1:c)./10^3,h(1:c),'color',col(i,:),'linewidth',2,'displayname','0');
                % calving front
                plot(x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
                % ice surface (c:ice_end)
                plot(x(c+1:ice_end)./10^3,h(c+1:ice_end),'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');
                plot(x(c+1:ice_end)./10^3,h(c+1:ice_end)-H(c+1:ice_end),'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');
                % floating bed
                plot(x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off');
                % bed elevation
                plot(x/10^3,hb,'k','linewidth',2,'HandleVisibility','off');
                % mean sea level
                plot([x(1),x(end)]/10^3,[0,0],'k--','HandleVisibility','off');
            ax2 = axes('Position',[0.38 0.1 0.28 0.8]); % ice speed
                hold on; grid on;
                set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
                title('Ice Speed'); legend('Location','northeast');
                xlim([0 70]); ylim([0 4000]); 
                xlabel('Distance Along Centerline (km)'); ylabel('Speed (m yr^{-1})');
                % ice speed (1:c)
                plot(x(1:c)./10^3,U(1:c).*3.1536e7,'color',col(i,:),'linewidth',2,'displayname','0');
                % ice speed (c:ice_end)
                plot(x(c:ice_end)./10^3,U(c:ice_end).*3.1536e7,'--','color',col(i,:),'linewidth',2,'HandleVisibility','off');
            ax3 = axes('Position',[0.7 0.1 0.28 0.8]); % calving front position
                hold on; grid on;
                set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
                title('Calving Front Position'); legend('Location','best');
                xlim([35 70]); ylim([0 10]);
                xlabel('Distance Along Centerline (km)'); ylabel('Year');
                % calving front position
                plot(x(c)/10^3,t(i)./3.1536e7,'.','markersize',15,'color',col(i,:),'displayname','0');
        elseif mod(i-1,round(length(t)/50))==0 % display every length(t)/50  
            figure(1); 
            if mod(i-1,round(length(t)/10))==0 % display in legend every length(t)/10 
                % ice surface (1:c)
                plot(ax1,x(1:c)/10^3,h(1:c),'-','color',col(i,:),'linewidth',2,'displayname',num2str(t(i)./3.1536e7));
                % ice speed (1:c)
                plot(ax2,x(1:c)/10^3,U(1:c).*3.1536e7,'-','Color',col(i,:),'linewidth',2,'DisplayName',num2str(t(i)./3.1536e7));  
                % calving front position
                plot(ax3,x(c)./10^3,t(i)./3.1536e7,'.','Color',col(i,:),'markersize',15,'displayname',num2str(t(i)./3.1536e7)); hold on;                
            else
                % ice surface (1:c)
                plot(ax1,x(1:c)/10^3,h(1:c),'-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
                % ice speed (1:c)
                plot(ax2,x(1:c)/10^3,U(1:c).*3.1536e7,'-','Color',col(i,:),'linewidth',2,'HandleVisibility','off'); hold on;  
                % calving front position
                plot(ax3,x(c)./10^3,t(i)./3.1536e7,'.','Color',col(i,:),'markersize',15,'HandleVisibility','off'); hold on;                
            end
            % ice surface (c+1:ice_end)
            plot(ax1,x(c+1:ice_end)./10^3,h(c+1:ice_end),'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');            
            % calving front
            plot(ax1,x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
            % floating bed (1:c)
            plot(ax1,x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off');
            % floating bed (c+1:ice_end)
            plot(ax1,x(c+1:ice_end)./10^3,h(c+1:ice_end)-H(c+1:ice_end),'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');
            % ice speed (c:ice_end)
            plot(ax2,x(c:ice_end)./10^3,U(c:ice_end).*3.1536e7,'--','color',col(i,:),'linewidth',2,'HandleVisibility','off');
        end

        % calculate the thickness required to remain grounded at each grid cell
        Hf = -(rho_sw./rho_i).*hb; %flotation thickness (m)
        % find the location of the grounding line and use a floating
        % geometry from the grounding line to the calving front
        gl = find(Hf-H>0,1,'first')-1; %grounding line location
        ice_end = find(H<=0,1,'first'); %end of ice-covered domain
        if isempty(ice_end) || ice_end>length(x) || ice_end<c
            ice_end = length(x);
            disp('ice end criteria not met.');
        end

        %calculate the glacier's surface elevation and slope
        h = hb+H; %h = surface elevation (m a.s.l.)
        h(gl+1:length(x)) = (1-rho_i/rho_sw).*H(gl+1:length(x)); %adjust the surface elevation of ungrounded ice to account for buoyancy
        dhdx = [(h(2:end)-h(1:end-1))./(x(2:end)-x(1:end-1)) 0]; % surface slope (unitless)
        if any(h<0) % surface cannot go below sea level
            H(c:end) = 0;
            h(c:end)=0;
        end

        % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
        Rxx = 2*nthroot((dUdx./(E.*A(1:length(dUdx)))),n)-sigma_b; % resistive stress (Pa)
        crev = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*fwd(f)); % crevasse penetration depth (m)
        c = find(h(1:ice_end-1)-crev(1:ice_end-1)<=0,1,'first'); % calving front located where the inland-most crevasse intersects sea level
        %if the crevasses never intersect sea level
        if isempty(c) == 1
            c = find(H<Hc,1,'first'); %set the calving front to a default minimum ice thickness value
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
            c=ci;
        end

        % calculate the effective pressure (ice overburden pressure minus water
        % pressure) assuming an easy & open connection between the ocean and
        % ice-bed interface
        sl = find(hb<=0,1,'first'); % find where the glacier base first drops below sea level
        N_ground = rho_i*g*H(1:sl); % effective pressure where the bed is above sea level (Pa)
        N_marine = rho_i*g*H(sl+1:length(x))+(rho_sw*g*hb(sl+1:length(x))); % effective pressure where the bed is below sea level (Pa)
        N = [N_ground N_marine];
        N(N<0)=1; % cannot have negative values

        % Solve for new velocity
        [U,~,vm,T] = U_convergence(x,U,U0,dUdx,dhdx,H,A,E,N,W,dx,c,ice_end,n,m,beta,rho_i,rho_sw,g);

        % calculate ice flux
        F = U.*H.*W; % ice flux (m^3 s^-1)
        F(isnan(F))=0;
        F(1) = F0; 

        % calculate the  change in ice thickness from continuity
        dHdt = -(1./W).*gradient(F,x);
        dH = dHdt.*dt;

        % implement SMB & SMR
        % calculate SMR where ice is ungrounded and fit a smoothing spline 
        % past the grounding line to avoid large spatial gradients
        smr = zeros(1,length(x));
            smr(gl+1:ice_end) = feval(fit([x(gl); x(c); x(ice_end)],...
                [0;smr0.*[1;1]],'poly1'),x(gl+1:ice_end));
            smr(ice_end+1:length(x)) = smr(ice_end); 
            smr = movmean(smr,10);
        smb = interp1(x0,smb0+Q0,x); 
            smb = movmean(smb,10);

        % new thickness (change from dynamics, SMB, & SMR)
        Hn = H+dH+(smb.*dt)+(smr.*dt);
        Hn(Hn < 0) = 0; % remove negative values
        H = Hn; % set as the new thickness value
        % smooth out the points near the ice divide
        H(1:10) = ones(1,length(H(1:10))).*nanmean(H(1:10),'all');

        % thickness & surface past calving front
        for j=c+1:length(xi)
            h(j) = h(j-1)-5; % decrease until reaching 0m
            H(j) = H(j-1)-50; % decrease until reaching 0m
            if H(j)>=h(j)-hb(j)
                H(j)=h(j)-hb(j); % can't go beneath bed elevation
            end
        end
        h(h<0)=0; % surface can't go below sea level
        H(H<0)=0; % no negative thicknesses

        % stop the model if it behaves unstably (monitored by ice thickness and speed)
        if max(H) > H_max
            disp(['Adjust dt']);
            break;
        end
        if mean(U) < 200/3.1536e7
            disp('Too slow!');
            break;
        end

        % find the precise location of the grounding line (where H=Hf)
        xf = x(find(Hf-H>0,1,'first')-1);

        %adjust the grid spacing so the grounding line is continuously tracked
        xl = round(xf/dx0); % number of ideal grid spaces needed to reach the grounding line
        dx = xf/xl; % new grid spacing (should be ~dx0)
        xn = 0:dx:L; % new distance vector

        %adjust the space-dependent variables to the new distance vector
        hb = interp1(x,hb,xn); hb(isnan(hb)) = hb(find(~isnan(hb),1,'last')); % glacier bed elevation (m)
        W = interp1(x,W,xn); W(isnan(W)) = W(find(~isnan(W),1,'last')); % glacier width (m)
        H = interp1(x,H,xn,'linear','extrap'); H(isnan(H)) = H(find(~isnan(H),1,'last')); % ice thickness (m)
        Hf = interp1(x,Hf,xn,'linear','extrap'); Hf(isnan(Hf)) = Hf(find(~isnan(Hf),1,'last'));% flotation thickness (m)
        U = interp1(x,U,xn,'linear','extrap'); U(isnan(U)) = U(find(~isnan(U),1,'last')); % speed (m s^-1)
        A = interp1(x,A,x,'linear','extrap'); A(isnan(A)) = A(find(~isnan(A),1,'last'));% rate factor (Pa^-n s^-1)
        beta = interp1(x,beta,xn,'linear','extrap'); beta(isnan(beta)) = beta(find(~isnan(beta),1,'last')); % basal roughness factor

        %find the location of the grounding line and end of the ice-covered domain for the adjusted data
        gl = find(Hf-H<0,1,'last');
        ice_end = find(H<=0,1,'first');
        if isempty(ice_end) || ice_end>length(xn) || ice_end<c
            ice_end = length(xn);
            disp('ice end criteria not met.')
        end

        %rename the distance vector
        x = xn; %distance from the divide (m)

        %calculate the new surface elevation and slope
        h = hb+H; % grounded ice surface elevation (m a.s.l.)
        h(gl+1:length(x)) = (1-rho_i/rho_sw).*H(gl+1:length(x)); % floating ice surface elevation (m a.s.l.)
        dhdx = [(h(2:end)-h(1:end-1))./(x(2:end)-x(1:end-1)) 0]; % surface slope (unitless)
        if any(h<0) % surface cannot go below sea level
            c = find(h<0,1,'first');
            H(c:end) = 0;
            h(c:end)=0;
        end
        H(H>=(h-hb))=h(H>=(h-hb))-hb(H>=(h-hb)); % thickness can't go beneath bed elevation

        % calculate new strain rate
        dUdx = [(U(2:end)-U(1:end-1))./(x(2:end)-x(1:end-1)) 0]; % strain rate

        % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
        Rxx = 2*nthroot((dUdx./(E.*A(1:length(dUdx)))),n)-sigma_b; % resistive stress (Pa)
        crev = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*fwd(f)); % crevasse penetration depth (m)
        c = find(h(1:ice_end-1)-crev(1:ice_end-1)<=0,1,'first'); % calving front located where the inland-most crevasse intersects sea level
        % if the crevasses never intersect sea level
        if isempty(c) == 1
            c = find(H<Hc,1,'first'); %set the calving front to a default minimum ice thickness value
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
            c=ci;
        end

        % calculate misfit in calving front position at each full model year
        if mod(t(i),3.1536e7)==0 && t(i)~=t_start
            c_misfit(f,t(i)./3.1536e7) = x(c)-interp1(termDate_obs,termx_obs,t(i)./3.1536e7+2009);
        end
        
    end 

end

% Calculate total misfit
misfit = mean(c_misfit,2);

% Display results
figure(11); clf
    set(gca,'fontsize',14,'fontname','Arial','linewidth',2);
    grid on; xlabel('fwd (m)'); ylabel('misfit'); hold on;
    for i=1:length(misfit)
        if any(isnan(c_misfit(i,:)))
            plot(fwd(i),misfit(i),'*r','markersize',15,'linewidth',2);
        else
            plot(fwd(i),misfit(i),'*b','markersize',15,'linewidth',2);            
        end
    end

% save results
if save_fwdbest
    cd([homepath,'inputs-outputs']);
    fwdbest = fwd(abs(misfit)==min(abs(misfit))); % m.w.e.
    save('fwdbest.mat','fwdbest');
    disp(['Best fwd = ',num2str(fwdbest),' m saved.']);
end
    
%% 6. Run sensitivity tests for SMB & SMR
% (Using 100 yr output as initialization)

% Run through simulated 2009-2019 with best solutions for parameters, 
%   then continue evolving the model for 10 years under new scenarios:
%   \Delta(SMB) &/or \Delta(SMR)

close all;

save_figure = 0;    % = 1 to save resulting figure
save_final = 0;     % = 1 to save final geometry and speed 

% time stepping (s)
    dt = 0.001*3.1536e7; 
    t_start = 0*3.1536e7; 
    t_end = 10*3.1536e7;    
    t = (t_start:dt:t_end);

% load optimal parameters and 100 yr output variables
    F0 = load('F0best.mat').F0best; % optimal inner boundary flux from step #2
    sigma_b = load('sigma_bbest.mat').sigma_bbest; % optimal back stress from step #3
    E = load('Ebest.mat').Ebest; % optimal enhancement factor from step #4
    fwd = load('fwdbest.mat').fwdbest; % optimal fresh water depth from step #5
    load('Crane_sensitivityTests_noChange.mat'); % load no change variables
    load('Crane_flowline_100yr_output.mat');

% set up changes in SMB & SMR
% note: decrease SMR in increments of 0.5 m a-1 (1.585e-8)
%   until reaching SMR found on at other Antarctic ice shelves: ~6 m a^-1 
%   (Adusumilli et al., 2020)
    delta_smb = 0/3.1536e7; % m/s change in SMB
    delta_smr = -6/3.1536e7; % m/s change in SMR
    
% initialize variables using final conditions from 100 yr scenario
    x=xi; h=hi; hb=hbi; W=Wi; H=Hi; A=Ai; beta=betai; U=Ui; dUdx=dUdxi; 
    ice_end=ice_endi; c=ci;

% run flowline model
for i=1:length(t)

    % plot geometry, speed, & calving front position        
    if t(i)==t_start
        col = parula(length(t)+20); %Color scheme for plots
        figure(1); clf
        set(gcf,'Position',[0 100 1300 400]); 
        ax1 = axes('Position',[0.05 0.1 0.28 0.8]); % glacier geometry
            hold on; grid on;
            set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
            title('Glacier Geometry'); legend('Location','northeast'); 
            xlim([0 70]); ylim([min(hb)-100 max(h)+200]);
            xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)');
            % ice surface (1:c)
            plot(x(1:c)./10^3,h(1:c),'color',col(i,:),'linewidth',2,'displayname','0');
            % calving front
            plot(x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
            % ice surface (c:ice_end)
            plot(x(c+1:ice_end)./10^3,h(c+1:ice_end),'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');
            plot(x(c+1:ice_end)./10^3,h(c+1:ice_end)-H(c+1:ice_end),'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');
            % floating bed
            plot(x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off');
            % bed elevation
            plot(x/10^3,hb,'k','linewidth',2,'HandleVisibility','off');
            % mean sea level
            plot([x(1),x(end)]/10^3,[0,0],'k--','HandleVisibility','off');
        ax2 = axes('Position',[0.38 0.1 0.28 0.8]); % ice speed
            hold on; grid on;
            set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
            title('Ice Speed'); legend('Location','northeast');
            xlim([0 70]); ylim([0 4000]); 
            xlabel('Distance Along Centerline (km)'); ylabel('Speed (m yr^{-1})');
            % ice speed (1:c)
            plot(x(1:c)./10^3,U(1:c).*3.1536e7,'color',col(i,:),'linewidth',2,'displayname','0');
            % ice speed (c:ice_end)
            plot(x(c:ice_end)./10^3,U(c:ice_end).*3.1536e7,'--','color',col(i,:),'linewidth',2,'HandleVisibility','off');
        ax3 = axes('Position',[0.7 0.1 0.28 0.8]); % calving front position
            hold on; grid on;
            set(gca,'FontSize',12,'linewidth',2,'fontweight','bold');
            title('Calving Front Position'); legend('Location','best');
            xlim([35 70]); ylim([0 20]);
            xlabel('Distance Along Centerline (km)'); ylabel('Year');
            % calving front position
            plot(x(c)/10^3,t(i)./3.1536e7,'.','markersize',15,'color',col(i,:),'displayname','0');
        figure(4); clf; 
        subplot(1,2,1); hold on;
            set(gca,'fontsize',12,'fontweight','bold','linewidth',2); grid on;
            xlabel('Distance Along Centerline (km)'); ylabel('SMB (m a^{-1})');
        subplot(1,2,2); hold on;
            set(gca,'fontsize',12,'fontweight','bold','linewidth',2); grid on;
            xlabel('Distance Along Centerline (km)'); ylabel('SMR (m a^{-1})');        
    elseif mod(i-1,round(length(t)/50))==0 % display every length(t)/50  
        figure(1); 
        if mod(i-1,round(length(t)/10))==0 % display in legend every length(t)/10 
            % ice surface (1:c)
            plot(ax1,x(1:c)/10^3,h(1:c),'-','color',col(i,:),'linewidth',2,'displayname',num2str(t(i)./3.1536e7));
            % ice speed (1:c)
            plot(ax2,x(1:c)/10^3,U(1:c).*3.1536e7,'-','Color',col(i,:),'linewidth',2,'DisplayName',num2str(t(i)./3.1536e7));  
            % calving front position
            plot(ax3,x(c)./10^3,t(i)./3.1536e7,'.','Color',col(i,:),'markersize',15,'displayname',num2str(t(i)./3.1536e7)); hold on;                
            figure(4); 
            subplot(1,2,1);
                plot(x./10^3,smb.*3.1536e7,'linewidth',2,'color',col(i,:));
            subplot(1,2,2);
                plot(x./10^3,smr.*3.1536e7,'linewidth',2,'color',col(i,:));
        else
            % ice surface (1:c)
            plot(ax1,x(1:c)/10^3,h(1:c),'-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
            % ice speed (1:c)
            plot(ax2,x(1:c)/10^3,U(1:c).*3.1536e7,'-','Color',col(i,:),'linewidth',2,'HandleVisibility','off'); hold on;  
            % calving front position
            plot(ax3,x(c)./10^3,t(i)./3.1536e7,'.','Color',col(i,:),'markersize',15,'HandleVisibility','off'); hold on;                
        end
        % ice surface (c+1:ice_end)
        plot(ax1,x(c+1:ice_end)./10^3,h(c+1:ice_end),'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');            
        % calving front
        plot(ax1,x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
        % floating bed (1:c)
        plot(ax1,x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off');
        % floating bed (c+1:ice_end)
        plot(ax1,x(c+1:ice_end)./10^3,h(c+1:ice_end)-H(c+1:ice_end),'--','color',col(i,:),'linewidth',1.5,'HandleVisibility','off');
        % ice speed (c:ice_end)
        plot(ax2,x(c:ice_end)./10^3,U(c:ice_end).*3.1536e7,'--','color',col(i,:),'linewidth',2,'HandleVisibility','off');
    end
    
    % calculate the thickness required to remain grounded at each grid cell
    Hf = -(rho_sw./rho_i).*hb; % flotation thickness (m)
    % find the location of the grounding line and use a floating
    % geometry from the grounding line to the calving front
    gl = find(Hf-H>0,1,'first')-1; % grounding line location
    ice_end = find(H<=0,1,'first'); % end of ice-covered domain
    if isempty(ice_end) || ice_end>length(x) || ice_end<c
        ice_end = length(x);
        disp('ice end criteria not met.');
    end

    %calculate the glacier's surface elevation and slope
    h = hb+H; %h = surface elevation (m a.s.l.)
    h(gl+1:length(x)) = (1-rho_i/rho_sw).*H(gl+1:length(x)); %adjust the surface elevation of ungrounded ice to account for buoyancy
    dhdx = [(h(2:end)-h(1:end-1))./(x(2:end)-x(1:end-1)) 0]; % surface slope (unitless)
    if any(h<0) % surface cannot go below sea level
        H(c:end) = 0;
        h(c:end)=0;
    end

    % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
    Rxx = 2*nthroot((dUdx./(E.*A(1:length(dUdx)))),n)-sigma_b; % resistive stress (Pa)
    crev = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*fwd); % crevasse penetration depth (m)
    c = find(h(1:ice_end-1)-crev(1:ice_end-1)<=0,1,'first'); % calving front located where the inland-most crevasse intersects sea level
    %if the crevasses never intersect sea level
    if isempty(c) == 1
        c = find(H<Hc,1,'first'); %set the calving front to a default minimum ice thickness value
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
        c=ci;
    end

    % calculate the effective pressure (ice overburden pressure minus water
    % pressure) assuming an easy & open connection between the ocean and
    % ice-bed interface
    sl = find(hb<=0,1,'first'); % find where the glacier base first drops below sea level
    N_ground = rho_i*g*H(1:sl); % effective pressure where the bed is above sea level (Pa)
    N_marine = rho_i*g*H(sl+1:length(x))+(rho_sw*g*hb(sl+1:length(x))); % effective pressure where the bed is below sea level (Pa)
    N = [N_ground N_marine];
    N(N<0)=1; % cannot have negative values

    % Solve for new velocity
    [U,~,vm,T] = U_convergence(x,U,U0,dUdx,dhdx,H,A,E,N,W,dx,c,ice_end,n,m,beta,rho_i,rho_sw,g);

    % calculate ice flux
    F = U.*H.*W; % ice flux (m^3 s^-1)
    F(isnan(F))=0;
    F(1) = F0; 

    % calculate the  change in ice thickness from continuity
    dHdt = -(1./W).*gradient(F,x);
    dH = dHdt.*dt;

    % implement SMB & SMR
    % calculate SMR where ice is ungrounded and fit a smoothing spline 
    % past the grounding line to avoid large spatial gradients
    if t(i)/3.1536e7>=10
        smr = zeros(1,length(x));
            smr(gl+1:ice_end) = feval(fit([x(gl); x(c); x(ice_end)],...
                [0;smr0+delta_smr.*[1;1]],'poly1'),x(gl+1:ice_end));
            smr(ice_end+1:length(x)) = smr(ice_end); 
        smb = interp1(x0,smb0+Q0+delta_smb,x); 
    else
        smr = zeros(1,length(x));
            smr(gl+1:ice_end) = feval(fit([x(gl); x(c); x(ice_end)],...
                [0;smr0.*[1;1]],'poly1'),x(gl+1:ice_end));
            smr(ice_end+1:length(x)) = smr(ice_end); 
            smr = movmean(smr,10);
        smb = interp1(x0,smb0+Q0,x); 
    end

    % new thickness (change from dynamics, SMB, & SMR)
    Hn = H+dH+(smb.*dt)+(smr.*dt);
    Hn(Hn < 0) = 0; % remove negative values
    H = Hn; % set as the new thickness value
    % smooth out the points near the ice divide
    H(1:10) = ones(1,length(H(1:10))).*nanmean(H(1:10),'all');

    % thickness & surface past calving front
%     for j=c+1:length(xi)
%         h(j) = h(j-1)-5; % decrease until reaching 0m
%         H(j) = H(j-1)-50; % decrease until reaching 0m
%         if H(j)>=h(j)-hb(j)
%             H(j)=h(j)-hb(j); % can't go beneath bed elevation
%         end
%     end
    h(h<0)=0; % surface can't go below sea level
    H(H<0)=0; % no negative thicknesses

    % stop the model if it behaves unstably (monitored by ice thickness and speed)
    if max(H) > H_max
        disp(['Adjust dt']);
        break;
    end
    if mean(U) < 200/3.1536e7
        disp('Too slow!');
        break;
    end

    % find the precise location of the grounding line (where H=Hf)
    xf = x(find(Hf-H>0,1,'first')-1);

    %adjust the grid spacing so the grounding line is continuously tracked
    xl = round(xf/dx0); % number of ideal grid spaces needed to reach the grounding line
    dx = xf/xl; % new grid spacing (should be ~dx0)
    xn = 0:dx:L; % new distance vector

    %adjust the space-dependent variables to the new distance vector
    hb = interp1(x,hb,xn); % glacier bed elevation (m)
    W = interp1(x,W,xn); W(isnan(W)) = W(find(~isnan(W),1,'last')); % glacier width (m)
    H = interp1(x,H,xn,'linear','extrap'); H(isnan(H)) = H(find(~isnan(H),1,'last')); % ice thickness (m)
    Hf = interp1(x,Hf,xn,'linear','extrap'); Hf(isnan(Hf)) = Hf(find(~isnan(Hf),1,'last'));% flotation thickness (m)
    U = interp1(x,U,xn,'linear','extrap'); U(isnan(U)) = U(find(~isnan(U),1,'last')); % speed (m s^-1)
    A = interp1(x,A,x,'linear','extrap'); A(isnan(A)) = A(find(~isnan(A),1,'last'));% rate factor (Pa^-n s^-1)
    beta = interp1(x,beta,xn,'linear','extrap'); beta(isnan(beta)) = beta(find(~isnan(beta),1,'last')); % basal roughness factor

    %find the location of the grounding line and end of the ice-covered domain for the adjusted data
    gl = find(Hf-H<0,1,'last');
    ice_end = find(H<=0,1,'first');
    if isempty(ice_end) || ice_end>length(xn) || ice_end<c
        ice_end = length(xn);
        disp('ice end criteria not met.')
    end

    %rename the distance vector
    x = xn; %distance from the divide (m)

    %calculate the new surface elevation and slope
    h = hb+H; % grounded ice surface elevation (m a.s.l.)
    h(gl+1:length(x)) = (1-rho_i/rho_sw).*H(gl+1:length(x)); % floating ice surface elevation (m a.s.l.)
    dhdx = [(h(2:end)-h(1:end-1))./(x(2:end)-x(1:end-1)) 0]; % surface slope (unitless)
    if any(h<0) % surface cannot go below sea level
        c = find(h<0,1,'first');
        H(c:end) = 0;
        h(c:end)=0;
    end
    H(H>=(h-hb))=h(H>=(h-hb))-hb(H>=(h-hb)); % thickness can't go beneath bed elevation

    % calculate new strain rate
    dUdx = [(U(2:end)-U(1:end-1))./(x(2:end)-x(1:end-1)) 0]; % strain rate

    % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
    Rxx = 2*nthroot((dUdx./(E.*A(1:length(dUdx)))),n)-sigma_b; % resistive stress (Pa)
    crev = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*fwd); % crevasse penetration depth (m)
    c = find(h(1:ice_end-1)-crev(1:ice_end-1)<=0,1,'first'); % calving front located where the inland-most crevasse intersects sea level
    % if the crevasses never intersect sea level
    if isempty(c) == 1
        c = find(H<Hc,1,'first'); %set the calving front to a default minimum ice thickness value
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
        c=ci;
    end

end 

% plot results on final time step
if t(i)==t_end
    h2=h; H2=H; x2=x; c2=c; gl2=gl; ice_end2=ice_end; % store final geometry
    figure(10); clf % sensitivity test changes
    hold on; grid on;
    set(gcf,'Position',[350 300 700 600]);
    set(gca,'FontSize',14,'linewidth',2,'fontweight','bold');
    xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)');
    legend('Location','east'); xlim([0 70]); ylim([-1200 1200]);
    title(['SMR = + ',num2str(round(delta_smr.*3.1536e7,1)),'m/a, SMB = + ',...
        num2str(round(delta_smb*3.1536e7,1)),'m/a']);
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
        plot(x1/10^3,hb,'k','linewidth',2,'HandleVisibility','off');
    % inset plot of terminus
    ax2 = axes('Position',[0.62 0.62 0.28 0.28]);
        hold on; grid on; set(gca,'linewidth',2,'fontweight','bold');
        title(['\Delta L = ',num2str(x2(c2)-x1(c1)),'m, ','\Delta H_{mean} = ',...
            num2str(round(mean(H2(1:c2))-mean(H1(1:c1)),1)),' m']);
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
    cd([homepath,'scripts/modelingWorkflow/results']);
    % Save figure for test
    figName = ['SMR',num2str(delta_smr),'_SMB',num2str(delta_smb)];
    saveas(gcf,[figName,'.png'],'png');
    saveas(gcf,[figName,'.fig'],'fig');
    disp(['Fig. 10 saved (.png & .fig) in: ',pwd]);
elseif ~ishandle(4)
    disp('figure does not exist.');
else
    disp('no figure saved.');
end

% save geometry
if save_final
    cd([homepath,'scripts/modelingWorkflow/results']);
    fileName=['SMR',num2str(delta_smr),'_SMB',num2str(delta_smb),'_geom.mat'];
    h2=h; H2=H; gl2=gl; c2=c; ice_end2=ice_end; U2=U; x2=x; 
    save(fileName,'h2','H2','c2','c2','ice_end2','U2','gl2','x2');
    disp('geometry saved.');
else
    disp('geometry not saved.');
end
