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
    addpath([homepath,'scripts/tuningSMB']); % add path to U_convergence

% Load Crane Glacier initialization variables
    load('Crane_flowline_initialization.mat');
    A0(end+1:length(x0)) = A0(end).*ones(1,length(x0)-length(A0)); 
    
% Load observations of dH to help tune SMB
    dH_obs = load('dHdt.mat').dHdt.dH_total; % (m) total change in thickness 2009-2018
    
% densities and g
    rho_i = 917; % ice density (kg m^-3)
    rho_sw = 1028; % ocean water density (kg m^-3)
    rho_fw = 1000; % fresh water density (kg m^-3)
    g = 9.81; % acceleration (m s^-2)

% time stepping (s)
    dt = 0.05*3.1536e7; 
    t_start = 0*3.1536e7; 
    t_end = 10*3.1536e7;    
    t = (t_start:dt:t_end);
    years = 2009:(2009+t_end/3.31536e7);

% stress parameters (unitless)
    m = 1; % basal sliding exponent
    n = 3; % flow law exponent
    E = 1; % enhancement factor
    fwd = 30; % fresh water depth in crevasses (m)
    
% maximum thickness cut-off to check for instability
    H_max = 2000; %maximum thickness (m)

% regrid the initialization data to match the desired grid spacing, 
    xi = x0; % rename initialization distance vector
    L = 57e3; % length of glacier (m)
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
    for i=1:c % thickness can't go beneath bed elevation
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

for i=1:length(t)
    
    % set up figures, plot geometries at t==0, then every t/10 iterations
    if t(i)==t_start
        col = parula(length(t)+10); %Color scheme for plots
        figure; % glacier geometry
            hold on; grid on;
            set(gcf,'Position',[0 50 500 400]);
            set(gca,'FontSize',14,'linewidth',2,'fontweight','bold'); 
            legend('Location','best'); xlim([0 65]); ylim([min(hb)-100 max(h)+100]);
            title('a) Glacier Geometry'); 
            xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)'); 
            % ice surface
            plot(x(1:c)./10^3,h(1:c),'color',col(i,:),'linewidth',2,'displayname','2009');
            % calving front
            plot(x(c)*[1,1]/10^3,[h(c)-H(c),h(c)],'.-','color',col(i,:),'linewidth',2,'HandleVisibility','off');
            % floating bed
            plot(x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off');
            % bed elevation
            plot(x/10^3,hb,'k','linewidth',2,'HandleVisibility','off'); 
            % mean sea level
            plot([x(1),x(end)]/10^3,[0,0],'k--','HandleVisibility','off'); 
        figure; % ice speed
            hold on; grid on; 
            set(gcf,'Position',[500 50 500 400]);
            set(gca,'FontSize',14,'linewidth',2,'fontweight','bold'); 
            title('b) Ice Speed Profile'); 
            xlabel('Distance Along Centerline (km)'); ylabel('Speed (m yr^{-1})'); 
            %legend('Location','northwest'); 
            plot(x(1:c)./10^3,U(1:c).*3.1536e7,'color',col(i,:),'linewidth',2,'displayname','2009');
        figure; % terminus position     
            hold on; grid on; 
            set(gcf,'Position',[1000 50 500 400]);
            set(gca,'FontSize',14,'linewidth',2,'fontweight','bold');
            title('c) Terminus Position'); ylabel('Year');
            xlabel('Distance Along Centerline (km)'); 
            %legend('Location','best'); 
            % 2009 terminus position
            plot(x(c)./10^3,2009,'*','markersize',10,'color',col(i,:),...
                'linewidth',1.5,'displayname','2009'); 
        figure; % dH
            set(gcf,'Position',[450 350 550 450]);
            hold on; grid on; legend;
            set(gca,'FontSize',14,'linewidth',2); grid on;
            xlabel('Distance Along Centerline (km)'); ylabel('dH (m)'); 
            title('dH (2009-2018)');
            plot(x0./10^3,dH_obs,'-m','linewidth',2,'displayname','observed');
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
            % floating bed
            plot(x(gl:c)/10^3,h(gl:c)-H(gl:c),'color',col(i,:),'linewidth',2,'HandleVisibility','off');
        end 
        
        figure(2); hold on; % Plot velocity every 10 time iterations
        if mod(i-1,round(length(t)/10))==0
            plot(x(1:c)/10^3,U(1:c).*3.1536e7,'-','Color',col(i,:),'linewidth',2,'DisplayName',num2str(t(i)./3.1536e7+2009)); hold on;
        else 
            plot(x(1:c)/10^3,U(1:c).*3.1536e7,'-','Color',col(i,:),'linewidth',2,'HandleVisibility','off'); hold on;           
        end
        
        figure(3); hold on; % Plot terminus position every 10 time iterations
        if mod(i-1,round(length(t)/10))==0
            plot(x(c)./10^3,t(i)./3.1536e7+2009,'*','Color',col(i,:),'markersize',15,'linewidth',1.5,'DisplayName',num2str(t(i)./3.1536e7+2009)); hold on;
        else 
            plot(x(c)./10^3,t(i)./3.1536e7+2009,'.','Color',col(i,:),'markersize',10,'linewidth',1.5,'HandleVisibility','off'); hold on;           
        end          
    end 

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
        crev = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*fwd); %crevasse penetration depth (m)
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
        crev = (Rxx./(rho_i.*g))+((rho_fw./rho_i).*fwd); %crevasse penetration depth (m)
        c = find(h(1:ice_end-1)-crev(1:ice_end-1)<=0,1,'first'); %calving front located where the inland-most crevasse intersects sea level
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
            c=dsearchn(transpose(x),x0(term_2009.x));            
        end 
      
end 