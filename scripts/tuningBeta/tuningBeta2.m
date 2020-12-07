%% flowlinemodel_Crane_tuningBeta

clear all; close all;
warning off; % turn off warnings (velocity coefficient matrix is close to singular)

%% define time and space independent variables
    
dx0 = 500:500:7e3; % desired grid spacings
    % Stress-coupling Length (SCL) =~ 5km

mu_RMSE = zeros(1,length(dx0));

figure_save = 1; % =1 to save figures for each spatial resolution

for z=1:length(dx0)
    
    use_binavg = 1;     % = 1 to use average within bins proportional to dx0
    poly_solve = 1;     % = 1 to solve for the best polynomial degree fit for beta

    % define home path in directory
    homepath = '/Users/raineyaberle/Desktop/Research/CraneGlacier_modeling/';
    cd([homepath,'inputs-outputs']);

    % add paths to required matlab functions
    addpath('/Users/raineyaberle/Desktop/Research/matlabFunctions');
    addpath([homepath,'scripts/tuningbeta']);

    % Load Crane Glacier initialization variables
    load('Crane_flowline_initialization.mat');
        A0 = feval(fit(x0',A0','poly1'),x0)';
        W0(isnan(W0)) = W0(find(~isnan(W0),1,'last'));
        smr = 4.75e-8;
        fwd = 30; % fresh water depth in crevasses (m)
        Hc = 400; % m
        H_max=2000; % m

        % Load 2009 terminus position
        x_cl = load('Crane_centerline.mat').x; 
        y_cl = load('Crane_centerline.mat').y;     
        termx = load('LarsenB_centerline.mat').centerline.termx;
        termy = load('LarsenB_centerline.mat').centerline.termy;    
        term = dsearchn([x_cl y_cl],[termx(5) termy(5)]);
            clear x_cl y_cl termx termy

        % densities and g
        rho_i = 917; % ice density (kg m^-3)
        rho_sw = 1028; % ocean water density (kg m^-3)
        g = 9.81; % acceleration (m s^-2)

        % time stepping (s)
        dt = 0.5*3.1536e7; 
        t_start = 0*3.1536e7; 
        t_end = 0.5*3.1536e7;    
        t = (t_start:dt:t_end);

        % stress parameters (unitless)
        m = 1; % basal sliding exponent
        n = 3; % flow law exponent
        E = 1; % enhancement factor

    % regrid the initialization data using the interp1 function to match the desired grid spacing
        x_i = x0; % rename initialization distance vector
        L = x0(term); % length of glacier (m)
        x = 0:dx0(z):L; % desired distance vector (m from ice divide)
        h = interp1(x_i,h0,x,'linear','extrap');
        hb = interp1(x_i,hb0,x,'linear','extrap');
        W = interp1(x_i,W0,x,'linear','extrap');
        U = interp1(x_i,U0,x,'linear','extrap');
        A = interp1(x_i,A0,x,'linear','extrap');
        SMB0 = interp1(x_i,smb0,x,'linear','extrap');
        dUdx = (U(2:end)-U(1:end-1))./(x(2:end)-x(1:end-1));
            dUdx = [0,dUdx];

        beta = zeros(length(t),length(x)); % (s^{1/m} m^{-1/m})

    % Use a staggered grid to calculate averages within bins
        if use_binavg
            xm = (x(1:end-1)+x(2:end))./2;
            for k=1:length(x)
                if k==1
                  h(k) = mean(h0(1:dsearchn(x0',xm(k))));
                  hb(k) = mean(hb0(1:dsearchn(x0',xm(k))));
                  W(k) = mean(W0(1:dsearchn(x0',xm(k))));
                  U(k) = mean(U0(1:dsearchn(x0',xm(k))));
                  A(k) = mean(A0(1:dsearchn(x0',xm(k))));
                  SMB0(k) = mean(smb0(1:dsearchn(x0',xm(k)))); 
                elseif k==length(x)
                  h(k) = mean(h0(dsearchn(x0',xm(k-1)):term));
                  hb(k) = mean(hb0(dsearchn(x0',xm(k-1)):term));
                  W(k) = mean(W0(dsearchn(x0',xm(k-1)):term));
                  U(k) = mean(U0(dsearchn(x0',xm(k-1)):term));
                  A(k) = mean(A0(dsearchn(x0',xm(k-1)):term));
                  SMB0(k) = mean(smb0(dsearchn(x0',xm(k-1)):term));
                else
                  h(k) = mean(h0(dsearchn(x0',xm(k-1)):dsearchn(x0',xm(k))));
                  hb(k) = mean(hb0(dsearchn(x0',xm(k-1)):dsearchn(x0',xm(k))));
                  W(k) = mean(W0(dsearchn(x0',xm(k-1)):dsearchn(x0',xm(k))));
                  U(k) = mean(U0(dsearchn(x0',xm(k-1)):dsearchn(x0',xm(k))));
                  A(k) = mean(A0(dsearchn(x0',xm(k-1)):dsearchn(x0',xm(k))));
                  SMB0(k) = mean(smb0(dsearchn(x0',xm(k-1)):dsearchn(x0',xm(k))));
                end
            end 
            H = h-hb;
        end 

        SMB=SMB0; dx=dx0(z); ice_end = length(x);

    % find the grounding line
    Hf = -(rho_sw./rho_i).*hb; % flotation thickness (m)
    gl = find(Hf-H>0,1,'first')-1; %grounding line location

    % add floating thickness past grounding line
    H(gl:end) = h(gl:end)*rho_sw/(rho_sw-rho_i); % buoyant thickness using surface
    H(H>h-hb) = h(H>h-hb)-hb(H>h-hb); % thickness cannot go beneath bed

    % set up figures
    col = parula(length(dx0)+5); %Color scheme for plots
    figure(1); % fig1 = glacier geometry
        set(gcf,'Position',[0 100 500 400]);
        set(gca,'FontSize',14,'FontName','Arial'); grid on;
        legend; xlim([0 60]); ylim([min(hb)-100 max(h)+100]);
        xlabel('Distance Along Centerline (km)'); ylabel('Elevation (m)'); 
        hold on; grid on;
        plot(x(end)*[1,1]/10^3,[hb(end),h(end)],'.-','color',col(z,:),'linewidth',2,'HandleVisibility','off');
        if z==1
            plot(x/10^3,hb,'k','linewidth',2,'HandleVisibility','off'); 
        end
        plot([x(end),x0(end)]/10^3,[0,0],'k--','displayname','Mean Sea Level','HandleVisibility','off'); 
    figure(2); % fig2 = ice speed
        hold on; grid on; set(gcf,'Position',[500 100 500 400]);
        set(gca,'FontSize',14,'FontName','Arial','linewidth',2); grid on;
        xlabel('Distance Along Centerline (km)'); ylabel('Speed (m yr^-^1)'); 
        legend('Location','northwest'); title('Ice Speed Profile');
    figure(3); % fig3 = beta results
    set(gcf,'Position',[0 300 1200 500]);
    subplot(1,3,1); % beta along x
        hold on; grid on;
        set(gca,'FontSize',14,'FontName','Arial','linewidth',2); grid on;
        xlabel('Distance Along Centerline (km)'); ylabel('\beta (s^{1/m} m^{-1/m})'); 
        title('Basal roughness factor');
        xlim([0 50]); ylim([-1e4 1e4]);
    subplot(1,3,2); % best poly fit
        hold on; grid on;
        set(gca,'FontSize',14,'FontName','Arial','linewidth',2); grid on; 
        xlabel('Distance Along Centerline (km)'); ylabel('\beta (s^{1/m} m^{-1/m})');
        title('\beta Polynomial Fits');
        hold on; grid on;
    subplot(1,3,3); % RMSE
        hold on; grid on; 
        set(gca,'FontSize',14,'FontName','Arial','linewidth',2); grid on;
        xlabel('dx (km)'); ylabel('RMSE (s^{1/m} m^{-1/m})');
        title('RMSE of poly fits');
        xlim([0 8]); ylim([0 6e9]);
        hold on; grid on; legend;       
        
    for i=1 % only run through one time increment
    
        % calculate the thickness required to remain grounded at each grid cell
        Hf = -(rho_sw./rho_i).*hb; %flotation thickness (m)
        % find the location of the grounding line and use a floating
        % geometry from the grounding line to the calving front
        gl = find(Hf-H>0,1,'first')-1; %grounding line location 

        %calculate the glacier's surface elevation and slope
        h = hb+H; %h = surface elevation (m a.s.l.)
        h(gl:length(x)) = (1-rho_i/rho_sw).*H(gl:length(x)); %adjust the surface elevation of ungrounded ice to account for buoyancy
        dhdx = [(h(2:end)-h(1:end-1))./(x(2:end)-x(1:end-1)) 0]; % surface slope (unitless)

        % find the calving front location (based on Benn et al., 2007 & Nick et al., 2010)
        Rxx = 2*nthroot((dUdx./(E*A(1:length(dUdx)))),n); %resistive stress (Pa)
        crev = (Rxx./(rho_i.*g))+((rho_sw./rho_i).*fwd); %crevasse penetration depth (m)
        c = find(h(1:length(x)-1)-crev(1:length(x)-1)<=0,1,'first'); %calving front located where the inland-most crevasse intersects sea level
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
            c=dsearchn(transpose(x),x0(term));            
        end 

        %calculate the effective pressure (ice overburden pressure minus water
        %pressure) assuming an easy & open connection between the ocean and
        %ice-bed interface
        sl = find(hb<=0,1,'first'); %find where the glacier base first drops below sea level
        N_ground = rho_i*g*H(1:sl); %effective pressure where the bed is above sea level (Pa)
        N_marine = rho_i*g*H(sl+1:length(x))+(rho_sw*g*hb(sl+1:length(x))); %effective pressure where the bed is below sea level (Pa)
        N = [N_ground N_marine];
        N(N<0)=1; %cannot have negative values

        % Solve for new beta
        [beta,U,vm,T,dUdx] = betaSolve(x,h,rho_i,g,m,n,A,H,W,dx0(z),c,ice_end,N,i,U0,x0,t(i),beta,E,rho_sw);

        % Plot geometries
            figure(1) 
            plot(x/10^3,h,'-','color',col(z,:),'linewidth',2,'displayname',num2str(dx0(z))); 
            plot([x(c)/10^3 x(c)/10^3],[h(c) h(c)-H(c)],'-','color',col(z,:),'linewidth',2,'HandleVisibility','off'); 
            plot(x(gl:end)/10^3,h(gl:end)-H(gl:end),'-','color',col(z,:),'linewidth',2,'HandleVisibility','off');        
        % Plot velocity
            figure(2)
            plot(x/10^3,U.*3.1536e7,'-','Color',col(z,:),'linewidth',2,'displayname',num2str(dx0(z))); hold on;           
        % Plot resulting beta
            figure(3);
            subplot(1,3,1);
            plot(x(1:c-1)./10^3,beta(i,1:c-1),'-','Color',col(z,:),'linewidth',2,'displayname',num2str(dx0(z))); hold on;
        
            % Solve for the best polynomial fit of beta using cross-validation
            if poly_solve==1
                deg = [0:5]; % polynomial degrees
                nMC=length(x); % number of Monte-Carlo simulations
                pTrain=0.9; % percentage of points to use as training data 
                replace=0; % do not replace values in MC simulations

                % Input into bestPolyFit function
                [RMSE,bestPoly] = bestPolyFit(x(1:end-1),beta(i,1:end-1),pTrain,deg,nMC,replace);
                betaFit = polyfit(x(1:end-1),beta(i,1:end-1),bestPoly);
                if length(betaFit)>1 && betaFit(2)<0
                    betaFit(2)=0;
                end
                betaFit = polyval(betaFit,x0);

                % Plot polynomial fit
                figure(3);
                subplot(1,3,2);
                plot(x0./10^3,betaFit,'color',col(z,:),'linewidth',2,...
                    'displayname',num2str(dx0(z)));
        
                mu_RMSE(z) = mean(RMSE);
                
                % Plot RMSE of polynomial fit
                figure(3);
                subplot(1,3,3); 
                plot(dx0(z),mean(RMSE),'*','color',col(z,:),'markersize',15,...
                    'linewidth',2,'displayname',num2str(dx0(z)));
                
            end 

    end 

end