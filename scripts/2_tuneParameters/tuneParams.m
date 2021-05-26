%% Script to tune the basal roughness factor beta, enhancement factor E, and the fresh water depth in crevasses fwd
% using observations of glacier speed, geometry, and terminus position.
% Adapted from Ellyn Enderlin's flowline model (Enderlin et al., 2013)
% Rainey Aberle
% 2020-2021

%% 1. define constants and spatial variables
  
clear all; close all;
warning off; % turn off warnings (velocity coefficient matrix is close to singular)
    
% define home path in directory
homepath = '/Users/raineyaberle/Desktop/Research/CraneGlacier_flowlinemodeling/';
addpath([homepath,'../matlabFunctions/cmocean_v2.0/cmocean/']);
addpath([homepath,'scripts/2_tuneParmeters/']);
addpath([homepath,'inputs-outputs/']);

% Load Crane Glacier initialization variables
load('Crane_flowlineModelInitialization.mat');
cl.X = load('Crane_centerline.mat').x; cl.Y = load('Crane_centerline.mat').y;
cl.x = zeros(1,length(cl.X));
for i=2:length(cl.X)
    cl.x(i) = sqrt((cl.X(i)-cl.X(i-1))^2+(cl.Y(i)-cl.Y(i-1))^2)+cl.x(i-1);
end
F0 = 0; % inner boundary flux (m^3/s)
%c0=210; % initial calving front location

% 2018 observed surface speed and elevation
U_2018 = load('Crane_centerlineSpeedsWidthAveraged_2007-2018.mat').U_widthavg(20).speed;
h_2018 = load('Crane_surfaceElevationObs.mat').h(36).surface;
% interpolate to x0
U_2018 = movmean(double(interp1(cl.x,U_2018,x0)),10);
h_2018 = interp1(cl.x,h_2018,x0);
% 2018 observed calving front position
termx = load('LarsenB_centerline.mat').centerline.termx(58); termy = load('LarsenB_centerline.mat').centerline.termy(58);
xcf_2018 = cl.x(dsearchn([cl.X cl.Y],[termx termy]));
clear termx termy
% 2019 observed calving front position
termx = load('LarsenB_centerline.mat').centerline.termx(61); termy = load('LarsenB_centerline.mat').centerline.termy(61);
xcf_2019 = cl.x(dsearchn([cl.X cl.Y],[termx termy]));
clear termx termy

% grid spacing
dx0 = 200; % m

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
sigma_b = 0; % back pressure (Pa) due to sea ice or sikkusak
fwd = 13; % fresh water depth in crevasses (m)

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
H0(gl0+1:length(x0))=h0(gl0+1:end)*rho_sw/(rho_sw-rho_i); % buoyant thickness using surface
H0(H0>=(h0-hb0))=h0(H0>=(h0-hb0))-hb0(H0>=(h0-hb0)); % thickness can't go beneath bed elevation
H0(c0+1:end) = 0; % zero ice thickness past calving front

% extend variables past where there are NaNs (set to last value)
A0(find(isnan(A0),1,'first'):end) = A0(find(isnan(A0),1,'first')-1);
U0(find(isnan(U0),1,'first'):end) = U0(find(isnan(U0),1,'first')-1);
dUdx0(find(isnan(dUdx0),1,'first'):end) = dUdx0(find(isnan(dUdx0),1,'first')-1);
hb0(find(isnan(hb0),1,'first'):end) = hb0(find(isnan(hb0),1,'first')-1);

W0(find(isnan(W0),1,'first'):end) = W0(find(isnan(W0),1,'first')-1);    

%% 2. solve for beta using 2009 observed surface speed

save_beta = 1;    % = 1 to save parameter solutions
save_figure = 1;   % = 1 to save figures

% initialize variables
x=x0; H=H0; hb=hb0; U=U0; dUdx=dUdx0; U0(c0+1:end)=NaN; A=A0;
fwd=11;

% load stress-coupling length
cd([homepath,'inputs-outputs/']);
ITRxx = load('Crane_SCL_results.mat').ITRxx;
x_scl = load('Crane_SCL_results.mat').x;
scl = zeros(1,length(ITRxx));
for i=2:length(ITRxx)
    scl(i) = x_scl(ITRxx(i))-x_scl(ITRxx(i-1));
end
SCL = nanmean(scl); % m
clear ITRxx x_scl scl

% solve for beta by gradient descent 
cd([homepath,'scripts/2_tuneParameters/']);
% start timer
tic
fh = @(betai)betaSolve(A0,A,betai,H,x,U,hb,n,E,m,dx0,rho_i,g,rho_sw,rho_fw,fwd,sigma_b,dUdx,c0,x0,hb0,W0,U_2018,h_2018,xcf_2018,smr0,smb0,Q0,H_max,U_min,F0);
% solve for beta at the resolution of the stress-coupling length
no_pts = round(x0(end)/SCL); % number of pts at which to solve for a beta value
beta0 = 1.1*ones(1,10); % initial guess for beta
[beta,~] = fmincon(fh,beta0,[],[],[],[],zeros(1,length(beta0)),1.5*ones(1,length(beta0)),[],optimoptions('fmincon','StepTolerance',1e-2,'MaxFunctionEvaluations',1e4));
% grab velocity from beta solution
[~,Un,xn,xcf,beta0x] = betaSolve(A0,A,beta,H,x,U,hb,n,E,m,dx0,rho_i,g,rho_sw,rho_fw,fwd,sigma_b,dUdx,c0,x0,hb0,W0,U_2018,h_2018,xcf_2018,smr0,smb0,Q0,H_max,U_min,F0);
beta = interp1(beta0x,beta,x0,'pchip'); beta(dsearchn(xn',xcf)+1:end)=NaN;
% stop timer
toc

% plot results
figure(1); clf; hold on;
set(gcf,'position',[200 200 900 700]);
yyaxis left; cla; hold on;
set(gca,'fontsize',18,'linewidth',2); grid on; legend('Location','northeast');
xlabel('Distance Along Centerline (km)');ylabel('\beta (s^{1/m} m^{-1/m})');
plot(x0/10^3,beta,'b','linewidth',2,'displayname','\beta');
yyaxis right; cla; ylabel('U (m a^{-1})');
plot(x0/10^3,U_2018*3.1536e7,'--k','linewidth',2,'displayname','U_{obs}');
plot(xn/10^3,Un*3.1536e7,'-k','linewidth',2,'displayname','U_{sol}');

if save_figure
    cd([homepath,'figures/']);
    figure(1);
    saveas(gcf,'betaSolution.png','png');
    disp('figure saved');
end

if save_beta
    cd([homepath,'inputs-outputs/']);
    beta0=beta;
    save('Crane_flowlineModelInitialization.mat','beta0','-append');
    save('Crane_betaSolution.mat','beta','Un','xn','xcf');
    disp('beta solution saved');
end

%% 3. tune fwd using 2018 speed and terminus position

% define settings
fwdbounds = [0 20]; % range of possible fwd values (m)
method = 'g';       % = 'g' to solve by gradient descent
                    % = 'b' to solve by brute force
plot_timeSteps = 0; % = 1 to plot geometry, speed, grounding line/terminus positions over time
save_fwd = 1;       % = 1 to save resulting fwd
save_figure = 0;    % = 1 to save figure resulting from brute force method

% initialize variables
beta0 = load('Crane_flowlineModelInitialization.mat').beta0;
A=A0; x=x0; h=h0; H=H0; hb=hb0; U=U0; dUdx=dUdx0; c=c0; gl=gl0;

% solve for fwd which minimizes the cost function J
cd([homepath,'scripts/2_tuneParameters/']);
if strcmp(method,'g')

    % start timer
    tic
    
    % create function handle for the beta solve function
    fh = @(fwd)fwdSolve(dUdx,U,E,A,m,n,rho_i,g,rho_fw,rho_sw,fwd,x0,c0,h,hb,x,H,dx0,hb0,W0,A0,beta0,plot_timeSteps,smr0,smb0,Q0,sigma_b,H_max,U_min,xcf_2018,U_2018);
    fwd0 = 6; % initial guesses fwd
    % solve for fwd which minimize the cost function J
    [fwd,fval] = fmincon(fh,fwd0,[],[],[],[],fwdbounds(1),fwdbounds(2));   
    
    % stop timer
    toc
    
    disp(['optimal fwd = ',num2str(fwd),' m']);
    
elseif strcmp(method,'b')
    
    % start timer
    tic
    
    % loop through possible fwd values
    fwd0 = fwdbounds(1):fwdbounds(2);
    J = NaN*ones(1,length(fwd0)); % initialize cost function
    for j=1:length(fwd0)
        fwd=fwd0(j);
        [J(j),~,~,~] = fwdSolve(dUdx,U,E,A,m,n,rho_i,g,rho_fw,rho_sw,fwd,x0,c0,h,hb,x,H,dx0,hb0,W0,A0,beta0,plot_timeSteps,smr0,smb0,Q0,sigma_b,H_max,U_min,xcf_2018,U_2018);
    end
    
    % save fwd with lowest cost
    Ibest = find(J==min(J(:)),1);
    fwd = fwd0(Ibest);
    % stop timer
    toc

     % plot results
    [~,x,U,xcf] = fwdSolve(dUdx,U,E,A,m,n,rho_i,g,rho_fw,rho_sw,fwd,x0,c0,h,hb,x,H,dx0,hb0,W0,A0,beta0,plot_timeSteps,smr0,smb0,Q0,sigma_b,H_max,U_min,xcf_2018);
    figure(3); clf; hold on;
    set(gcf,'position',[200 300 1000 500]);
    subplot(1,2,1); hold on;
    set(gca,'fontsize',18,'linewidth',2); grid on;
    plot(fwd0,J,'.','markersize',15);
    plot(fwd,J(Ibest),'*r','markersize',15,'linewidth',2);
    xlabel('FWD (m)'); ylabel('Cost');
    subplot(1,2,2); hold on;
    set(gca,'fontsize',18,'linewidth',2); grid on; legend('location','northwest');
    plot(x/10^3,U*3.1536e7,'-k','linewidth',2,'displayname','U_{sol}');
    plot(x0(1:dsearchn(x0',xcf))/10^3,U_2018(1:dsearchn(x0',xcf))*3.1536e7,'--k','linewidth',2,'DisplayName','U_{obs}');
    xlabel('Distance Along Centelrine (m)'); ylabel('m a^{-1}');    
    
else
    disp('adjust method');
end

if save_fwd
    cd([homepath,'inputs-outputs/']);
    fwd0=fwd;
    save('Crane_flowlineModelInitialization.mat','fwd0','-append');
    disp('fwd saved');
end

if save_figure
    cd([homepath,'figures/']);
    if ishandle(figure(3))
        saveas(gcf,'fwdSolution.png','png');
        disp('figure 3 saved');
    else
        disp('figure 3 does not exist');
    end
end

