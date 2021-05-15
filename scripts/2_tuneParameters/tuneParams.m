%% Script to tune the basal roughness factor and the fresh water depth in crevasses
% using observations of glacier speed, geometry, and terminus position.
% Adapted from Ellyn Enderlin's flowline model (Enderlin et al., 2013)
% Rainey Aberle
% 2020-2021

%% define constants and spatial variables
  
clear all; close all;
warning off; % turn off warnings (velocity coefficient matrix is close to singular)
    
% define home path in directory
homepath = '/Users/raineyaberle/Desktop/Research/CraneGlacier_flowlinemodeling/';
addpath([homepath,'../matlabFunctions/cmocean_v2.0/cmocean/']);
addpath([homepath,'scripts/tuneParmeters/']);
addpath([homepath,'inputs-outputs/']);

% Load Crane Glacier initialization variables
load('Crane_flowlineModelInitialization.mat');
cl.X = load('Crane_centerline.mat').x; cl.Y = load('Crane_centerline.mat').y;
cl.x = zeros(1,length(cl.X));
for i=2:length(cl.X)
    cl.x(i) = sqrt((cl.X(i)-cl.X(i-1))^2+(cl.Y(i)-cl.Y(i-1))^2)+cl.x(i-1);
end

% 2019 observed terminus position
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
fwd = 20; % fresh water depth in crevasses (m)

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

% extend variables past where there are NaNs (set to last value)
A0(find(isnan(A0),1,'first'):end) = A0(find(isnan(A0),1,'first')-1);
U0(find(isnan(U0),1,'first'):end) = U0(find(isnan(U0),1,'first')-1);
dUdx0(find(isnan(dUdx0),1,'first'):end) = dUdx0(find(isnan(dUdx0),1,'first')-1);
hb0(find(isnan(hb0),1,'first'):end) = hb0(find(isnan(hb0),1,'first')-1);
W0(find(isnan(W0),1,'first'):end) = W0(find(isnan(W0),1,'first')-1);    

%% solve for basal roughness factor, beta

save_beta = 1;    % = 1 to save parameter solutions
save_figure = 1;   % = 1 to save figures

% initialize variables
x=x0; H=H0; hb=hb0; U=U0; dUdx=dUdx0; A=A0; dx=dx0; U0(c0+1:end)=NaN; gl=gl0;

% solve for beta by gradient descent 
% start timer
tic
fh = @(beta)betaSolve(beta,H,x,U,hb,n,E,m,dx0,rho_i,g,rho_sw,sigma_b,dUdx,c0,x0,hb0,W0,A0,U0,h0);
beta0 = 2*ones(1,length(x)); % initial guess for beta
[beta,~] = fmincon(fh,beta0,[],[],[],[],0,10,[],optimoptions('fmincon','StepTolerance',1e-4,'MaxFunctionEvaluations',1e6));
% grab velocity from beta solution
beta(gl:end)=0;
[J,Un,xn] = betaSolve(beta,H,x,U,hb,n,E,m,dx0,rho_i,g,rho_sw,sigma_b,dUdx,c0,x0,hb0,W0,A0,U0,h0);
% stop timer
toc

% plot final results
figure(1); clf; hold on;
set(gcf,'position',[440 300 700 450]);
set(gca,'fontsize',14,'linewidth',2); grid on; legend('Location','northwest');
xlabel('distance along centerline (km)');ylabel('\beta (s^{1/m} m^{-1/m})');
plot(x(1:c0)/10^3,beta(1:c0),'linewidth',2,'displayname','\beta');
yyaxis right; ylabel('U (m a^{-1})');
plot(x0/10^3,U0*3.1536e7,'--k','linewidth',2,'displayname','U_{obs}');
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
    disp('beta solution saved')
end

%% fine tune with A

save_A = 1;    % = 1 to save parameter solutions
save_figure = 1;   % = 1 to save figures

% initialize variables
x=x0; H=H0; hb=hb0; U=U0; dUdx=dUdx0; dx=dx0; U0(c0+1:end)=NaN; gl=gl0;

% solve for beta by gradient descent 
% start timer
tic
fh = @(A)betaSolve(beta,H,x,U,hb,n,E,m,dx0,rho_i,g,rho_sw,sigma_b,dUdx,c0,x0,hb0,W0,A,U0,h0);
[A,fval] = fmincon(fh,A0,[],[],[],[],0,1,[],optimoptions('fmincon','StepTolerance',1e-324,'MaxFunctionEvaluations',1e6));
% grab velocity from A solution
[J,Un,xn] = betaSolve(beta,H,x,U,hb,n,E,m,dx0,rho_i,g,rho_sw,sigma_b,dUdx,c0,x0,hb0,W0,A,U0,h0);
% stop timer
toc

% plot final results
figure(2); clf; hold on;
set(gcf,'position',[440 300 700 450]);
set(gca,'fontsize',14,'linewidth',2); grid on; legend('Location','northwest');
xlabel('distance along centerline (km)');ylabel('A');
plot(x(1:c0)/10^3,A(1:c0),'linewidth',2,'displayname','A');
yyaxis right; ylabel('U (m a^{-1})');
plot(x0/10^3,U0*3.1536e7,'--k','linewidth',2,'displayname','U_{obs}');
plot(xn/10^3,Un*3.1536e7,'-k','linewidth',2,'displayname','U_{sol}');

if save_figure
    cd([homepath,'figures/']);
    figure(2);
    saveas(gcf,'ASolution.png','png');
    disp('figure saved');
end

if save_A
    cd([homepath,'inputs-outputs/']);
    save('Crane_flowlineModelInitialization.mat','A','-append');
    disp('A solution saved');
end

%% tune fresh water depth in crevasses, fwd

% define settings
fwdbounds = [5 30]; % range of possible fwd values (m)
method = 'g';       % = 'g' to solve by gradient descent
                    % = 'b' to solve by brute force
plot_timeSteps = 0; % = 1 to plot geometry, speed, grounding line/terminus positions over time
save_fwd = 1;       % = 1 to save resulting fwd
save_figure = 1;    % = 1 to save figure resulting from brute force method

% initialize variables
x=x0; h=h0; H=H0; hb=hb0; U=U0; dUdx=dUdx0; c=c0; gl=gl0;
beta0 = load('Crane_flowlineModelInitialization.mat').beta;
A0 = load('Crane_flowlineModelInitialization.mat').A0;

% solve for fwd which minimizes the cost function J
if strcmp(method,'g')

    % start timer
    tic
    
    % create function handle for the beta solve function
    fh = @(fwd)fwdSolve(dUdx,U,E,A,m,n,rho_i,g,rho_fw,rho_sw,fwd,x0,c0,h,hb,x,H,dx0,hb0,W0,A0,beta0,plot_timeSteps,smr0,smb0,Q0,sigma_b,H_max,U_min,xcf_2019);
    fwd0 = 20; % initial guesses fwd
    % solve for fwd which minimize the cost function J
    [fwd,fval] = fmincon(fh,fwd0,[],[],[],[],fwdbounds(1),fwdbounds(2));   
    
    % stop timer
    toc
    
    disp(['optimal fwd = ',num2str(fwd),' m']);
    
elseif strcmp(method,'b')
    
    % start timer
    tic
    
    % loop through possible fwd values
    fwd0 = fwdbounds(1):0.5:fwdbounds(2);
    J = NaN*ones(1,length(fwd0)); % initialize cost function
    for j=1:length(fwd0)
        fwd=fwd0(j);
        J(j) = fwdSolve(dUdx,U,E,A,m,n,rho_i,g,rho_fw,rho_sw,fwd,x0,c0,h,hb,x,H,dx0,hb0,W0,A0,beta0,plot_timeSteps,smr0,smb0,Q0,sigma_b,H_max,U_min,xcf_2019);
    end
    
    % save fwd with lowest cost
    Ibest = find(J==min(J(:)),1);
    fwd = fwd0(Ibest);
    % stop timer
    toc
    
    % plot 
    figure(3); clf; hold on;
    set(gca,'fontsize',14,'linewidth',2); grid on;
    plot(fwd0,J,'.','markersize',10);
    plot(fwd,J(Ibest),'*r','markersize',15,'linewidth',2);
    xlabel('fwd (m)'); ylabel('misfit (m)');
    title(['Optimal fwd = ',num2str(fwd),' m']);
    
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






