%% Script to tune the basal roughness factor beta, enhancement factor E, and the fresh water depth in crevasses fwd
% using observations of glacier speed, geometry, and terminus position.
% Adapted from Ellyn Enderlin's flowline model (Enderlin et al., 2013)
% Rainey Aberle
% 2020-2021

%% 1. define paths to functions and observed conditions
  
clear all; close all;
warning off; % turn off warnings (velocity coefficient matrix is close to singular)
    
% define home path in directory
homepath = '/Users/raineyaberle/Desktop/Research/CraneModeling/CraneGlacier_flowlinemodeling/';
addpath([homepath,'matlabFunctions/cmocean_v2.0/cmocean/']); % path to cmocean color palettes
addpath([homepath,'scripts/2_tuneParameters/']); % path to tuning functions
addpath([homepath,'scripts/']); % path to U_convergence
addpath([homepath,'inputs-outputs/']);

% load centerline coordinates and initial conditions
cl.X = load('Crane_centerline.mat').x; cl.Y = load('Crane_centerline.mat').y;
cl.x = zeros(1,length(cl.X));
for i=2:length(cl.X)
    cl.x(i) = sqrt((cl.X(i)-cl.X(i-1))^2 + (cl.Y(i)-cl.Y(i-1))^2) + cl.x(i-1);
end
x0 = load('flowlineModelInitialization.mat').x0;

% 2018 observed surface speed and elevation
U_2018 = load('centerlineSpeedsWidthAveraged_2007-2018.mat').U_widthavg(20).speed;
h_2018 = load('surfaceElevationObs.mat').h(36).surface;
% interpolate to x0
U_2018 = movmean(double(interp1(cl.x,U_2018,x0)),10);
h_2018 = interp1(cl.x,h_2018,x0);
% 2018 observed calving front position
termx = load('LarsenB_centerline.mat').centerline.termx(58); 
termy = load('LarsenB_centerline.mat').centerline.termy(58);
xcf_2018 = cl.x(dsearchn([cl.X cl.Y],[termx termy]));
clear termx termy
% 2019 observed calving front position
termx = load('LarsenB_centerline.mat').centerline.termx(61); termy = load('LarsenB_centerline.mat').centerline.termy(61);
xcf_2019 = cl.x(dsearchn([cl.X cl.Y],[termx termy]));
clear termx termy

%% 2. solve for beta using observed surface speed

save_beta = 1;    % = 1 to save parameter solutions
save_figure = 1;   % = 1 to save figures

% initial guess for depth of fresh water in crevasses
DFW0 = 0; % m

% solve for beta by gradient descent 
% start timer
tic
fh = @(betai)betaSolve(homepath,x0,betai,DFW0,U_2018,xcf_2018);
% solve for beta at the resolution of the stress-coupling length
beta0 = [0.5 0.7 1.4 1.2 2.1*ones(1,7)]; % initial guess for beta
[beta,~] = fmincon(fh,beta0,[],[],[],[],0*ones(1,length(beta0)),5*ones(1,length(beta0)),[],optimoptions('fmincon','StepTolerance',1e-2,'MaxFunctionEvaluations',1e4));
% grab velocity from beta solution
[~,Un,xn,xcf,beta0x] = betaSolve(homepath,x0,beta0,DFW0,U_2018,xcf_2018);
beta = interp1(beta0x',beta',x0,'pchip'); beta(dsearchn(xn',xcf)+1:end)=NaN;
% stop timer
toc

% test a beta value
% beta0 = [0.5 0.7 1.4 1.2 2.1*ones(1,7)]; % initial guess for beta
% % grab velocity from beta solution
% [~,Un,xn,xcf,beta0x] = betaSolve(homepath,x0,beta0,DFW0,U_2018,xcf_2018);
% beta = interp1(beta0x',beta0',x0,'pchip'); beta(dsearchn(xn',xcf)+1:end)=NaN;

% plot results
figure(3); clf; hold on;
set(gcf,'position',[200 200 900 700],'defaultAxesColorOrder',[0 0 1; 0 0 0]);
yyaxis left; cla; hold on; ylim([min(beta)-0.1 max(beta)+0.1]);
set(gca,'fontsize',18,'linewidth',2); grid on; legend('Location','northwest');
xlabel('Distance Along Centerline (km)');ylabel('\beta (s^{1/m} m^{-1/m})');
plot(x0/10^3,beta,'b','linewidth',2,'displayname','\beta');
yyaxis right; cla; ylabel('U (m a^{-1})'); ylim([0 900]);
plot(x0(1:dsearchn(x0',xcf_2018))/10^3,U_2018(1:dsearchn(x0',xcf_2018))*3.1536e7,'--k','linewidth',2,'displayname','U_{obs}');
plot(xn/10^3,Un*3.1536e7,'-k','linewidth',2,'displayname','U_{sol}');

if save_figure
    cd([homepath,'figures/']);
    figure(3);
    saveas(gcf,'betaSolution.png','png');
    disp('figure saved');
end

if save_beta
    cd([homepath,'inputs-outputs/']);
    beta0=beta;
    save('flowlineModelInitialization.mat','beta0','-append');
    save('betaSolution.mat','beta','Un','xn','xcf');
    disp('beta solution saved');
end

%% 3. tune DFW using 2019 terminus position

% define settings
DFWBounds = [0 5]; % range of possible fwd values (m)
method = 'g';       % = 'g' to solve by gradient descent
                    % = 'b' to solve by brute force
plotTimeSteps = 0;  % = 1 to plot geometry, speed, grounding line/terminus positions over time
plotMisfits = 1;    % = 1 to plot misfit in model year 2018 with observed conditions
save_DFW = 1;       % = 1 to save resulting DFW
save_figure = 1;    % = 1 to save figure resulting from brute force method

% load beta0 from previous step
beta0 = load('flowlineModelInitialization.mat').beta0; 

% solve for fwd which minimizes the cost function J
if strcmp(method,'g')

    % start timer
    tic
    
    % create function handle for the beta solve function
    fh = @(DFW)DFWSolve(homepath,plotTimeSteps,plotMisfits,beta0,DFW,xcf_2019);
    DFW0 = 0; % initial guess DFW
    % solve for fwd which minimize the cost function J
    [DFW,fval] = fmincon(fh,DFW0,[],[],[],[],DFWBounds(1),DFWBounds(2));   
    
    % stop timer
    toc
    
    disp(['optimal DFW = ',num2str(DFW),' m']);
    
elseif strcmp(method,'b')
    
    % start timer
    tic
    
    % loop through possible fwd values
    DFW0 = DFWBounds(1):DFWBounds(2);
    J = NaN*ones(1,length(DFW0)); % initialize cost function
    for j=1:length(DFW0)
        DFW=DFW0(j);
        [J(j),~,~,~] = DFWSolve(homepath,plotTimeSteps,plotMisfits,beta0,DFW,xcf_2019);
    end
    
    % save fwd with lowest cost
    Ibest = find(J==min(J(:)),1);
    DFW = DFW0(Ibest);
    % stop timer
    toc

     % plot results
    [~,x,U,xcf] = DFWSolve(homepath,plotTimeSteps,plotMisfits,beta0,DFW,xcf_2019);
    figure(3); clf; hold on;
    set(gcf,'position',[200 300 1000 500]);
    subplot(1,2,1); hold on;
    set(gca,'fontsize',18,'linewidth',2); grid on;
    plot(DFW0,J,'.k','markersize',15);
    plot(DFW,J(Ibest),'*r','markersize',15,'linewidth',2);
    xlabel('d_{fw} (m)'); ylabel('Cost');
    subplot(1,2,2); hold on;
    set(gca,'fontsize',18,'linewidth',2); grid on; legend('location','northwest');
    plot(x/10^3,U*3.1536e7,'-k','linewidth',2,'displayname','U_{sol}');
    plot(x0(1:dsearchn(x0',xcf))/10^3,U_2018(1:dsearchn(x0',xcf))*3.1536e7,'--k','linewidth',2,'DisplayName','U_{obs}');
    xlabel('Distance Along Centelrine (m)'); ylabel('m a^{-1}');    
    
else
    disp('adjust method');
end

if save_DFW
    cd([homepath,'inputs-outputs/']);
    DFW0=DFW;
    save('flowlineModelInitialization.mat','DFW0','-append');
    disp('DFW saved');
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

