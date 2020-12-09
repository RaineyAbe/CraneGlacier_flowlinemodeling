%% Calculate cumulative strain along Crane Glacier centerline
% RKA 12/2020
% Script to calculate the cumulative strain along the Crane Glacier
% centerline which will be used to adjust the rate factor 

clear all; close all;

% Define homepath in directory
homepath = '/Users/raineyaberle/Desktop/Research/CraneGlacier_modeling/';
cd([homepath,'inputs-outputs']);

% Add path to functions
addpath([homepath,'../matlabFunctions/']);
addpath([homepath,'../matlabFunctions/hugheylab-nestedSortStruct']);

% Load Crane Centerline
cl.x = load('Crane_centerline.mat').x; cl.y = load('Crane_centerline.mat').y;
    % define x as distance along centerline
    x = zeros(1,length(cl.x));
    for i=2:length(cl.x)
        x(i) = sqrt((cl.x(i)-cl.x(i-1))^2+(cl.y(i)-cl.y(i-1))^2)+x(i-1);
    end

% Load observations of speed 2007-2017
U = load('Crane_CenterlineSpeeds_2007-2017.mat').U; 

% Load calculated rate factor
A = load('Crane_RateFactorA.mat').A;

% Set up figures
figure(1); clf
    hold on; grid on; set(gca,'fontsize',14,'linewidth',2);
    xlabel('distance along centerline (km)'); title('Cumulative Strain Rate');
    ylabel('$$ \Sigma ( \dot{\eta} ) $$','Interpreter','latex','fontname','Arial','fontsize',18);
    legend('Location','best');

% Calculate cumulative strain rate over time (starting in 2007)
%   eta_dot = cumsum( (\delta x / U)*(dUdx) )

eta_dot_cum = zeros(length(U),length(cl.x)); % cumulative strain rate (unitless)
eta_dot = zeros(length(U),length(cl.x)); % strain rate (1/s)
dx = [0 x(2:end)-x(1:end-1)]; % (m)

Iu = [2 4 6 8 9 14:19]; % index of annual velocities to use (2007-2017)
col = parula(length(Iu)+1); % color scheme for plotting

for i=1:length(Iu)
    
    % strain rate
    eta_dot(i,:) = [0 abs(U(Iu(i)).speed_linearextrap(2:end)-U(Iu(i)).speed_linearextrap(1:end-1))'...
        ./(x(2:end)-x(1:end-1))].*3.1536e7; % (1/yr)
    
    % cumulative strain rate
    if i==1
        % strain rate at year zero
        eta_dot_cum(i,:) = eta_dot(i,:);
    else
        % strain rate at subsequent years
        for j=1:length(x)
            if j==1               
                eta_dot_cum(i,j) = eta_dot_cum(i-1,j)+eta_dot(i,j);
            end
            eta_dot_cum(i,j) = (x(j)/U(Iu(i)).speed_linearextrap(j))*eta_dot(i,j)+eta_dot_cum(i-1,j);
        end
    end
    
    % Plot resulting strain rate
    figure(1);
    plot(x./10^3,eta_dot_cum(i,:),'color',col(i,:),'displayname',num2str(round(U(Iu(i)).date)),'linewidth',2);
    
end
