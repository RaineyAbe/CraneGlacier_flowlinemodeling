%% Calculate cumulative strain along Crane Glacier centerline
% RKA 12/2020
% Script to calculate the cumulative strain along the Crane Glacier
% centerline which will be used to adjust the rate factor 

clear all; close all;

save_A_adj = 1; % = 1 to save adjusted rate factor 

% Define homepath in directory
homepath = '/Users/raineyaberle/Desktop/Research/CraneGlacier_flowlinemodeling/';
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
set(gcf,'Position',[200 100 1000 400]);
    subplot(1,2,1); 
    hold on; grid on; set(gca,'fontsize',14,'linewidth',2); title('Cumulative Strain Rate');
    ylabel('$$ \Sigma ( \dot{\eta} ) $$','Interpreter','latex','fontname','Arial','fontsize',18);
    legend('Location','best');

% 1. Calculate cumulative strain rate over time (starting in 2007)
%   eta_dot = cumsum( (\delta x / U)*(dUdx) )

    Iu = [2 4 6 8 9 14:19]; % index of annual velocities to use (2007-2017)
    col = parula(length(Iu)+1); % color scheme for plotting

    eta_dot_cum = zeros(length(Iu),length(cl.x)); % cumulative strain rate (unitless)
    eta_dot = zeros(length(Iu),length(cl.x)); % strain rate (1/s)
    dx = [0 x(2:end)-x(1:end-1)]; % (m)
    
    for i=1:length(Iu)
        
        % strain rate (1/yr)
        eta_dot(i,1) = (U(Iu(i)).speed_linearextrap(2)-U(Iu(i)).speed_linearextrap(1))'...
            ./(x(2)-x(1)).*3.1536e7; % forward difference
        eta_dot(i,2:end-1) = (U(Iu(i)).speed_linearextrap(3:end)-U(Iu(i)).speed_linearextrap(1:end-2))'...
            ./(x(3:end)-x(1:end-2)).*3.1536e7; % central difference
        eta_dot(i,end) = (U(Iu(i)).speed_linearextrap(end-1)-U(Iu(i)).speed_linearextrap(end))'...
            ./(x(end-1)-x(end)).*3.1536e7; % backward difference
        % cumulative strain (1/yr)
        eta_dot_cum(i,:) = cumsum(eta_dot(i,:).*dx./(U(Iu(i)).speed_linearextrap'.*3.1536e7));

        % Plot resulting strain rate
        figure(1);
        subplot(1,2,1);
        plot(x./10^3,eta_dot_cum(i,:),'color',col(i,:),'displayname',num2str(round(U(Iu(i)).date)),'linewidth',2);

    end

% 2. Calculate relative weights of cumulative strain rate at each year
    
    coeff = eta_dot_cum;
    
% 3. Adjust rate factor using coefficients
    
    % Set up subplot
    figure(1);
    subplot(1,2,2); hold on; grid on; set(gca,'fontsize',14,'linewidth',2);
    ylabel('A_{adjusted} (Pa^{-n} s^{-1})');
    title('Adjusted Rate Factor');
    
    for i=1:length(Iu)
        A_adj(i).A_adj = A'.*(1+coeff(i,:));
        A_adj(i).date = U(Iu(i)).date;
        A_adj(i).x = x;
        A_adj(i).A_adj(135:end) = A_adj(i).A_adj(134); % use the last point for the rest of the centerline
        plot(x./10^3,A_adj(i).A_adj,'color',col(i,:),'displayname',num2str(round(U(Iu(i)).date)),'linewidth',2);        
    end
    
% 4. Save adjusted rate factor
if save_A_adj
    cd([homepath,'inputs-outputs']);
    save('Crane_adjustedAnnualRateFactor_2009-2019.mat','A_adj');
    disp(['Adjusted rate factor saved in: ',pwd]);
end 
    
    


