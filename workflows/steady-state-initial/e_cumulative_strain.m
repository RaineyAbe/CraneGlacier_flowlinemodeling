%% Calculate cumulative strain along Crane Glacier centerline
% RKA 12/2020
% Script to calculate the cumulative strain along the Crane Glacier
% centerline which will be used to adjust the rate factor 

clear all; close all;

save_A_adj = 0; % = 1 to save adjusted rate factor 

% Define homepath in directory
homepath = '/Users/raineyaberle/Research/MS/CraneGlacier_flowlinemodeling/';
cd([homepath,'inputs-outputs']);

% Add path to functions
addpath([homepath,'functions/']);
addpath([homepath,'functions/hugheylab-nestedSortStruct']);

% Load Crane Centerline
cl.x = load('Crane_centerline.mat').x; cl.y = load('Crane_centerline.mat').y;
    % define x as distance along centerline
    x = zeros(1,length(cl.x));
    for i=2:length(cl.x)
        x(i) = sqrt((cl.x(i)-cl.x(i-1))^2+(cl.y(i)-cl.y(i-1))^2)+x(i-1);
    end

% Load observations of speed 2007-2018
U = load('observed_surface_speeds.mat').U; 

% Load calculated rate factor
A = load('modeled_rate_factor.mat').A;

% Set up figure
figure(1); clf;
set(gcf,'Position',[200 100 900 700]);
subplot(2,1,1); hold on;
grid on; set(gca,'fontsize',14,'linewidth',2); 
ylabel('Strain history','fontname','Arial','fontsize',18);
legend('Location','eastoutside');

% 1. Calculate cumulative strain rate over time (starting in 2007)
%   eta_dot = cumsum( (\delta x / U)*(dUdx) )

    Iu = 17:22; % index of annual velocities to use (1994-2017)
    col = parula(length(Iu)+1); % color scheme for plotting

    eta_dot_cum = zeros(length(Iu),length(cl.x)); % cumulative strain rate (unitless)
    eta_dot = zeros(length(Iu),length(cl.x)); % strain rate (1/s)
    dx = [0 x(2:end)-x(1:end-1)]; % (m)
    
    for i=1:length(Iu)
        
        % strain rate (1/yr)
        eta_dot(i,1) = (U(Iu(i)).U_width_ave(2)-U(Iu(i)).U_width_ave(1))'...
            ./(x(2)-x(1)).*3.1536e7; % forward difference
        eta_dot(i,2:end-1) = (U(Iu(i)).U_width_ave(3:end)-U(Iu(i)).U_width_ave(1:end-2))...
            ./(x(3:end)-x(1:end-2)).*3.1536e7; % central difference
        eta_dot(i,end) = (U(Iu(i)).U_width_ave(end-1)-U(Iu(i)).U_width_ave(end))...
            ./(x(end-1)-x(end)).*3.1536e7; % backward difference
        % cumulative strain (1/yr)
        eta_dot_cum(i,:) = cumsum(eta_dot(i,:).*dx./(U(Iu(i)).U_width_ave.*3.1536e7),'omitnan');

        % Plot resulting strain rate
        figure(1);
        plot(x(1:135)./10^3,eta_dot_cum(i,1:135),'color',col(i,:),'displayname',U(Iu(i)).date,'linewidth',1);

    end
    
    % use mean cumulative strain rate to adjust A
    eta_dot_cum_m = nanmean(eta_dot_cum);
    plot(x(1:135)/10^3,eta_dot_cum_m(1:135),'-k','linewidth',2,'displayname','mean');

% 2. Calculate relative weights of cumulative strain rate at each year
    
    eta_dot_cum_m_norm = (eta_dot_cum_m-min(eta_dot_cum_m))/(max(eta_dot_cum_m)-min(eta_dot_cum_m));

% 3. Adjust rate factor using coefficients
    
    % Calculate adjusted rate factor A_adj
    A_adj = A'.*(1+eta_dot_cum_m_norm);
    A_adj(135:end) = A_adj(134); % use the last point for the rest of the centerline
    % Plot
    subplot(2,1,2); hold on; 
    xlabel('Distance along centerline (km)'); ylabel('Rate factor [Pa^{-3} a^{-1}]');
    grid on; legend('Location','eastoutside'); set(gca,'fontsize',14,'linewidth',2);
    plot(x(1:135)/10^3,A_adj(1:135)*3.1536e7,'-k','linewidth',2,'DisplayName','A_{adj}');
    plot(x(1:135)/10^3,A(1:135)*3.1536e7,'--k','linewidth',2,'DisplayName','A');
    
% 4. Save adjusted rate factor and figure
if save_A_adj
    save([homepath,'inputs-outputs/modeled_rate_factor.mat'],'A_adj','eta_dot_cum','-append');
    disp('Adjusted rate factor saved');
    saveas(figure(1),[homepath,'figures/modeled_rate_factor.png'],'png');
    disp('figure 1 saved');
end 
    
    


