%% Create model initialization file
% RKA
% Fall 2020
% Script to save initialization variables for glacier flowline model using
% 2009 conditions for most variables

% Required inputs along glacier centerline: 
%   1. Crane centerline coordinates
%   2. ice surface elevation
%   3. glacier bed elevation
%   4. glacier width
%   5. ice surface speed
%   6. rate factor, A
%   7. basal roughness factor, beta
%   8. surface mass balance, smb

clear all; close all;

homepath = '/Users/raineyaberle/Desktop/Research/CraneGlacier_modeling/inputs-outputs';
cd(homepath);

% 1. Crane centerline coordinates
x_cl = load('Crane_centerline.mat').x; y_cl = load('Crane_centerline.mat').y; 
if size(x_cl)==[186 1]
    x_cl=x_cl';
end
if size(y_cl)==[186 1]
    y_cl=y_cl';
end
    % define x as distance along centerline
    x0 = zeros(1,length(y_cl));
    for i=2:length(x0)
        x0(i) = sqrt((x_cl(i)-x_cl(i-1))^2+(y_cl(i)-y_cl(i-1))^2)+x0(i-1);
    end

% 2. ice surface elevation
h0 = load('Crane_SurfaceObservations_2009-2018.mat').h(1).surface;
if size(h0)==[186 1]
    h0=h0';
end
h0(isnan(h0))=0; % replace NaNs with zeros

% 3. glacier bed elevation
hb0 = load('Crane_ObservedBed_Tate.mat').hb.hb0;
if size(hb0)==[186 1]
    hb0=hb0';
end

% 4. glacier width
W0 = load('Crane_CalculatedWidth.mat').width.W;
if size(W0)==[186 1]
    W0=W0';
end

% 5. ice surface speed
U = load('Crane_CenterlineSpeeds_2007-2017.mat').U(13).speed;
if size(U)==[186 1]
    U=U';
end
    % use first point in 2013 speed profile
    % (first real data point near ice divide - speed does not seem to
    % change much over time in this region of the glacier)
    U(1) = load('Crane_CenterlineSpeeds_2007-2017.mat').U(4).speed(1);
    % fit smoothing spline to fill in gaps
    U0 = feval(fit(x0(~isnan(U))',U(~isnan(U))','linearinterp'),x0);

% 6. rate factor, A
A0 = load('Crane_RateFactorA.mat').A;
if size(A0)==[186 1]
    A0=A0';
end

% 7. basal roughness factor, beta
beta = load('Crane_CalculatedBeta.mat').beta.beta;
x_beta = load('Crane_CalculatedBeta.mat').beta.x;
beta0 = interp1(x_beta,beta,x0);

% 8. surface mass balance w/ uncertainty, smb and smb_err
smb = load('Crane_downscaledSMB_2009-2016.mat').SMB(1).smb_interp; % m/a
smb_err = load('Crane_downscaledSMB_2009-2016.mat').SMB(1).sigma_smb; % m/a
smb0 = [smb' smb(end).*ones(1,length(x0)-length(smb))]./3.1536e7; % m/s
smb0_err = [smb_err smb_err(end).*ones(1,length(x0)-length(smb_err))]./3.1536e7; % m/s

% 9. tributary ice volume flux
Q = load('TributaryFlux.mat').tribFlux.Q; % m/a
Q_err = load('TributaryFlux.mat').tribFlux.Q_err; % m/a
x_Q = load('TributaryFlux.mat').tribFlux.x; % m along centerline
Q0 = interp1(x_Q,Q,x0)./3.1536e7; % m/s
Q0_err = interp1(x_Q,Q_err,x0)./3.1536e7; % m/s

% Save resulting variables
save('Crane_flowline_initialization.mat','h0','hb0','W0','A0','U0',...
    'beta0','x0','smb0','smb0_err','Q0','Q0_err');
disp(['initialization variable saved in: ',pwd]);


