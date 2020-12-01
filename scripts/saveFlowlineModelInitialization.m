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
x0 = load('Crane_centerline.mat').x; y0 = load('Crane_centerline.mat').y; 
if size(x0)==[186 1]
    x0=x0';
end
if size(y0)==[186 1]
    y0=y0';
end
    % define x as distance along centerline
    x = zeros(1,length(y0));
    for i=2:length(x)
        x(i) = sqrt((x0(i)-x0(i-1))^2+(y0(i)-y0(i-1))^2)+x(i-1);
    end

% 2. ice surface elevation
h0 = load('Crane_OIBObservations_2009-2018.mat').Crane_OIB(1).surf;
if size(h0)==[186 1]
    h0=h0';
end
    % Use linear interp to fill in gaps
    h0 = feval(fit(x(~isnan(h0))',h0(~isnan(h0))','linearinterp'),x);

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
U0 = load('Crane_CenterlineSpeeds_2007-2017.mat').U(1).speed;
if size(U0)==[186 1]
    U0=U0';
end

% 6. rate factor, A
A0 = load('Crane_RateFactorA.mat').A;
if size(A0)==[186 1]
    A0=A0';
end

% 7. basal roughness factor, beta
beta = load('Crane_CalculatedBeta.mat').beta.beta;
xbeta = load('Crane_CalculatedBeta.mat').beta.x;
beta0 = interp1(xbeta,beta,x);

% 8. surface mass balance w/ uncertainty, smb and smb_err
smb = load('Crane_downscaledSMB_2009-2016.mat').SMB(1).smb_linear; % m/a
smb_err = load('Crane_downscaledSMB_2009-2016.mat').SMB(1).sigma_smb; % m/a
smb0 = feval(fit(x0(1:length(smb))',smb,'poly1'),x0)./3.1536e7; % m/s
smb0_err = [smb_err smb_err(end).*ones(1,length(x0)-length(smb_err))]./3.1536e7; % m/s

% 9. tributary ice volume flux
Q = load('TributaryFlux.mat').tribFlux.Q; % m/a
Q_err = load('TributaryFlux.mat').tribFlux.Q_err; % m/a
Q_x = load('TributaryFlux.mat').tribFlux.x; % m along centerline
Q0 = interp1(Q_x,Q,x)./3.1536e7; % m/s
Q0_err = interp1(Q_x,Q_err,x)./3.1536e7; % m/s

% Save resulting variables
x0=x;
save('Crane_flowline_initialization.mat','h0','hb0','W0','A0','U0',...
    'beta0','x0','smb0','smb0_err','Q0','Q0_err');
disp(['initialization variable saved in: ',pwd]);


