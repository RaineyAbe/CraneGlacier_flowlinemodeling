%% Create model initialization file
% Rainey aberle
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
%   9. submarine melting rate, smr
%   10. tributary ice volume flux, Q
%   11. calving front location, c

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
    if size(U0)==[186,1]
        U0=U0';
    end

% 6. rate factor, A & adjusted rate factor
A_adj = load('Crane_AdjustedAnnualRateFactor_2009-2019.mat').A_adj;
A = load('Crane_RateFactorA.mat').A;
%A0 = polyval(polyfit(x0,A_adj(1,:),1),x0);
A0 = polyval(polyfit(x0,A,1),x0);
if size(A0)==[186 1]
    A0=A0';
end

% 7. basal roughness factor, beta
beta = load('betabest.mat').betabest.beta;
betax = load('betabest.mat').betabest.x;
beta0 = interp1(betax,beta,x0); %beta0=movmean(beta0,5);

% 8. surface mass balance w/ uncertainty, smb and smb_err
smb = load('Crane_downscaledSMB_2009-2016.mat').SMB(1).smb_interp; % m/a
smb_err = load('Crane_downscaledSMB_2009-2016.mat').SMB(1).sigma_smb; % m/a
smb0 = [smb' smb(end).*ones(1,length(x0)-length(smb))]./3.1536e7; % m/s
smb0_err = [smb_err smb_err(end).*ones(1,length(x0)-length(smb_err))]./3.1536e7; % m/s

% 9. submarine melting rate, smr
%   Dryak and Enderlin (2020), Crane iceberg melt rates:
%       2013-2014: 0.70 cm/d = 8.1e-8 m/s
%       2014-2015: 0.51 cm/d = 5.29e-8 m/s
%       2015-2016: 0.46 cm/d = 4.77e-8 m/s
%       2016-2017: 0.08 cm/d = 0.823e-8 m/s
%       Mean melt rate (2013-17) = 4.75e-8 m/s = 1.5 m/yr
%   Adusumilli et al. (2018):
%       Larsen C basal melt rate (1994-2016) = 0.5+/-1.4 m/a = 1.59e-8 m/s
%       Larsen C net mass balance (1994-2016)= -0.4+/-1.3 m/a = 1.27e-8 m/s
smr0 = -1.5/3.1536e7; % m/s SMR

% 10. tributary ice volume flux
Q = load('TributaryFlux.mat').tribFlux.Q; % m/a
Q_err = load('TributaryFlux.mat').tribFlux.Q_err; % m/a
x_Q = load('TributaryFlux.mat').tribFlux.x; % m along centerline
Q0 = interp1(x_Q,Q,x0)./3.1536e7; % m/s
Q0_err = interp1(x_Q,Q_err,x0)./3.1536e7; % m/s

% 10. calving front location
termx = load('LarsenB_centerline.mat').centerline.termx;
termy = load('LarsenB_centerline.mat').centerline.termy;
termdate = load('LarsenB_centerline.mat').centerline.termdate;
c0 = dsearchn([x_cl' y_cl'],[termx(5),termy(5)]);

% Save resulting variables
save('Crane_flowline_initialization.mat','h0','hb0','W0','A0','U0',...
    'beta0','x0','smb0','smr0','Q0','c0');
disp(['initialization variable saved in: ',pwd]);


