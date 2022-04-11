%% Script to save model initialization parameters for steady-state model simulation
% Rainey Aberle
% Spring 2020
%
% Required inputs:
%   1. Centerline coordinates, x0
%   2. Surface elevation, h0(x)
%   3. Bed elevation, b0(x)
%   4. Width, W0(x)
%   5. Surface speed, U0(x)
%   6. Rate factor, A0(x)
%   7. Surface mass balance, SMB0(x)
%   8. Maximum submarine melt rate, SMR0
%   9. Tributary ice volume flux, Q0(x)
%   10. Calving front location (index of x), c0
% -------------------------------------------------------------------------

clear all; close all;

% Define homepath
homepath = '/Users/raineyaberle/Research/MS/CraneGlacier_flowlinemodeling/';
cd([homepath,'inputs-outputs/']);

% Modify settings
save_initial = 1;   % = 1 to save initialization file
regrid = 1;         % = 1 to regrid to dx resolution
L = 100e3;          % length of model domain [m]
dx = 200;           % model grid spacing [m]

% -------------------------------------------------------------------------
% 1. Centerline coordinates, x0
% -------------------------------------------------------------------------
cl.x = load('Crane_centerline.mat').x; cl.y = load('Crane_centerline.mat').y; 
if size(cl.x)==[186 1]
    cl.x=cl.x';
end
if size(cl.y)==[186 1]
    cl.y=cl.y';
end
% define x as distance along centerline
x0 = zeros(1,length(cl.y));
for i=2:length(x0)
    x0(i) = sqrt((cl.x(i)-cl.x(i-1))^2+(cl.y(i)-cl.y(i-1))^2)+x0(i-1);
end

% -------------------------------------------------------------------------
% 2. Surface elevation, h0(x)
% -------------------------------------------------------------------------
h0 = load('surfaceElevationObs_1996-2018.mat').h(14).h_centerline;
if size(h0)==[186 1]
    h0=h0';
end
h0(isnan(h0))=0; % replace NaNs with zeros
h0 = movmedian(h0, 5);

% -------------------------------------------------------------------------
% 3. glacier bed elevation, b0(x)
% -------------------------------------------------------------------------
b0 = load('bedElevation_widthAveraged.mat').b_adj;
b0(1:2) = b0(3);
% smooth slightly
b0=movmean(b0,20);
if size(b0)==[186 1]
    b0=b0';
end
% adjust bed elevation to account for grounding line location
% initial grounding line ~5km inland of fjord_end (Rebesco et al., 2006)
fjord_end = shaperead([homepath,'data/terminus/fjord_end.shp']);
Ifjord_end = find(cl.x > polyxpoly(cl.x, cl.y, fjord_end.X, fjord_end.Y), 1, 'first');
gl0 = find(x0 >= x0(Ifjord_end)-5e3, 1, 'first');
% calculate required thickness to be perfectly grounded at gl0 using
% surface h0
% rho_i*Hi = rho_sw*H_sw --> rho_i*H = rho_sw*(-b) --> H = (rho_sw/rho_i * h) / (rho_sw/rho_i - 1)
% rho_sw = 1028; % kg/m^3
% rho_i = 917; % kg/m^3
% Hf0 = (rho_sw/rho_i * h0(gl0)) / (rho_sw/rho_i - 1);
% b0_gl0 = h0(gl0) - Hf0;

% increase elevation near terminus
% b0(136:end) = b0(136:end)+100;
% % smooth the transition
% b0(125:135) = interp1([x0(124);x0(136)],[b0(124);b0(136)],x0(125:135));

% -------------------------------------------------------------------------
% 4. glacier width, W0(x)
% -------------------------------------------------------------------------
W0 = load('calculatedWidth.mat').width.W;
if size(W0)==[186 1]
    W0=W0';
end

% -------------------------------------------------------------------------
% 5. surface speed, U0(x)
% -------------------------------------------------------------------------
U0 = double(load('surfaceSpeeds_widthAveraged_1994-2018').U(22).U_width_ave);
if size(U0)==[186 1]
    U0=U0';
end
U0 = movmedian(U0, 5);

% -------------------------------------------------------------------------
% 6. rate factor, A0(x)
% -------------------------------------------------------------------------
A0 = load('adjustedRateFactor.mat').A_adj;
if size(A0)==[186 1]
    A0=A0';
end

% -------------------------------------------------------------------------
% 8. surface mass balance, SMB0(x)
% -------------------------------------------------------------------------
% Use the mean value at each point along the centerline for 2009-2019 
SMB0 = load('downscaledClimateVariables_2009-2019.mat').SMB.downscaled_average_linear./3.1536e7; % m/s
%RO0 = load('downscaledClimateVariables_2009-2019.mat').RO.downscaled_average_linear./3.1536e7; % m/s
RO0 = load('downscaledClimateVariables_2009-2019.mat').SM.downscaled_average_linear.*0.05./3.1536e7; % m/s
    % From Vaughan (2006): ﻿the upper estimates for runoff in 2000 and 2050 are equivalent 
    % to only 11% and 27% of the yearly total snow accumulation for the Antarctic Peninsula, 
    % respectively (given by areas H to J in Vaughan et al., 1999)

% -------------------------------------------------------------------------
% 9. submarine melting rate, SMR0(x)
% -------------------------------------------------------------------------
%   Dryak and Enderlin (2020), Crane iceberg melt rates:
%       2013-2014: 0.70 cm/d = 8.1e-8 m/s
%       2014-2015: 0.51 cm/d = 5.29e-8 m/s
%       2015-2016: 0.46 cm/d = 4.77e-8 m/s
%       2016-2017: 0.08 cm/d = 0.823e-8 m/s
%       Mean melt rate (2013-17) = 4.75e-8 m/s = 1.5 m/yr
%       Max melt rate (2013-17) = 8.1e-8 m/s = 2.55 m/yr
%   Adusumilli et al. (2018):
%       Larsen C basal melt rate (1994-2016) = 0.5 +/- 1.4 m/a = 1.59e-8 m/s
%       Larsen C net mass balance (1994-2016)= -0.4 +/- 1.3 m/a = 1.27e-8 m/s
SMR0 = -1.5/3.1536e7; % m/s SMR - mean found at Crane

% -------------------------------------------------------------------------
% 10. Tributary ice volume flux, Q0(x)
% -------------------------------------------------------------------------
Q = load('tributaryFluxes.mat').tribFlux.Q; % m/a
Q_err = load('tributaryFluxes.mat').tribFlux.Q_err; % m/a
x_Q = load('tributaryFluxes.mat').tribFlux.x; % m along centerline
Q0 = interp1(x_Q,Q,x0)./3.1536e7; % m/s
Q0_err = interp1(x_Q,Q_err,x0)./3.1536e7; % m/s

% -------------------------------------------------------------------------
% 10. Calving front location (index of x), c0
% -------------------------------------------------------------------------
% use 2002 observed terminus position 
term = load([homepath,'inputs-outputs/terminusPositions_2002-2019.mat']).term;
c0 = dsearchn(x0', term.x(1));

% end of fjord location
f0 = find(cl.x >= -2402990, 1, 'first');

% -------------------------------------------------------------------------
% Regrid spatial variables and extend to length of model domain
% -------------------------------------------------------------------------
if regrid
    xi = 0:dx:L; % new grid vector
    h0 = interp1(x0,h0,xi); h0(find(isnan(h0),1,'first'):end) = h0(find(isnan(h0),1,'first')-1);
    b0 = interp1(x0,b0,xi); b0(find(isnan(b0),1,'first'):end) = b0(find(isnan(b0),1,'first')-1);
    W0 = interp1(x0,W0,xi); W0(find(isnan(W0),1,'first'):end) = W0(find(isnan(W0),1,'first')-1);
    U0 = interp1(x0,U0,xi); U0(find(isnan(U0),1,'first'):end) = U0(find(isnan(U0),1,'first')-1);
    A0 = interp1(x0,A0,xi); A0(find(isnan(A0),1,'first'):end) = A0(find(isnan(A0),1,'first')-1);
    SMB0 = interp1(x0,SMB0,xi); SMB0(find(isnan(SMB0),1,'first'):end) = SMB0(find(isnan(SMB0),1,'first')-1);
    RO0 = interp1(x0,RO0,xi); RO0(find(isnan(RO0),1,'first'):end) = RO0(find(isnan(RO0),1,'first')-1);
    Q0 = interp1(x0,Q0,xi); Q0(find(isnan(Q0),1,'first'):end) = Q0(find(isnan(Q0),1,'first')-1);
    c0 = dsearchn(xi',x0(c0));
    f0 = dsearchn(xi',x0(f0));
    x0=xi;
end

% -------------------------------------------------------------------------
% Save parameters to file
% -------------------------------------------------------------------------
if save_initial
    save('modelInitialization_preCollapse.mat','h0','b0','W0','A0','U0',...
        'x0','SMB0','RO0','SMR0','Q0','c0','-append');
    disp('model initialization parameters saved to file');
end

%% Estimate pre-collapse glacier ice backstress between 2002 terminus position and the end of the fjord
% lateral resistance = 2*H/W * (5*U/E*A*W)^(1/n)

% floating thickness from c0:f0
% rho_i*H_i = rho_sw*H_sw --> H = -rho_sw/(rho_i - rho_sw) * h
rho_i = 917; % kg/m^3
rho_sw = 1000; % kg/m^3
H = -rho_sw/(rho_i-rho_sw) * h0(c0:f0);

% lateral resistance
n = 3;
Rxy = 2*H./W0(c0:f0) .* nthroot(5.*U0(c0:f0) ./ (A0(c0:f0).*W0(c0:f0)), n);
Rxy_sum = sum(Rxy); % Pa
% longitudinal stress
dUdx(c0) = (U0(c0+1)-U0(c0))./(x0(c0+1)-x0(c0)); % forward difference
dUdx(c0+1:f0-1) = (U0(c0+2:f0)-U0(c0:f0-2))./(x0(c0+2:f0)-x0(c0:f0-2)); % central difference
dUdx(f0) = (U0(f0)-U0(f0-1))/(x0(f0)-x0(f0-1)); % backward difference at c
vm = (A0(c0:f0).^(-1/n)).*(abs(dUdx(c0:f0))).^((1-n)/n);
Rxx = (-2/(dx^2)).*(H.*vm);
Rxx(1)=Rxx(2); % remove infinite value
Rxx_sum = sum(Rxx); % Pa
disp(['Sum Rxy = ',num2str(Rxy_sum/10^3),' kPa']);
disp(['Sum Rxx = ',num2str(Rxx_sum/10^3),' kPa']);

% plot
figure(1); clf;
hold on; set(gca,'fontsize',12,'linewidth',1);
plot(x0(c0:f0)/10^3, Rxy/10^3, '-b', 'linewidth',2, 'displayname', 'R_{xy}');
plot(x0(c0:f0)/10^3, Rxx/10^3, '-m', 'linewidth', 2, 'displayname', 'R_{xx}');
grid on; legend;
xlabel('distance along centerline [km]');
ylabel('lateral resistance [kPa]');

figure(2); clf
subplot(1,2,1); 
hold on; set(gca,'fontsize',12,'linewidth',1);
plot(x0/10^3, U0.*3.1536e7, '-k','linewidth',2);
xlabel('distance along centerline [km]');
ylabel('speed [m/y]');
grid on
subplot(1,2,2);
hold on; set(gca,'fontsize',12,'linewidth',1);
plot(x0/10^3, h0, '-b','linewidth',2);
plot(x0(c0:f0)/10^3, h0(c0:f0)-H,'-c','linewidth',2);
plot(x0/10^3, b0,'-k', 'linewidth',2);
xlabel('distance along centerline [km]');
ylabel('elevation [m]');
grid on
