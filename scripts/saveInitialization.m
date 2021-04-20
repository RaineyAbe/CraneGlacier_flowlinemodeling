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

% Define homepath
homepath = '/Users/raineyaberle/Desktop/Research/CraneGlacier_flowlinemodeling/';
cd([homepath,'inputs-outputs/']);

save_initial = 1; % = 1 to save initialization file
regrid = 0;       % = 1 to regrid to 200m resolution

L = 70e3; % length of model domain

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
h0 = load('Crane_surfaceElevationObs.mat').h(1).surface;
if size(h0)==[186 1]
    h0=h0';
end
h0(isnan(h0))=0; % replace NaNs with zeros

% 3. glacier bed elevation
hb0 = load('Crane_observedBed_Tate.mat').hb.hb0;
if size(hb0)==[186 1]
    hb0=hb0';
end

% 4. glacier width
W0 = load('Crane_calculatedWidth.mat').width.W;
if size(W0)==[186 1]
    W0=W0';
end

% 5. ice surface speed
U = load('Crane_centerlineSpeeds_2007-2017.mat').U(13).speed;
if size(U)==[186 1]
    U=U';
end
    % use first point in 2013 speed profile
    % (first real data point near ice divide - speed does not seem to
    % change much over time in this region of the glacier)
    U(1) = load('Crane_centerlineSpeeds_2007-2017.mat').U(4).speed(1);
    % fit smoothing spline to fill in gaps
    U0 = feval(fit(x0(~isnan(U))',U(~isnan(U))','linearinterp'),x0);
    if size(U0)==[186,1]
        U0=U0';
    end

% 6. rate factor, A & adjusted rate factor
A0 = load('Crane_adjustedAnnualRateFactor_2009-2019.mat').A_adj(1).A_adj;
%A0 = polyval(polyfit(x0,A0,1),x0);
%A = load('Crane_rateFactorA.mat').A;
%A0 = polyval(polyfit(x0,A,1),x0);
if size(A0)==[186 1]
    A0=A0';
end

% 7. basal roughness factor, beta
beta0 = load('optimalBeta.mat').optBeta.beta;
beta0x = load('optimalBeta.mat').optBeta.x;
%beta0 = interp1(beta0x,beta0,x0);

% 8. surface mass balance w/ uncertainty, smb and smb_err
smb = load('Crane_downscaledSMB_2009-2019.mat').SMB(11).smb_interp; % m/a
smb_err = load('Crane_downscaledSMB_2009-2019.mat').SMB(1).sigma_smb; % m/a
smb0 = [smb' smb(end).*ones(1,length(x0)-length(smb))]./3.1536e7; % m/s
% replace NaN values with the last SMB value near the terminus
smb0(isnan(smb0)) = smb0(find(isnan(smb0),1,'first')-1); 
smb0_err = [smb_err smb_err(end).*ones(1,length(x0)-length(smb_err))]./3.1536e7; % m/s
% increase smb slope to incorporate potential runoff
for i=1:length(smb0)
    if i<135
        smb0(i) = smb0(i)-15/3.1536e7;%-5/3.1536e7/x0(135)*x0(i);
    else
        smb0(i) = smb0(i)-15/3.1536e7;        
    end
end

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
Q = load('tributaryFluxes.mat').tribFlux.Q; % m/a
Q_err = load('tributaryFluxes.mat').tribFlux.Q_err; % m/a
x_Q = load('tributaryFluxes.mat').tribFlux.x; % m along centerline
Q0 = interp1(x_Q,Q,x0)./3.1536e7; % m/s
Q0_err = interp1(x_Q,Q_err,x0)./3.1536e7; % m/s

% 10. calving front location
termx = load('LarsenB_centerline.mat').centerline.termx;
termy = load('LarsenB_centerline.mat').centerline.termy;
termdate = load('LarsenB_centerline.mat').centerline.termdate;
c0 = dsearchn([x_cl' y_cl'],[termx(5),termy(5)]);

% regrid spatial variables and extend to length of model domain
if regrid
    xi = 0:200:L; % new grid vector
    h0 = interp1(x0,h0,xi); h0(find(isnan(h0),1,'first'):end) = h0(find(isnan(h0),1,'first')-1);
    hb0 = interp1(x0,hb0,xi); hb0(find(isnan(hb0),1,'first'):end) = hb0(find(isnan(hb0),1,'first')-1);
    W0 = interp1(x0,W0,xi); W0(find(isnan(W0),1,'first'):end) = W0(find(isnan(W0),1,'first')-1);
    U0 = interp1(x0,U0,xi); U0(find(isnan(U0),1,'first'):end) = U0(find(isnan(U0),1,'first')-1);
    A0 = interp1(x0,A0,xi); A0(find(isnan(A0),1,'first'):end) = A0(find(isnan(A0),1,'first')-1);
    beta0 = interp1(beta0x,beta0,xi); beta0(find(isnan(beta0),1,'first'):end) = beta0(find(isnan(beta0),1,'first')-1);
    smb0 = interp1(x0,smb0,xi); smb0(find(isnan(smb0),1,'first'):end) = smb0(find(isnan(smb0),1,'first')-1);
    Q0 = interp1(x0,Q0,xi); Q0(find(isnan(Q0),1,'first'):end) = Q0(find(isnan(Q0),1,'first')-1);
    c0 = dsearchn(xi',x0(c0));
    x0=xi;
else
    beta0 = interp1(beta0x,beta0,x0);
end

% Save resulting variables
if save_initial
    save('Crane_flowlinemodelInitialization.mat','h0','hb0','W0','A0','beta0','U0',...
        'x0','smb0','smr0','Q0','c0');
    disp(['initialization variable saved in: ',pwd]);
end


