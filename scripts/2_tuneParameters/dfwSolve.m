function [J,x,U,xcf] = DFWSolve(homepath,plotTimeSteps,plotMisfits,beta0,DFW0,xcf_2018)
% Rainey Aberle, 2021
% Function to tune the basal roughness factor beta using 2018 observed
% speed along the glacier centerline.
%
% INPUTS:   homepath = path to "CraneGlacier_flowlinemodeling" in directory 
%                      (string ending with "/")
%           plotTimeSteps = option to plot the model geometry and speed at
%                           several model time steps (0 = no, 1 = yes).
%           plotMisfits = option to plot the misfits between observed and
%                         modeled surface elevation and speed in model 
%                         year 2018 (0 = no, 1 = yes).
%           x0 = initial centerline grid vector (m along centerline)
%           beta0 = basal roughness factor ()
%           DFW0 = fresh water depth in crevasses (m)
%           xcf_2018 = 2018 observed calving front location (m along centerline)
%
% OUTPUTS:  J = cost, calculated using the observed and modeled speed (unitless)
%           U = modeled speed along the centerline in 2018 (m/s)
%           x = model grid in year 2018 (m along centerline)
%           xcf = modeled calving front position in year 2018 (m along centerline)
%
% NOTE(s):  Requires flowlineModel.m and U_convergence.m

% define time stepping (s)
dt = 0.001*3.1536e7;
t_start = 0*3.1536e7;
t_end = 9*3.1536e7;

% run flowline model
try
    
    [x,U,h,hb,H,gl,c,xcf,dUdx,~,~,~] = flowlineModel(homepath,plotTimeSteps,plotMisfits,dt,t_start,t_end,beta0,DFW0,0,0,0);
    
     % calculate cost of resulting terminus position
    % modified from Morlighem et al., 2010; Larour et al., 2012; Kyrke-Smith et al., 2018
    U_err = 33/3.1536e7; % m/s
    h_err = 22; % m
    xcf_err = 15; % m
    J = abs((xcf-xcf_2018-xcf_err)/xcf_2018);%+... % surface elevation misfit term
    %nanmean(sqrt((U-interp1(x0,U_2018,x)).^2)-U_err)/nanmean(U_2018);%+... % speed misfit term
    
catch
    J=NaN; x=NaN; U=NaN; xcf=NaN;
    disp([num2str(DFW0),' m failed']);
end

end