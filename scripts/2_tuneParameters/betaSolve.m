function [J,U,x,xcf,beta0x] = betaSolve(homepath,x0,beta0,DFW0,U_2018,xcf_2018)
% Rainey Aberle, 2021
% Function to tune the basal roughness factor beta using 2018 observed
% speed along the glacier centerline.
%
% INPUTS:   homepath = path to "CraneGlacier_flowlinemodeling" in directory 
%                      (string ending with "/")
%           x0 = initial centerline grid vector (m along centerline)
%           beta0 = initial guess for beta at several points along the centerline
%               Does not have to be at the same spatial resolution as x0. 
%               beta0 is interpolated onto the x0 grid.
%           U_2018 = 2018 observed surface speed, width-averaged (m/s)
%           xcf_2018 = 2018 observed calving front location (m along centerline)
%
% OUTPUTS:  J = cost, calculated using the observed and modeled speed (unitless)
%           U = modeled speed along the centerline in 2018 (m/s)
%           x = model grid in year 2018 (m along centerline)
%           xcf = modeled calving front position in year 2018 (m along centerline)
%           beta0x = initial grid for the input beta0 vector (m along centerline)
%
% NOTE(s):  Requires flowlineModel.m and U_convergence.m

% interpolate beta0 to grid spacing
beta0x = 0:round(x0(end)/(length(beta0))):round(x0(end)/(length(beta0)))*(length(beta0)-1); 
beta = interp1(beta0x,beta0,x0,'pchip');
    
% define time stepping (s)
dt = 0.001*3.1536e7;
t_start = 0*3.1536e7;
t_end = 9*3.1536e7;

plotTimeSteps = 0; % = 1 to plot time steps in flowline model
plotMisfits = 0; % = 1 to plot misfits with 2018 observed conditions

% Run flowline model
addpath([homepath,'scripts/']);
[x,U,h,hb,H,gl,c,xcf,dUdx,~,~,~] = flowlineModel(homepath,plotTimeSteps,plotMisfits,dt,t_start,t_end,beta,DFW0,0,0,0);

% calculate cost of parameter solutions
% modified from Morlighem et al., 2010; Larour et al., 2012; Kyrke-Smith et al., 2018
U_err = 33/3.1536e7; % m/s
h_err = 22; % m
K = log(beta); K(~isfinite(K))=0;
J = nanmean(sqrt((U-interp1(x0(1:dsearchn(x0',xcf_2018)),U_2018(1:dsearchn(x0',xcf_2018)),x)).^2)-U_err)/nanmean(U_2018(1:dsearchn(x0',xcf_2018)))+... % speed misfit term
    nanmean(abs(gradient(gradient(K)))); % regularization term to penalize changes in beta gradient        
    %nanmean(sqrt((h-interp1(x0,h_2018,x)).^2)-h_err)./nanmean(h_2018); % surface elevation misfit term

end