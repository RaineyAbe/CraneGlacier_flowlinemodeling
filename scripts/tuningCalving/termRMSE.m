function term_rmse = termRMSE(x,c,t,termx_obs,termDate_obs)
% Function to calculate the RMSE for each iteration of tuningCalving
% INPUTS:   x = longitudinal grid along glacier centerline (m)
%           c = terminus / calving front location (index on x)
%           t = model time (seconds since 2009)
%           termx_obs = observed terminus position over time (m along
%               centerline)
%           termDate_obs = dates correponding to observed terminus
%               positions (termx_obs) (decidate)
% OUTPUTS:  termRMSE = root mean square error (RMSE) of the observed 
%               terminus location and the modeled terminus location (c)

% use a polynomial fit for observed terminus position dates
    % make observed terminus conditions column vectors
    if size(termx_obs) == [1 length(termx_obs)]
        termx_obs = termx_obs';
    end
    if size(termDate_obs) == [1 length(termDate_obs)]
        termDate_obs = termDate_obs';
    end    
    termDate_obs_int = 2009:0.1:2020; % a time vector to evaluate the polynomial
    termx_obs_fit = feval(fit(termDate_obs,termx_obs,'poly2'),termDate_obs_int);

% convert time format to decidate (s -> yr + 2009)
t_decidate = t./3.1536e7 + 2009; % decidate

% interpolate observed terminus location at model time
termx_obs_current = interp1(termDate_obs_int,termx_obs_fit,t_decidate); % (m)

% compute the RMSE of the modeled vs. observed terminus position
term_rmse = sqrt((x(c) - termx_obs_current)^2);

end 