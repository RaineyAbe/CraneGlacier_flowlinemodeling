function varargout = wgs2ps(varargin)
%--------------------------------------------------------------------------
% WGS2PS
% -------------------------------------------------------------------------
% SYNTAX:
% [x,y] = WGS2PS(lon,lat)       
% [x,y] = WGS2PS(...,'PropertyName',propertyvalue)
% [x,y,k] = WGS2PS(...)
% [x,y,k,alpha] = WGS2PS(...)
%
% DESCRIPTION:
% Converts from geographic coordinates (longitude, latitude) to polar
% stereographic coordinates using an ellipsoid.  The Polar Stereograpic 
% projection is a conformal projection, i.e. an angle on the sphere remains
% the same in the projection plane. However, the projection distort areas.  
%
% INPUT:  
% lon, lat           Vectors containing geographic longitude and latitude
%                    given in degrees.  
% 'StandardParallel' The standard Parallel for the polar stereographic
%                    projection. Default is 70 degrees north.
% 'StandardMeridian' The Standard Meridian for the polar stereographic
%                    projection. This is the meridian for the positive y
%                    direction. Default is the -45 meridian
% 'Ellipsoid'        Ellipsoid model used in the transformation. Valid
%                    options are:
%                        'WGS-84'
%                        'International'
%                    Default is WGS-84.   
% 'FalseNorthing'    Offset between origin of transformed coordinate system
%                    and origin of projection. Often used to avoid negative
%                    coordinate values. Default is 0.
% 'FalseEasting'     Offset between origin of transformed coordinate system
%                    and origin of projection. Often used to avoid negative
%                    coordinate values. Default is 0.
%
% OUTPUT:
% x,y           Vectors containing transformed polar stereographic 
%               coordinates.
% k             Scaling factor for each coordinate pair.
% alpha         Area distortion ratio, i.e. the ratio between the area of a
%               small rectangle in the grid and its corresponing area on
%               the globe. 
%
% REMARKS:
%  Universal Polar Stereographic projection is a special kind using a polar
%  scaling factor of 0.994 which correspond to a standard parallel of
%  -81.114527778 and using a False easting and northing of 2000000.
%  Standard meridian is the Greenwich meridian. The Universal Polar
%  Stereographic projecttion complements the Universal Transverese Mercator
%  (UTM) projection for polar region and is only applicable between 
%  90 degrees south and 80 degrees south, and 90 degrees north and
%  84 degrees north.
% 
% VERSION:  1.0
% 
% IMH:
% LINE 196: SWITCHED Y SIGN FOR NORTH
% REMOVED TRANSPOSE OF OUTPUTS
%
% AUTHOR:
% Rickard Pettersson, Oct 2004
% Department of Physics, St. Olaf College, MN, USA
% petterss@stolaf.edu
%
% REVISION HISTORY:
%
% COMPABILITY:
% MATLAB 7 (R14)
%
% DEPENDENCIES:
%
% SUBROUTINES:
%
% REFRENCES:
% Snyder J., 1987. Map projections -A working manual. USGS professinal
% paper, 1395, 383p.
%--------------------------------------------------------------------------
minarg = 2;
maxarg = 12;

% Check input
error(nargchk(minarg,maxarg,nargin));

% requisite arguments
lon = varargin{1};
lat = varargin{2};

% Default option values
FN = 0;                 % False northing in meter
FE = 0;                 % False easting in meter
slat = 70;             % Standard parallel in degrees
slon = -45;               % Standard meridian in degrees
% WGS84 Ellipsoid
a = 6378173;            % major axis of ellipsoid
f = 1/298.257223563;    % Flattening of ellipsoid

% Alteration of default values
for n = (minarg+1:2:nargin-1);
        % Checks that every value is a numeric value and every property
        % is a sting
        if ~isstr(varargin{n})
            error('Invalid parameter/value pair arguments.');                 
        % These properties take string arguments otherwise it must be
        % numeric value
        elseif ~isnumeric(varargin{n+1})... 
                & ~strcmpi(varargin{n},'ellipsoid')
            msg = sprintf(['Bad value for property: '...
                           '''' varargin{n} '''\n unknown option.']); 
            error(msg);
        end
       
        % Sets property value
        switch lower(varargin{n})
            case 'falsenorthing'
                FN = varargin{n+1};
            case 'falseeasting'
                FE = varargin{n+1};
            case 'standardparallel'
                slat = varargin{n+1};
            case 'standardmeridian'
                slon = varargin{n+1};
            case 'ellipsoid'
                    switch lower(varargin{n+1})
                        case 'wgs-84'
                            % WGS84 Ellipsoid
                            a = 6378173;
                            f = 1/298.257223563;
                        case 'international'
                            % International ellipsoid
                            a = 6378388;
                            f = 1/297;
                        otherwise
                            msg = sprintf(['Bad value for property: '...
                            '''' string '''\n unknown ellipsoid.']); 
                            error(msg);
                    end
        end
end

% Make sure that both the standard parallel and the data are on the same
% hemisphere
if sign(slat) ~= sign(lat)
    error('Standard parallel and input must be on the same hemisphere');
end

% Convert input to radians
deg2rad = pi/180;
lon = deg2rad .* lon;
lat = deg2rad .* lat;
slon = deg2rad .* slon;
slat = deg2rad .* slat;

% Calculate some constants
b  = a*(1-f);                 % Minor axis of ellipsoid
e  = sqrt(1-(b^2/a^2));       % Eccentricity if ellipsoid
e2 = f*(2-f);                 % Squared eccentricity of ellipsoid

% Determine hemisphere 
if (slat < 0)
  sn   = -1.0;
  lat  = -lat;
  slat = -slat;
else
  sn    = 1.0;
  lat   = lat;
  slat  = slat;
end

% for all values to transform (eqs. 15-9 and 14-15, Snyder, 1987)
t = ((1 - e * sin(lat)) ./ (1 + e * sin(lat))).^(e*0.5);
t = tan((pi/4) - 0.5.*lat) ./t;
m = cos(lat) ./ sqrt(1 - e2 * sin(lat).^2);

% for center latitude (where the scale factor is 1)
tc = ((1 - e*sin(slat)) / (1 + e*sin(slat)))^(e*0.5);
tc = tan((pi/4) - 0.5*slat) /tc;
mc = cos(slat)./sqrt(1 - e2*sin(slat).^2);

% Eq. 21-34, Snyder, 1987
rho = a * mc .* t ./ tc;

% Calculate the coordinates (eqs. 21-30 and 21-31, Snyder, 1987)
x = (sin(lon-slon)).*rho + FE;
y = (cos(lon-slon)).*rho + FN;

% Scale factor at each location (eq.21-32, Snyder, 1987) 
k = rho ./ (a .* m);

% Scale factor at the pole (eq. 21-35). For some reason I need to remove
% 'a' from the denominator to get correct values. 
i = find(lat == pi/2);
k(i) = (1/2)*mc*sqrt((((1 + e)^(1 + e))*((1 - e)^(1 - e))))/tc;

% Area distortion ratio, i.e. the ratio between the area of a small
% rectangle in the grid and its corresponing area on the globe.
alpha = ((1 + sin(lat))./(1 + sin(slat))).^2;


if sn == 1
    y = -y; %SWITCH SIGN FOR NORTH (IMH) 
end


% Sets output
if nargout < 2
    varargout{1} = [x y];
else
    varargout{1} = x;
    varargout{2} = y;
    varargout{3} = k;
    varargout{4} = alpha;
end