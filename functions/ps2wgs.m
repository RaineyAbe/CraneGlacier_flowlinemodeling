function varargout = ps2wgs(varargin)
%--------------------------------------------------------------------------
% PS2WGS
% -------------------------------------------------------------------------
% SYNTAX:
% [lon,lat] = PS2WGS(x,y)       
% [lon,lat] = PS2WGS(...,'PropertyName',propertyvalue)
%
% DESCRIPTION:
% Converts from polar sterographic coordinates to geographic coordinates
% (latitude, longitude).  
%
% INPUT:  
% x,y                Vectors containing polar stereographic coordinates
%                    given in meters.  
% 'StandardParallel' The standard Parallel for the polar stereographic
%                    projection. Default is 70 degree north.
% 'StandardMeridian' The Standard Meridian for the polar stereographic
%                    projection. This is the meridian for the positive y
%                    direction. Default is the -45 meridian.
% 'Ellipsoid'        Ellipsoid model for the stereographic projection. 
%                    Valid options are:
%                        'WGS-84'
%                        'International'
%                    Default is WGS-84.   
% 'FalseNorthing'    Offset between origin of coordinate system
%                    and origin of projection. Often used to avoid negative
%                    coordinate values. Default is 0.
% 'FalseEasting'     Offset between origin of coordinate system
%                    and origin of projection. Often used to avoid negative
%                    coordinate values. Default is 0.
% 'Threshold'        Error threshold limit for itteration given in radians.
%                    Default is 1e-15 radians.
%
% OUTPUT:
% lon,lat            Vectors containing transformed geographic coordinates 
%                    (longitude and latitude).
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
% AUTHOR:
% Rickard Pettersson, Nov 2004
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
x = varargin{1};
y = varargin{2};

% Default option values
FN = 0;                 % False northing in meter
FE = 0;                 % False easting in meter
slat = 70;             % Standard parallel in degrees
slon = -45;               % Standard meridian in degrees
threshold = 1e-15;      % Threshold for the itteration (in radians)
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
%                             a = 6378137; f = 1/298.257223563; %EPSG:3031 (Antarctic PS)
                        case 'threshold'
                            threshold = varargin{n+1};
                        otherwise
                            msg = sprintf(['Bad value for property: '...
                            '''' string '''\n unknown ellipsoid.']); 
                            error(msg);
                    end
        end
end

% Conversion factors
rad2deg = 180/pi;
deg2rad = pi/180;                            

% Convert input to radians
slon = slon*deg2rad;
slat = slat*deg2rad;

% Calculate some constants for the ellipsoid
b  = a*(1-f);                 % Minor axis of ellipsoid
e  = sqrt(1-(b^2/a^2));       % Eccentricity if ellipsoid
e2 = f*(2-f);                 % Squared eccentricity of ellipsoid

% Remove false easting and northing
x = FE - x;
y = FN - y;

% Eq. 20-18, Snyder, 1987
rho = sqrt(x.^2 + y.^2);

%	First use a approximate inverse calculation. Then we use a itteration
%	of the forward calculation (which is more exact) to correct this
%	approximation. 
    
% Decide on which hemsiphere the Standard latitude is located and turn
% the standard latitude to positive value for the calculations, i.e. do
% the calculations on the northern hemisphere
if (slat < 0)
    sn = -1;
    slat = -slat;
    slon = -slon;
    x = -x;
    y = -y;
else
    sn = 1;
end

y=-y;
% Eq 14-15, 15-9 and 21-40, Snyder, 1987
mc = cos(slat)/sqrt(1 - e2*sin(slat)^2);
tc = tan(pi/4 - 0.5*slat)/((1 - e*sin(slat))/(1 + e*sin(slat)))^(e/2);
    
t = rho.*tc./(a*mc);
    
% Eq 7-13, Snyder 1987
chi = pi/2 - 2*atan(t);
    
% Series expansion for the first approximation, eq. 3-5, Snyder 1987
lat = chi + (0.5*e^2 + 5*e^4/24 + e^6/12 + e^8/360)*sin(2*chi) + ...
            (7*e^4/48 + 29*e^6/240 + 811*e^8/11520)*sin(4*chi) + ...
            (7*e^6/120 + 81*e^8/1120)*sin(6*chi) + ...
            (4279*e^8/161280)*sin(8*chi);
    
% Calculate the first approximation of the longitude
lon = atan2(sn.*x,-sn.*y) + slon;

% Convert longitude to the range -180 to 180
i = find(lon < -pi);
lon(i) = lon(i) + 2*pi;
    
i = find(lon > pi);
lon(i) = lon(i) - 2*pi;

% Store the apporixmated latitude
approxlat = lat;

% Using the approximate result as a starting point, iterate to improve the
% accuracy of the inverse solution
tmp = lat;
while abs(max(tmp) - max(lat)) > threshold
    tmp = lat;
    tmp2 = tmp;
    lat = pi/2 - 2*atan(t*((1 - e*sin(lat))/(1 + e*sin(lat)))^(e/2));

    % Check for convergence, so no infinte loop is created 
    i = find((tmp - lat) > (tmp2 - lat));
    if ~isempty(i)
        warning('Calculations did not converge: Using the approximated value');
        lat(i) = approxlat(i);
    end
 
end

% Convert back to proper hemsiphere and convert do decimal degrees.
lon = sn.*lon*rad2deg;
lat = sn.*lat*rad2deg;
% 
%IMH: CORRECT FOR SEMETRIC ERROR
if sn == 1
    lon = lon-(2.*(lon + 45));
    
else
    if lon > 0
        lon = 180-lon;
    else
        lon = -180-lon;
    end
end

% Make sure the output is a row vector
lon = reshape(lon,length(lon),1);
lat = reshape(lat,length(lat),1);

% Sets output
if nargout < 2
    varargout{1} = [lon lat];
else
    varargout{1} = lon;
    varargout{2} = lat;
end
