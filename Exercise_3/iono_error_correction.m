% Positioning and Location Based Services
% A.A. 2024/2025
% 3rd Exercise: Ionospheric delay computation
%
% @author: Marianna Alghisi

function [corr] = iono_error_correction(lat, lon, az, el, time_rx, ionoparams, sbas)
 
% SYNTAX:
%   [corr] = iono_error_correction(lat, lon, az, el, time_rx, ionoparams, sbas);
%
% INPUT:
%   lat = receiver latitude    [degrees]
%   lon = receiver longitude   [degrees]
%   az  = satellite azimuth    [degrees]
%   el  = satellite elevation  [degrees]
%   time_rx    = receiver reception time
%   ionoparams = ionospheric correction parameters
%   sbas = SBAS corrections    < optional if available >
%
% OUTPUT:
%   corr = ionospheric error correction [m]
%
% DESCRIPTION:
%   Computation of the pseudorange correction due to ionospheric delay.
%   Klobuchar model or SBAS ionosphere interpolation.

iono_model = 1;
s_light = 299792458;

% If in ATM_model section in config file iono is set to 0 it does not use ionospheric model
if (iono_model == 0)
    corr = zeros(size(el));
else
    % Initialization
    corr = zeros(size(el));
    
    % If ionosphere parameters are available and SBAS corrections are disabled/not available
    if ((nargin == 6) && (sum(abs(ionoparams)) > 0)) || ...
            ((nargin > 6) && (sum(abs(ionoparams)) > 0) && isempty(sbas))
        % Apply Klobuchar ionosphere model
        corr = klobuchar_model(lat, lon, az, el, time_rx, ionoparams);
    elseif ((nargin > 6) && ~isempty(sbas))
        % Apply SBAS interpolated ionospheric delay (where possible)
        corr = sbas_iono_interp(lat, lon, az, el, sbas);
        
        % Detect if some satellites could not be corrected by SBAS
        not_corr = isnan(corr);
        if (any(not_corr))
            % Apply Klobuchar model where SBAS corrections are not available
            corr(not_corr) = klobuchar_model(lat, lon, az(not_corr), el(not_corr), time_rx, ionoparams);
        end
    end
end

% -------------------------------------------------------------------------
% Nested function: Klobuchar model
% -------------------------------------------------------------------------
function [delay] = klobuchar_model(lat, lon, az, el, time_rx, ionoparams)
    % Initialization
    delay = zeros(size(el));

    % Ionospheric parameters
    a0 = ionoparams(1);
    a1 = ionoparams(2);
    a2 = ionoparams(3);
    a3 = ionoparams(4);
    b0 = ionoparams(5);
    b1 = ionoparams(6);
    b2 = ionoparams(7);
    b3 = ionoparams(8);

    % Elevation from 0 to 90 degrees
    el = abs(el);

    % Conversion to semicircles
    lat = lat / 180;
    lon = lon / 180;
    az = az / 180;
    el = el / 180;

    f = 1 + 16 * (0.53 - el).^3;
    psi = (0.0137 ./ (el + 0.11)) - 0.022;

    phi = lat + psi .* cos(az * pi);
    phi(phi > 0.416)  =  0.416;
    phi(phi < -0.416) = -0.416;

    lambda = lon + ((psi .* sin(az * pi)) ./ cos(phi * pi));
    ro = phi + 0.064 * cos((lambda - 1.617) * pi);

    t = lambda * 43200 + time_rx;
    t = mod(t, 86400);

    a = a0 + a1 * ro + a2 * ro.^2 + a3 * ro.^3;
    a(a < 0) = 0;

    p = b0 + b1 * ro + b2 * ro.^2 + b3 * ro.^3;
    p(p < 72000) = 72000;

    x = (2 * pi * (t - 50400)) ./ p;

    % Ionospheric delay
    index = abs(x) < 1.57;
    delay(index) = s_light * f(index) .* (5e-9 + a(index) .* (1 - (x(index).^2) / 2 + (x(index).^4) / 24));

    index = abs(x) >= 1.57;
    delay(index) = s_light * f(index) * 5e-9;
end

% -------------------------------------------------------------------------
% Nested function: SBAS ionospheric interpolation
% -------------------------------------------------------------------------
function [delay] = sbas_iono_interp(lat, lon, az, el, sbas)
    % Initialization
    delay = NaN(size(el));

    % Convert to radians
    lat = lat * pi / 180;
    lon = lon * pi / 180;
    az  = az * pi / 180;
    el  = el * pi / 180;

    for i = 1:length(az)
        % Ionosphere pierce point coordinates and slant factor
        [latpp, lonpp, fpp] = iono_pierce_point(lat, lon, az(i), el(i));

        % Find the nodes of the cell that contains the ionosphere pierce point
        [igp4, tv] = sel_igp(latpp, lonpp, sbas.igp, sbas.lat_igp, sbas.lon_igp);

        % Interpolate the ionospheric vertical delay if valid grid cell exists
        if (sum(isfinite(igp4)) >= 3)
            [ivd_pp] = interp_ivd(igp4, sbas.igp, sbas.ivd, latpp, lonpp, tv);
            delay(i) = fpp * ivd_pp; % m
        end
    end
end

end
