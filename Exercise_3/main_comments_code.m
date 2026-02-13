% Positioning and Location Based Services
% A.A. 2024/2025
% 3rd Exercise: Ionospheric Delay Computation
%
% @author: Marianna Alghisi

%% Load Parameters
line1 = 'GPSA 7.4506D-09 1.4901D-08 -5.9605D-08 -1.1921D-07 IONOSPHERIC CORR';
line2 = 'GPSB 9.2160D+04 1.3107D+05 -6.5536D+04 -5.2429D+05 IONOSPHERIC CORR';
ionoparams = [cell2mat(textscan(line1, '%*s %f %f %f %f %*s')) ...
              cell2mat(textscan(line2, '%*s %f %f %f %f %*s'))];

%% Zenithal Ionospheric Error Maps
% Initialize values for the zenithal cycle

% Satellite elevation
elev = 90; % degrees

% Time steps for computation
T1 = [0 6 12 18] * 3600; % seconds

% Latitude and longitude ranges
lat1 = -80:0.5:80; % degrees
long1 = -180:0.5:180; % degrees

% Azimuth (90 degrees is not relevant for zenithal map computation)
azimuth1 = 0;

% Initialize matrix for ionospheric corrections
ionospheric_map1 = zeros(length(long1), length(lat1), length(T1));

% Time, phi, and lambda cycles
for timeIdx = 1:length(T1)
    for phiIdx = 1:length(lat1)
        for lamIdx = 1:length(long1)
            ionospheric_map1(lamIdx, phiIdx, timeIdx) = iono_error_correction(lat1(phiIdx), long1(lamIdx), azimuth1, elev, T1(timeIdx), ionoparams);
        end
    end
end

% Plot zenithal ionospheric error maps
[phi_grid, lambda_grid] = meshgrid(lat1, long1);
figure(1);
title('Ionospheric Error Maps');
for i = 1:length(T1)
    subplot(2, 2, i);
    geoshow(phi_grid, lambda_grid, ionospheric_map1(:, :, i), 'DisplayType', 'texturemap', 'FaceAlpha', 0.5);
    hold on;
    geoshow('landareas.shp', 'FaceColor', 'none');
    title(['Time = ', num2str(T1(i) / 3600), ':00']);
    xlabel('Longitude [deg]');
    ylabel('Latitude [deg]');
    xlim([-180 180]);
    ylim([-80 80]);
    colormap(jet);
end
hp4 = get(subplot(2, 2, 4), 'Position');
colorbar('Position', [hp4(1) + hp4(3) + 0.028, hp4(2), 0.03, hp4(2) + hp4(3) * 2.1]);


%% Polar Ionospheric Error Maps for Milan
% Milan position (latitude and longitude in degrees)
phi2 = 45 + 28 / 60 + 38.28 / 3600; % degrees
lambda2 = 9 + 10 / 60 + 53.40 / 3600; % degrees

% Elevation and azimuth ranges
elevation2 = 0:0.5:90; % degrees
azimuth2 = -180:0.5:180; % degrees

% Time steps for polar maps
time2 = [0 12] * 3600; % seconds

% Initialize matrix for Milan's ionospheric corrections
Iono_map2 = zeros(length(elevation2), length(azimuth2), length(time2));

% Compute ionospheric error corrections for Milan polar maps
for tIdx = 1:length(time2)
    for eIdx = 1:length(elevation2)
        for azIdx = 1:length(azimuth2)
            Iono_map2(eIdx, azIdx, tIdx) = iono_error_correction(phi2, lambda2, azimuth2(azIdx), elevation2(eIdx), time2(tIdx), ionoparams);
        end
    end
end

% Plot polar ionospheric error maps
[Az, El] = meshgrid(azimuth2, elevation2);
for i = 1:length(time2)
    figure(i + 1);
    title(['Ionospheric Error Polar Map for Milan Observer, Time = ', num2str(time2(i) / 3600), ':00']);
    axesm('eqaazim', 'MapLatLimit', [0 90]);
    axis off;
    framem on;
    gridm on;
    mlabel on;
    plabel on;
    setm(gca, 'MLabelParallel', 0);
    geoshow(El, Az, Iono_map2(:, :, i), 'DisplayType', 'texturemap', 'FaceAlpha', 0.6);
    colormap(jet);
    hcb = colorbar('eastoutside');
    set(get(hcb, 'XLabel'), 'String', 'Legend');
end
