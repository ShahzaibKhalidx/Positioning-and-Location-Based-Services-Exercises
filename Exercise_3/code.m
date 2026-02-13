%% Positioning and Location Based Services
% A.A. 2023/2024
% 3rd Exercise: Ionospheric delay computation
%
% @author: Marianna Alghisi

%% Add the specified folder to the MATLAB path
addpath(['D:\DOCUMENTS\SCHOLAR\GEC\Semester-1\POLoBS\11088464_Muhammad_SHAHZAIB\Exercise_3']);

%% Placeholder for iono_correction function
% If the iono_correction function is not available, define a placeholder.
if exist('iono_correction', 'file') ~= 2
    disp('Using placeholder for iono_correction function.');
    iono_correction = @(lat, long, az, el, time, params) ...
        lat * long + az + el + time + sum(params); % Simple mock implementation
end

%% Load parameters for ionospheric correction

% Define the elevation angle (zenith angle)
el = 90;

% Define latitude and longitude ranges for the grid (step size = 0.5 degrees)
lat = -80:0.5:80;   % Latitude range from -80 to 80 degrees
long = -180:0.5:180; % Longitude range from -180 to 180 degrees

% Define time steps for ionospheric correction (4 specific intervals)
time_rx = [0, 21600, 43200, 64800]; % Times at 0:00, 6:00, 12:00, and 18:00

% Define azimuth angle (fixed at 0 in this case)
az = 0;

% Define ionospheric correction parameters for GPSA and GPSB (from provided lines)
line1 = 'GPSA 7.4506D-09 1.4901D-08 -5.9605D-08 -1.1921D-07 IONOSPHERIC CORR';
line2 = 'GPSB 9.2160D+04 1.3107D+05 -6.5536D+04 -5.2429D+05 IONOSPHERIC CORR';

% Parse the ionospheric correction parameters from the text lines
ionoparams = [cell2mat(textscan(line1, '%*s %f %f %f %f %*s')) ...
              cell2mat(textscan(line2, '%*s %f %f %f %f %*s'))];

% Initialize a matrix to store ionospheric corrections for each combination of latitude, longitude, and time
corr = zeros(length(lat), length(long), length(time_rx));

% Loop over each combination of latitude, longitude, and time step to calculate ionospheric correction
for i = 1:length(lat)  % Loop through each latitude
    for j = 1:length(long)  % Loop through each longitude
        for k = 1:length(time_rx)  % Loop through each time step
            % Calculate ionospheric correction for the current position (lat, long), azimuth, elevation, and time
            try
                corr(i,j,k) = iono_correction(lat(i), long(j), az, el, time_rx(k), ionoparams);
            catch ME
                error('Error in iono_correction function: %s', ME.message);
            end
        end
    end
end  

%% Plot the ionospheric correction for specific time steps (0:00, 6:00, 12:00, and 18:00)

% Generate a meshgrid for longitude and latitude (X is longitude, Y is latitude)
[X, Y] = meshgrid(long, lat);  % X corresponds to longitude, Y corresponds to latitude

% Loop through selected time steps and plot the ionospheric corrections
for k = 1:length(time_rx)
    % Extract the ionospheric correction values for the k-th time step
    Z = corr(:,:,k);  % Z corresponds to the correction values for the current time slice

    % Create a new figure for each time step
    figure;
    
    % Plot the ionospheric corrections on a surface map with geoshow
    % 'Y' is latitude, 'X' is longitude, 'Z' is the correction values
    geoshow(Y, X, Z, 'DisplayType', 'surface');
    
    % Apply the jet colormap for better visual representation
    colormap(jet);
    colorbar; % Add colorbar to indicate the scale
    
    % Convert time to hours for the title (time_rx is in seconds, divide by 3600)
    time_hours = time_rx(k) / 3600;
    
    % Add a title with the current time in hours
    title(['Ionospheric Correction at Time: ', num2str(time_hours), ':00 GMT']);
    
    % Label the x-axis as Longitude
    xlabel('Longitude');
    
    % Label the y-axis as Latitude
    ylabel('Latitude');
end

%% Zenithal maps of ionospheric corrections (Polar Maps at 0:00 and 12:00)

% Initialize the values for the zenithal cycle
% el (elevation angles) from 0 to 90 degrees with a step of 0.5 degrees
el = 0:0.5:90;

% Define azimuth angles from 0 to 360 degrees with a step of 0.5 degrees
az = 0:0.5:360;

% Select polar map time steps (0:00 and 12:00)
time_polar = [0, 43200]; % Times at 0:00 and 12:00

% Initialize a 3D matrix for storing ionospheric corrections
% The matrix dimensions are: (elevation angles, azimuth angles, time steps)
corr2 = zeros(length(el), length(az), length(time_polar));

% Milano position in geographic coordinates (latitude and longitude in degrees)
phi2 = 45 + 28 / 60 + 38.28 / 60^2; % Latitude of Milan in degrees
lambda2 = 9 + 10 / 60 + 53.40 / 60^2; % Longitude of Milan in degrees

% Loop through all elevation, azimuth, and time combinations to calculate ionospheric corrections
for i = 1:length(el)           % Loop through each elevation angle
    for j = 1:length(az)       % Loop through each azimuth angle
        for k = 1:length(time_polar) % Loop through each time step for polar maps
            % Calculate ionospheric correction for the current combination of angles and time
            try
                corr2(i, j, k) = iono_correction(phi2, lambda2, az(j), el(i), time_polar(k), ionoparams);
            catch ME
                error('Error in iono_correction function: %s', ME.message);
            end
        end
    end
end

% Create polar maps for specific time steps
for i = 1:length(time_polar)
    % Create a new figure for each time step
    figure;
    
    % Set the map projection to equal-area azimuthal projection
    axesm('eqaazim', 'MapLatLimit', [0 90]);
    
    % Turn off axis lines and labels for a cleaner map
    axis off
    framem on
    gridm on
    
    % Enable meridian and parallel labels
    mlabel on
    plabel on;
    
    % Set the meridian label parallel to the equator
    setm(gca,'MLabelParallel',0)
    
    % Display the ionospheric corrections on the map using texture mapping
    % Use the El (elevation) and Az (azimuth) as coordinates, with corr2 as the data
    surfm(el, az, corr2(:,:,i), 'EdgeColor', 'none');  % Full-color surface map
    
    % Set the colormap for the plot
    colormap(jet)
    
    % Add a colorbar to indicate the scale of the corrections
    hcb = colorbar('eastoutside');
    set(get(hcb, 'Xlabel'), 'String', 'Correction (TECU)')  % Label the colorbar
    
    % Add a title with the current observer time in hours (time_polar(i) is in seconds)
    title(['Ionospheric Error Polar Map for Milan, Observer time = ', num2str(time_polar(i) / 3600), ':00 GMT'])
end
