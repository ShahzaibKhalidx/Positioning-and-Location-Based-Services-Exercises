# -*- coding: utf-8 -*-
"""
Positioning and Location Based Services
A.A. 2023/2024
3rd Exercise: Ionospheric delay computation

@author: Marianna Alghisi
"""

'''
Goals: 
   1) 4 zenithal maps of ionospheric error corrections:
   - Elevation: 90°
   - Latitude: [-80°, 80°, step = 0.5°]
   - Longitude: [-180°, 180°, step = 0.5°]
   - time: [00:00, 06:00, 12:00, 18:00]
   
   2) 2 polar maps of ionospheric error corrections:
    - Observer located in Milan
    - Elevation: [0, 90°, step = 0.5°]
    - Azimuth: [-180°, 180°, step = 0.5°]
    time: [00:00, 12:00]
'''

# Import required libraries
import numpy as np
import matplotlib.pyplot as plt
import ionoCorrection as ic

# Ionospheric correction parameters:
alpha = [7.4506*10**(-9), 1.4901*10**(-8), -5.9605*10**(-8), -1.1921*10**(-7)]
beta = [9.2160*10**(4), 1.3107*10**(5), -6.5536*10**(4), -5.2429*10**(5)]

'''
1) ZENITHAL MAPS
'''
# Initialization of the parameters: define inputs for the Zenithal maps (elevation, azimuth, time, lat, lon)

# Loop on time, latitude and longitude --> compute for each point the Ionospheric delay
# TIP: store latitude, longitude and iono_delay in list objects!

# SUGGESTION FOR PLOTTING

timePrint = ['00:00', '06:00', '12:00', '18:00']

for i in range(4):
    IONO = #list of iono delays for considered epoch
    results = pd.DataFrame()
    results['latitude'] = #list of latitude
    results['longitude'] = #list of longitude
    results['iono_delay'] = IONO
    gdf = gpd.GeoDataFrame(results, geometry=gpd.points_from_xy(results.longitude, results.latitude), crs = 3857)
    world = gpd.read_file('world/world.shp')
    fig, ax = plt.subplots (figsize = (15,15))
    world.boundary.plot(ax=ax, color='black')
    ax.set(xlabel='Longitude', ylabel='Latitude', title='Zenithal map of ionospheric delay at '+ str(timePrint[i]))
    gdf.plot(column='iono_delay', ax = ax, marker='o', legend=True)

'''
2) POLAR MAPS
'''
# Definition of Milano's coordinates
lat_mi = 45 + 28/60 + 38.28/60**2
lon_mi = 9 + 10/60 + 53.4/60**2

# Inizialization of the parameters for the loop: time, elevation, azimuth

# Loop on the parameters ---> compute for each point the Ionospheric delay

# SUGGESTION FOR PLOTTING

timePrint = ['00:00', '12:00']

for i in [0,1]:
    t = timePrint[i]
    fig = plt.figure(figsize = (15,15))
    ax = fig.add_subplot(projection='polar')
    plt.scatter(azimuth, elevation, c=all_iono_delays[i], cmap='brg',  alpha=0.75, label=all_iono_delays[i])
    ax.set_title('Ionospheric Error Polar Map for Milan Observer time = '+str(t))
    plt.colorbar(label='Ionospheric Delay')
    plt.show()