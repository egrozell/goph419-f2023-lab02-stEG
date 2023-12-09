import numpy as np
import matplotlib.pyplot as plt
from linalg import linalg_interp as l

#import air and water density data
air_data = np.loadtxt('../data/air_density_vs_temp_eng_toolbox.txt')
wtr_data = np.loadtxt('../data/water_density_vs_temp_usgs.txt')
#solving spline interpolation functions
#air density data
x1 = air_data[:,0]
y1 = air_data[:,1]
spline1 = l.spline_function(x1, y1, order = 1)
spline2 = l.spline_function(x1, y1, order = 2)
spline3 = l.spline_function(x1, y1, order = 3)

#water density data
x2 = wtr_data[:,0]
y2 = wtr_data[:,1]
spline4 = l.spline_function(x2, y2, order = 1)
spline5 = l.spline_function(x2, y2, order = 2)
spline6 = l.spline_function(x2, y2, order = 3)

#obtaining 100 equally spaced temperature values between min and max values
air_T = np.linspace(np.min(x1),np.max(x1), 100)
wtr_T = np.linspace(np.min(x2), np.max(x2), 100)

fig, sub_fig =plt.subplots(2,3, figsize = (15,10))

#order = 1
y_air_1 = spline1(air_T)
sub_fig[0][0].plot(air_T,y_air_1,'r')
sub_fig[0][0].plot(x1,y1, 'ro')
sub_fig[0][0].set_xlabel('Temperature (°C)')
sub_fig[0][0].set_ylabel('Water density (g/cm^3)')
sub_fig[0][0].set_title('Water density vs temp (linear spline function)')

#order =2
y_air_2 = spline2(air_T)
sub_fig[0][1].plot(air_T,y_air_2,'b')
sub_fig[0][1].plot(x1,y1, 'bs')
sub_fig[0][1].set_xlabel('Temperature (°C)')
sub_fig[0][1].set_ylabel('Water density (g/cm^3)')
sub_fig[0][1].set_title('Water density vs temp (quadratic spline function)')

#order = 3
y_air_3 = spline3(air_T)
sub_fig[0][2].plot(air_T,y_air_3,'tab:blue')
sub_fig[0][2].plot(x1,y1, 'c^')
sub_fig[0][2].set_xlabel('Temperature (°C)')
sub_fig[0][2].set_ylabel('Water density (g/cm^3)')
sub_fig[0][2].set_title('Water density vs temp (cubic spline function)')

#order 1
y_wtr_1 = spline4(wtr_T)
sub_fig[1][0].plot(wtr_T,y_wtr_1,'r')
sub_fig[1][0].plot(x2,y2, 'ro')
sub_fig[1][0].set_xlabel('Temperature (°C)')
sub_fig[1][0].set_ylabel('Air density (kg/m^3)')
sub_fig[1][0].set_title('Air density vs temp (linear spline function)')

#order =2
y_wtr_2 = spline5(wtr_T)
sub_fig[1][1].plot(wtr_T,y_wtr_2,'b')
sub_fig[1][1].plot(x2,y2, 'bs')
sub_fig[1][1].set_xlabel('Temperature (°C)')
sub_fig[1][1].set_ylabel('Air density (kg/m^3)')
sub_fig[1][1].set_title('Air density vs temp (quadratic spline function)')

#order = 3
y_wtr_3 = spline6(wtr_T)
sub_fig[1][2].plot(wtr_T,y_wtr_3,'tab:blue')
sub_fig[1][2].plot(x2,y2, 'c^')
sub_fig[1][2].set_xlabel('Temperature (°C)')
sub_fig[1][2].set_ylabel('Air density (kg/m^3)')
sub_fig[1][2].set_title('Air density vs temp (cubic spline function)')

plt.savefig('../figures/density vs temp graphs')
