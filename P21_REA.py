import pyvista as pv
from vtkmodules.all import *
import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon, MultiPolygon
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import scipy.stats as stats


def transl_matrix(r: float, n_sides: float) -> np.ndarray:

    theta = 2*np.pi / n_sides

    angles = np.arange(0, 2*np.pi, theta)

    n_angles = len(angles)

    matrix = np.zeros((n_angles, 4, 4))

    x = r*np.sin(angles)
    y = r * np.cos(angles)
    z = np.zeros_like(x)

    matrix[:] = np.identity(4)

    matrix[:, 0, -1] = x
    matrix[:, 1, -1] = y
    matrix[:, 2, -1] = z

    return matrix


def rot_matr(centers, angle):
    theta = np.deg2rad(angle)
    n_centers = len(centers)
    matrix = np.zeros((n_centers, 4, 4))
    for i, center in enumerate(centers):
        x = center[0]
        y = center[1]

        T = np.array([
         [1, 0, 0, -x],
         [0, 1, 0, -y],
         [0, 0, 1, 0],
         [0, 0, 0, 1]
        ]).T
        R = np.array([
         [np.cos(theta), -np.sin(theta), 0, 0],
         [np.sin(theta), np.cos(theta), 0, 0],
         [0, 0, 1, 0],
         [0, 0, 0, 1]
        ]).T
        T_1 = np.array([
         [1, 0, 0, x],
         [0, 1, 0, y],
         [0, 0, 1, 0],
         [0, 0, 0, 1]
        ]).T

        matrix[i] = T.dot(R).dot(T_1)

    return matrix


def centers(start_point: np.ndarray, stop_point: np.ndarray, r: float, n_sides: float) -> np.ndarray:

    theta_rad = 2*np.pi / n_sides
    d = 2*r*np.cos(theta_rad/2)
    dx = d*np.sin(theta_rad/2)
    dy = d*np.cos(theta_rad/2)
    x = np.arange(start_point[0], stop_point[0], dx)

    if n_sides == 3:
        y = np.arange(start_point[1], stop_point[1], 3 * dy)
    else:
        y = np.arange(start_point[1], stop_point[1], 2*dy)

    xv, yv = np.meshgrid(x, y)

    yv[:, 1::2] += dy
    if n_sides == 3:
        xv[1::2, :] += dx

    zv = np.zeros_like(xv)
    ones = np.ones_like(xv)
    centers_arr = np.hstack((xv.reshape(-1, 1), yv.reshape(-1, 1), zv.reshape(-1, 1), ones.reshape(-1, 1)))
    return centers_arr

# x_dim = 300 # grid dimension along x
# y_dim = 300 # grid dimension along y
#
# x_res = 20 # the resolution is the n of cells along the direction (so for one cell of 15x15 we have 20 cells)
# y_res = 20
#
# grid = pv.Plane(i_size=x_dim,
#                 j_size=y_dim,
#                 i_resolution=x_res,
#                 j_resolution=y_res)


boundary = gpd.read_file("Interpretation-boundary.shp")
fractures = gpd.read_file("FN_set_1.shp")


# boundary.plot()
xmin, ymin, xmax, ymax = boundary.total_bounds


start = np.array([xmin, ymin, 0], dtype=float)
stop = np.array([xmax, ymax, 0], dtype=float)

n_sides = 6

# print(np.linspace(1, 15, 20))

P21_all_data = pd.DataFrame(columns=['scan_area_size', 'P21'])

for l in np.linspace(1, 26, 26):

    rad = l / (2 * np.sin(np.pi / n_sides))

    h = rad * np.cos(np.pi / n_sides)

    area_target = np.round(n_sides * (l * h) / 2, 6)

    print(f'Calculating on polygon of n={n_sides}, length={l} and area={area_target}')

    centers_array = centers(start, stop, rad, n_sides)

    trans_matrix = transl_matrix(rad, n_sides)

    points = np.einsum('kmn, in->ikm', trans_matrix, centers_array, optimize='optimal')[:, :, :-1]

    polygons = np.empty(len(points), dtype=object)

    for index, p in enumerate(points):
        polygons[index] = Polygon(p)

    grid = gpd.GeoDataFrame({'geometry': polygons}, crs=boundary.crs)

    grid_clip = gpd.clip(grid, boundary)

    grid_clip['area'] = np.round(grid_clip.geometry.area, 6)

    grid_clip = grid_clip[grid_clip['area'] == area_target]
    p21 = np.zeros_like(grid_clip.geometry.area)

    for line, square in enumerate(grid_clip.geometry):

        fractures_clip = gpd.clip(fractures, square)
        total_length = np.sum(fractures_clip.geometry.length.values)
        p21[line] = total_length/area_target

    #grid_clip['p21'] = p21 #activate these three lines if you want to see all the scan area grids
    #grid_clip.plot(column='p21')
    #plt.show()

    # print(grid_clip)
    # fractures_clip.plot()

    x_vals = np.repeat(l, len(p21))

    #plt.plot(x_vals, p21, 'ko')
    #plt.xlabel('Cell side length [m]')
    #plt.ylabel('P21 [m^-1]')

    #Update P21 dataframe

    for i in range(len(p21)):
        # Create a new row
        new_row = pd.DataFrame({'scan_area_size': [x_vals[i]], 'P21': [p21[i]]})

        # Concatenate the new row
        P21_all_data = pd.concat([P21_all_data, new_row], ignore_index=True)

#plt.show()

#get delta IQR values

IQR_results = []
wisker_result = []
unique_scan_area_size = P21_all_data['scan_area_size'].unique()

for value in unique_scan_area_size:
    IQR_first = P21_all_data.loc[P21_all_data['scan_area_size'] == value, 'P21'].quantile(0.75) - P21_all_data.loc[P21_all_data['scan_area_size'] == value, 'P21'].quantile(0.25)
    IQR_second = P21_all_data.loc[P21_all_data['scan_area_size'] == value+1, 'P21'].quantile(0.75) - P21_all_data.loc[P21_all_data['scan_area_size'] == value+1, 'P21'].quantile(0.25)
    delta_IQR = IQR_second - IQR_first
    print("delta_IQR = ", delta_IQR)
    IQR_results.append({'scan_area_size': value + 0.5, 'delta_IQR': delta_IQR}) # sommo 0.5 alla scan area size per posizionare i valori di delta IQR a metà tra i boxplot, questo va bene per il passo di campionamento che utilizziamo adesso
IQR_results = pd.DataFrame(IQR_results)
IQR_results.dropna()
print(IQR_results)

#get delta wiskers

for value in unique_scan_area_size:
    IQR = P21_all_data.loc[P21_all_data['scan_area_size'] == value, 'P21'].quantile(0.75) - P21_all_data.loc[P21_all_data['scan_area_size'] == value, 'P21'].quantile(0.25)
    lower_wisker = P21_all_data.loc[P21_all_data['scan_area_size'] == value, 'P21'].quantile(0.25) - 1.5*IQR
    upper_wisker = P21_all_data.loc[P21_all_data['scan_area_size'] == value, 'P21'].quantile(0.75) + 1.5*IQR
    delta_wisker = upper_wisker - lower_wisker
    wisker_result.append({'scan_area_size': value, 'delta_wiskers': delta_wisker})
wisker_results = pd.DataFrame(wisker_result)
print(wisker_results)

#label per frequenza assouta
x_label = P21_all_data['scan_area_size'].unique()
y_label = np.full(len(P21_all_data['scan_area_size'].unique()), -0.15) #absolute frequency value position
label = [f'{x}' for x in P21_all_data.groupby('scan_area_size').count()['P21'].to_list()]

# create box plots, delta IQR polt and delta wiskers plot
fig = plt.figure(figsize=(8, 12))  # Overall figure size
gs = gridspec.GridSpec(3, 1, height_ratios=[2.5, 1, 1]) #set relative dimension of the three plots

ax1 = fig.add_subplot(gs[0])
P21_all_data.boxplot(column='P21', by='scan_area_size', grid=False, ax=ax1)
for idx in range(len(label)):
    ax1.text(x_label[idx], y_label[idx], label[idx])
ax1.set_title("P21 REA")
ax1.set
ax1.set_xticklabels([])
ax1.set_ylabel("P21 [m⁻¹]")
ax1.set_xlabel("")
ax1.set_xlim(0, 27)
#ax1.grid(True)
ax1.axhline(y=0.715, color='r', linestyle='--', label='Fixed Line (y=10)') #Outcrop mean p21, outcrop area is equivalent to an hexagon with a 65.42m side
plt.suptitle('')

ax2 = fig.add_subplot(gs[1])
IQR_results.plot(x='scan_area_size', y='delta_IQR', kind='line', ax=ax2)
#ax2.set_title("delta IQR")
ax2.set_xticklabels([])
ax2.set_xlabel("")
ax2.set_xlim(0, 27)
ax2.set_ylabel("delta IQR []")
#ax2.grid(True)

ax3 = fig.add_subplot(gs[2])
wisker_results.plot(x='scan_area_size', y='delta_wiskers', kind='line', ax=ax3, color='red')
#ax3.set_title("delta wiskers")
ax3.set_xlabel("scan area edge length [m]")
ax3.set_xlim(0, 27)
ax3.set_ylabel("delta wiskers []")
#ax2.grid(True)

plt.show()
