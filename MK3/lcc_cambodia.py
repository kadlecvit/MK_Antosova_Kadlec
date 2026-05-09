from math import *
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd

SHP_PATH = r"C:\Users\kadle\OneDrive - CUNI\A-Dokumenty\.Uni\_4.MGR\Matkarto\3_uloha\SHP\Cambodia.shp"

R = 1

#Pole
uk = 16.4267 * pi/180
vk = 104.3286 * pi/180

#Northernmost point
u1 = 14.458760 * pi/180
v1 = 103.659814 * pi/180

#Southernmost point
u2 = 10.422705 * pi/180
v2 = 104.470327 * pi/180

#Transformation to the oblique aspect
s1 = asin(sin(u1) * sin(uk) + cos(u1) * cos(uk) * cos(vk - v1))
s2 = asin(sin(u2) * sin(uk) + cos(u2) * cos(uk) * cos(vk - v2))

# Constant c
cn = log(cos(s1)) - log(cos(s2))
cd = log(tan(s2/2 + pi/4)) - log(tan(s1/2 + pi/4))
c = cn / cd

# s0, rho0, rho1, rho2
s0 = asin(c)
rho0_n = 2 * R * cos(s0) * cos(s1) * (tan(s1/2 + pi/4))**c
rho0_d = c * (cos(s0) * (tan(s0/2 + pi/4))**c + cos(s1) * (tan(s1/2 + pi/4))**c)
rho0 = rho0_n / rho0_d

rho1 = rho0 * ((tan(s0/2 + pi/4)) / (tan(s1/2 + pi/4)))**c
rho2 = rho0 * ((tan(s0/2 + pi/4)) / (tan(s2/2 + pi/4)))**c

#Scales
m1 = (c * rho1) / (R * cos(s1))
m2 = (c * rho2) / (R * cos(s2))
m0 = (c * rho0) / (R * cos(s0))

#Distortions
ny1 = (m1 - 1) * 1000
ny2 = (m2 - 1) * 1000
ny0 = (m0 - 1) * 1000

print (ny1, ny2, ny0)


#Convert geographic coordinates to map coordinates
def uv_to_sd(u, v):
    s = np.arcsin(np.sin(u) * np.sin(uk) + np.cos(u) * np.cos(uk) * np.cos(vk - v))
    d = np.arctan2(
        np.cos(u) * np.sin(v - vk),
        np.sin(u) * np.cos(uk) - np.cos(u) * np.sin(uk) * np.cos(vk - v)
    )
    return s, d


#Load country boundary
gdf = gpd.read_file(SHP_PATH)
pts = []
for geom in gdf.geometry:
    if geom.geom_type == 'Polygon':
        pts.extend([[y, x] for x, y in geom.exterior.coords])
    elif geom.geom_type == 'MultiPolygon':
        for poly in geom.geoms:
            pts.extend([[y, x] for x, y in poly.exterior.coords])

u_b = np.array([radians(p[0]) for p in pts])
v_b = np.array([radians(p[1]) for p in pts])

s_b, d_b = uv_to_sd(u_b, v_b)
d_center = np.mean(d_b)

#Country boundary coordinates
def sd_to_xy(s, d):
    rho = rho0 * (np.tan(s0/2 + pi/4) / np.tan(s/2 + pi/4))**c
    theta = c * (d - d_center)
    x = rho * np.sin(theta)
    y = rho0 - rho * np.cos(theta)
    return x, y


bx, by = sd_to_xy(s_b, d_b)

# Mirroring X-axis
bx = -bx


#Draw map
fig, ax = plt.subplots(figsize=(8, 8))

#Graticule 
buf = radians(2)
u_min, u_max = u_b.min() - buf, u_b.max() + buf
v_min, v_max = v_b.min() - buf, v_b.max() + buf

for v_g in np.arange(floor(degrees(v_min) / 2) * 2, ceil(degrees(v_max) / 2) * 2 + 2, 2):
    u_g = np.linspace(u_min, u_max, 300)
    v_g_arr = np.full_like(u_g, radians(v_g))
    s_g, d_g = uv_to_sd(u_g, v_g_arr)
    gx, gy = sd_to_xy(s_g, d_g)
    ax.plot(-gx, gy, 'k-', lw=0.3)

for u_g in np.arange(floor(degrees(u_min) / 2) * 2, ceil(degrees(u_max) / 2) * 2 + 2, 2):
    v_g = np.linspace(v_min, v_max, 300)
    u_g_arr = np.full_like(v_g, radians(u_g))
    s_g, d_g = uv_to_sd(u_g_arr, v_g)
    gx, gy = sd_to_xy(s_g, d_g)
    ax.plot(-gx, gy, 'k-', lw=0.3)


#Equideformates
X, Y = np.meshgrid(
    np.linspace(bx.min() - 0.05, bx.max() + 0.05, 500),
    np.linspace(by.min() - 0.05, by.max() + 0.05, 500))

# Backward transformation
RHO = np.sqrt(X**2 + (rho0 - Y)**2)
tan_half = np.tan(s0/2 + pi/4) * (rho0 / RHO)**(1.0 / c)
S_g = 2 * np.arctan(tan_half) - pi/2

M = (c * RHO) / (R * np.cos(S_g))
NU = (M - 1) * 1000

ny_bnd = max(abs(ny1), abs(ny2), abs(ny0))
step   = round(ny_bnd/4, 1)
levels = np.round(np.arange((m0 - 1)*1000 - step*2, ny_bnd + step*2, step), 2)
if 0.0 not in levels:
    levels = np.sort(np.append(levels, 0.0))

cs = ax.contour(X, Y, NU, levels=levels, colors='red', linewidths=0.8)
ax.clabel(cs, fmt='%.2f', fontsize=7)


# Country boundary
ax.plot(bx, by, 'b-', lw=1.5)

ax.set_aspect('equal')
plt.show()