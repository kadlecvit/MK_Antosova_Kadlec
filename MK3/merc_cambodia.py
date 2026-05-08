from math import *
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd

SHP_PATH = r"C:\Users\vanto\Desktop\Skola\UK\Magistr\Matkarto\Úloha_3\SHP\Cambodia.shp"

#Points on the quator
u1 = 12.591177838062*pi/180
v1 = 106.243264727677*pi/180
u2 = 12.424203418023*pi/180
v2 = 103.722519547293*pi/180

#Northernmost point
u3 = 14.447780000349*pi/180
v3 = 103.659150000327*pi/180

#Southernmost point
u4 = 10.543094006121*pi/180
v4 = 106.296125694414*pi/180

#Pole
vk = atan2(tan(u1)*cos(v2)-tan(u2)*cos(v1), tan(u2)*sin(v1)-tan(u1)*sin(v2))
uk = atan(-1/tan(u2)*cos(vk-v2))

#Transformation to the oblique aspect
s1 = asin(sin(u1)*sin(uk)+cos(u1)*cos(uk)*cos(vk-v1))
s2 = asin(sin(u2)*sin(uk)+cos(u2)*cos(uk)*cos(vk-v2))
s3 = asin(sin(u3)*sin(uk)+cos(u3)*cos(uk)*cos(vk-v3))
s4 = asin(sin(u4)*sin(uk)+cos(u4)*cos(uk)*cos(vk-v4))

#True parallel
s0 = acos(2*cos(s3)/(1+cos(s3)))

#Scales
m1=cos(s0)/cos(s1)
m2=cos(s0)/cos(s2)
m3=cos(s0)/cos(s3)
m4=cos(s0)/cos(s4)

#Distortions
ny1=(m1-1)*1000
ny2=(m2-1)*1000
ny3=(m3-1)*1000
ny4=(m4-1)*1000
print(ny1, ny2, ny3, ny4)

#Convert geographic coordinates to map coordinates
def uv_sd(u, v):
    s = np.arcsin(np.sin(u) * np.sin(uk) + np.cos(u) * np.cos(uk) * np.cos(vk - v))
    d = np.arctan2(np.cos(u) * np.sin(v - vk), np.sin(u) * np.cos(uk) - np.cos(u) * np.sin(uk) * np.cos(vk - v))
    return s,d 


#Load country boundary
gdf = gpd.read_file(SHP_PATH)
pts = []
for geom in gdf.geometry:
    if geom.geom_type == 'Polygon':
        pts.extend([[y, x] for x, y in geom.exterior.coords])

u_b = np.array([radians(p[0]) for p in pts])
v_b = np.array([radians(p[1]) for p in pts])
s_b, d_b = uv_sd(u_b, v_b)
d_center = np.mean(d_b)

#Country boundary coordinates
bx = cos(s0) * (d_b - d_center)
by = -np.log(np.tan(s_b/2 + pi/4))

#Draw map
fig, ax = plt.subplots(figsize=(6, 11))

#Graticule 
buf = radians(2)
for v_g in np.arange(v_b.min()-buf, v_b.max()+2*buf, radians(2)):
    u_g = np.linspace(u_b.min()-buf, u_b.max()+buf, 200)
    sg, dg = uv_sd(u_g, np.full_like(u_g, v_g))
    gx = cos(s0) * (dg - d_center)
    gy = -np.log(np.tan(sg/2 + pi/4))
    ax.plot(gx, gy, 'k-', lw=0.3)

for u_g in np.arange(u_b.min()-buf, u_b.max()+2*buf, radians(2)):
    v_g = np.linspace(v_b.min()-buf, v_b.max()+buf, 200)
    sg, dg = uv_sd(np.full_like(v_g, u_g), v_g)
    gx = cos(s0) * (dg - d_center)
    gy = -np.log(np.tan(sg/2 + pi/4))
    ax.plot(gx, gy, 'k-', lw=0.3)

#Equideformates
X, Y = np.meshgrid(np.linspace(bx.min()-0.05, bx.max()+0.05, 400),
                   np.linspace(by.min()-0.05, by.max()+0.05, 400))
S_g = 2*np.arctan(np.exp(-Y)) - pi/2
NU  = (cos(s0)/np.cos(S_g) - 1) * 1000

ny_bnd = max(abs(ny1), abs(ny2), abs(ny3), abs(ny4))
step   = round(ny_bnd/4, 1)
levels = np.round(np.arange((cos(s0)-1)*1000 - step, ny_bnd + step*1.5, step), 1)

cs = ax.contour(X, Y, NU, levels=levels, colors='red', linewidths=0.8)
ax.clabel(cs, fmt='%.1f', fontsize=7)

# Country boundary
ax.plot(bx, by, 'b-', lw=1.5)

ax.set_aspect('equal')
plt.show()