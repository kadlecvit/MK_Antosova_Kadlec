from numpy import *
from pyproj import *
from matplotlib.pyplot import *

def project(proj_name, R_z, lat, lon, lat0, lon0):
    # Create new projection given by proj_name
    my_proj =  Proj(proj=proj_name, R=R_z, lat_1 = lat0, lon_0 = lon0)

    # Project point calculation
    [X,Y] = my_proj(lon, lat)

    #Distortions
    dist = my_proj.get_factors(lon, lat)
    a = dist.tissot_semimajor
    b = dist.tissot_semiminor

    return X, Y, a, b 

def graticule(lat_min, lon_min, lat_max, lon_max, Dlat, Dlon, dlat, dlon, R, lat0, lon0, proj_name):
    #Create graticule of the given map projection
    #Create meridians
    lat_mer = arange(lat_min, lat_max + dlat/2, dlat)
    lon_mer = arange(lon_min, lon_max + Dlon/2, Dlon)
    
    #Create parallels
    lat_par = arange(lat_min, lat_max + Dlat/2, Dlat)
    lon_par = arange(lon_min, lon_max + dlon/2, dlon)

    #Project meridians (each meridian separately, NaN-separated to avoid joining lines)
    Xm_list, Ym_list = [], []
    for lon_m in lon_mer:
        lats = lat_mer
        lons = full_like(lats, lon_m)
        Xs, Ys, _, _ = project(proj_name, R, lats, lons, lat0, lon0)
        Xs, Ys = split_jumps(Xs, Ys)
        Xm_list.append(Xs); Xm_list.append(array([nan]))
        Ym_list.append(Ys); Ym_list.append(array([nan]))
    Xm = concatenate(Xm_list)
    Ym = concatenate(Ym_list)

    #Project parallels (each parallel separately, NaN-separated to avoid joining lines)
    Xp_list, Yp_list = [], []
    for lat_p in lat_par:
        lons = lon_par
        lats = full_like(lons, lat_p)
        Xs, Ys, _, _ = project(proj_name, R, lats, lons, lat0, lon0)
        Xs, Ys = split_jumps(Xs, Ys)
        Xp_list.append(Xs); Xp_list.append(array([nan]))
        Yp_list.append(Ys); Yp_list.append(array([nan]))
    Xp = concatenate(Xp_list)
    Yp = concatenate(Yp_list)
    
    return Xm, Ym, Xp, Yp

def split_jumps(X, Y, factor = 5):
    #Insert NaN where consecutive points jump (antimeridian / projection discontinuity)
    d = sqrt(diff(X)**2 + diff(Y)**2)
    finite = d[isfinite(d) & (d > 0)]
    if len(finite) == 0:
        return X, Y
    threshold = factor * median(finite)
    idx = where(d > threshold)[0] + 1
    X_out, Y_out = list(X), list(Y)
    for i in reversed(idx):
        X_out.insert(i, nan)
        Y_out.insert(i, nan)
    return array(X_out), array(Y_out)

#Define projection
#proj_name = "sinu"
#proj_name = "bonne"
#proj_name = "eck5"
#proj_name = "wintri"
proj_name = "aitoff"

R = 6380000
lat0 = 10
lon0 = -60

#Define projection grid
lat_min = -80
lat_max = 80
lon_min = -180
lon_max = 180
Dlat = 10
Dlon = 10
dlat = 0.1 * Dlat
dlon = 0.1 * Dlon
nlat = 100
nlon = 100

#Create intervals
lat = linspace(lat_min, lat_max, nlat)
lon = linspace(lon_min, lon_max, nlon)

#Create  meshgrid
latg, long = meshgrid(lat, lon)

#Project meshgrid
X, Y, a, b = project(proj_name, R, latg, long, lat0, lon0)

#Airy local
h2_a = 0.5*((a-1)**2+(b-1)**2)

#Complex local
h2_c = 0.5*(abs(a-1)+abs(b-1)) + a/b - 1

#Airy global
H2_a = mean(h2_a)

#Complex global
H2_c = mean(h2_c)

#Airy weighted global
w = cos(latg * pi /180)
H2_aw = sum(w*h2_a)/sum(w)

#Complex weighted global
H2_cw = sum(w*h2_c)/sum(w) 

print(H2_a, H2_c, H2_aw, H2_cw)

#Draw continents
continents = loadtxt(r"C:\Users\kadle\OneDrive - CUNI\A-Dokumenty\.Uni\_4.MGR\Matkarto\4_uloha\Kontinenty\amer_s.txt")

#Extract coordinates
latc = continents[:, 0]
lonc = continents[:, 1]

#Project points
Xc, Yc, ac, bc = project(proj_name, R, latc, lonc, lat0, lon0)

#Draw points
plot(Xc, Yc, linewidth = 2)

#Create meridians and parallels
Xm, Ym, Xp, Yp = graticule(lat_min, lon_min, lat_max, lon_max, Dlat, Dlon, dlat, dlon, R, lat0, lon0, proj_name)

#PLot meridians and parallels
plot(Xm, Ym, color = 'black', linewidth = 0.5)
plot(Xp, Yp, color = 'black', linewidth = 0.5)

#Variable map scale
S = 100000000
Sv = S/a

#Mask cells across projection discontinuities (antimeridian)
d0 = sqrt(diff(X, axis=0)**2 + diff(Y, axis=0)**2)
d1 = sqrt(diff(X, axis=1)**2 + diff(Y, axis=1)**2)
threshold = 5 * median(concatenate([d0.ravel(), d1.ravel()]))
mask = zeros_like(Sv, dtype=bool)
mask[:-1, :] |= (d0 > threshold)
mask[1:, :]  |= (d0 > threshold)
mask[:, :-1] |= (d1 > threshold)
mask[:, 1:]  |= (d1 > threshold)
Sv = ma.masked_array(Sv, mask=mask)

#Create contour lines
dS = arange(20000000, 200000000, 10000000)
contours = contour(X, Y, Sv, levels = dS, colors = 'red')

#Create contour labels
clabel(contours, inline_spacing = -20)

show()