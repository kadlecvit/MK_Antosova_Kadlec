import pyproj

crs_geo = pyproj.CRS.from_epsg(4326)

# Gnomonic projection with north pole at Lat = 90°, Lon = 0°
crs_gnom = pyproj.CRS.from_proj4('+proj=gnom +lat_0=90 +lon_0=0 +datum=WGS84')

p = pyproj.Proj(crs_gnom)

gnom_f = p.get_factors(36, 52.6226, False, True)

print(f"Meridional scale: {gnom_f.meridional_scale}")
print(f"Parallel scale: {gnom_f.parallel_scale}")
print(f"Tissot's indicatrix - semi-major axis a: {gnom_f.tissot_semimajor}")
print(f"Tissot's indicatrix - semi-minor axis b: {gnom_f.tissot_semiminor}")
print(f"Angle between meridian and parallel: {gnom_f.meridian_parallel_angle}°")
print(f"Angular distortion: {gnom_f.angular_distortion}°")
print(f"Areal scale: {gnom_f.areal_scale}")
print(f"Meridian convergence: {gnom_f.meridian_convergence}°")
