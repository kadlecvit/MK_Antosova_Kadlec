from math import *
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

SHP_PATH = r"C:\Users\kadle\OneDrive - CUNI\A-Dokumenty\.Uni\_4.MGR\Matkarto\3_uloha\SHP\Vietnam.shp"


def load_polygon_points(shp_path):
    gdf = gpd.read_file(shp_path)
    points = []
    for geom in gdf.geometry:
        if geom is None: continue
        if geom.geom_type == 'Polygon':
            points.extend([[y, x] for x, y in geom.exterior.coords])
        elif geom.geom_type == 'MultiPolygon':
            for poly in geom.geoms:
                points.extend([[y, x] for x, y in poly.exterior.coords])
    return points


def reduce_to_convex_hull(points_geo):
    arr  = np.array(points_geo)
    hull = ConvexHull(arr)
    return arr[hull.vertices].tolist()


# Smallest enclosing circle (Welzl's algorithm)
# Note: implementation of Welzl's algorithm adapted with the help of Claude Opus 4.6

def circumcenter_3pts(P1, P2, P3):
    ax, ay = P1;  bx, by = P2;  cx, cy = P3
    D = 2*(ax*(by-cy) + bx*(cy-ay) + cx*(ay-by))
    ux = ((ax**2+ay**2)*(by-cy) + (bx**2+by**2)*(cy-ay) + (cx**2+cy**2)*(ay-by)) / D
    uy = ((ax**2+ay**2)*(cx-bx) + (bx**2+by**2)*(ax-cx) + (cx**2+cy**2)*(bx-ax)) / D
    center = np.array([ux, uy])
    return center, np.linalg.norm(center - np.array(P1))


def point_in_circle(pt, center, radius, tol=1e-10):
    return np.linalg.norm(np.array(pt) - center) <= radius + tol


def welzl(pts, boundary, n):
    if n == 0 or len(boundary) == 3:
        if len(boundary) == 0: return np.array([0.0, 0.0]), 0.0
        if len(boundary) == 1: return np.array(boundary[0], dtype=float), 0.0
        if len(boundary) == 2:
            c = (np.array(boundary[0]) + np.array(boundary[1])) / 2
            return c, np.linalg.norm(np.array(boundary[1]) - np.array(boundary[0])) / 2
        return circumcenter_3pts(*boundary)
    p = pts[n-1]
    center, radius = welzl(pts, boundary, n-1)
    if point_in_circle(p, center, radius):
        return center, radius
    return welzl(pts, boundary + [p], n-1)


def solve_sec(points_geo):
    pts = [list(p) for p in points_geo]
    center, radius = welzl(pts, [], len(pts))
    return center[0], center[1], radius


def run_stereo(shp_path):
    #Load polygon and reduce to convex hull
    points_geo = load_polygon_points(shp_path)
    hull_pts   = reduce_to_convex_hull(points_geo)

    #Pole K = center of smallest enclosing circle
    uk, vk, psi_j_deg = solve_sec(hull_pts)
    psi_j = radians(psi_j_deg)
    s_j   = pi/2 - psi_j

    #Multiplicative constant
    mju  = (2*cos(psi_j/2)**2) / (1 + cos(psi_j/2)**2)

    #Undistorted parallel
    psi0 = 2*acos(sqrt(mju))
    s0   = pi/2 - psi0

    #Scales
    m_pol = mju / cos(0)**2
    m1    = mju / cos(psi_j/2)**2
    m0    = mju / cos(psi0/2)**2

    #Distortions
    ny_pol = (m_pol-1)*1000
    ny1    = (m1-1)*1000
    ny0    = (m0-1)*1000

    print(ny1, ny_pol, ny0)
    return uk, vk, mju, psi_j, psi0, ny1, ny_pol, ny0


#Convert geographic coordinates to map coordinates
def uv_sd(u, v, uk_r, vk_r):
    s = np.arcsin(np.sin(u)*sin(uk_r) + np.cos(u)*cos(uk_r)*np.cos(vk_r-v))
    d = np.arctan2(np.cos(u)*np.sin(v-vk_r),
                   np.sin(u)*cos(uk_r) - np.cos(u)*sin(uk_r)*np.cos(vk_r-v))
    return s, d


def plot(shp_path, uk, vk, mju, psi_j, psi0, ny1, ny_pol, ny0):
    uk_r, vk_r = radians(uk), radians(vk)

    #Load country boundary
    gdf = gpd.read_file(shp_path)
    pts = []
    for geom in gdf.geometry:
        if geom is None: continue
        if geom.geom_type == 'Polygon':
            pts.extend([[y, x] for x, y in geom.exterior.coords])
        elif geom.geom_type == 'MultiPolygon':
            for poly in geom.geoms:
                pts.extend([[y, x] for x, y in poly.exterior.coords])

    #Country boundary coordinates
    u_b = np.array([radians(p[0]) for p in pts])
    v_b = np.array([radians(p[1]) for p in pts])
    s_b, d_b = uv_sd(u_b, v_b, uk_r, vk_r)

    psi_b = pi/2 - s_b
    rho_b = 2 * mju * np.tan(psi_b/2)
    bx =  rho_b * np.sin(d_b)
    by =  rho_b * np.cos(d_b)

    fig, ax = plt.subplots(figsize=(8, 8))

    #Graticule
    buf = radians(1)
    for v_g in np.arange(v_b.min()-buf, v_b.max()+2*buf, radians(1)):
        u_g = np.linspace(u_b.min()-buf, u_b.max()+buf, 200)
        sg, dg = uv_sd(u_g, np.full_like(u_g, v_g), uk_r, vk_r)
        pg = pi/2 - sg
        rg = 2*mju*np.tan(pg/2)
        ax.plot(rg*np.sin(dg), rg*np.cos(dg), 'k-', lw=0.3, alpha=0.4)
    for u_g in np.arange(u_b.min()-buf, u_b.max()+2*buf, radians(1)):
        v_g = np.linspace(v_b.min()-buf, v_b.max()+buf, 200)
        sg, dg = uv_sd(np.full_like(v_g, u_g), v_g, uk_r, vk_r)
        pg = pi/2 - sg
        rg = 2*mju*np.tan(pg/2)
        ax.plot(rg*np.sin(dg), rg*np.cos(dg), 'k-', lw=0.3, alpha=0.4)

    #Equideformates
    mg = 0.01
    X, Y = np.meshgrid(np.linspace(bx.min()-mg, bx.max()+mg, 400),
                       np.linspace(by.min()-mg, by.max()+mg, 400))
    RHO = np.sqrt(X**2 + Y**2)
    PSI = 2*np.arctan(RHO / (2*mju))
    NU  = (mju / np.cos(PSI/2)**2 - 1) * 1000

    ny_bnd = abs(ny1)
    step   = round(ny_bnd/4, 2)
    levels = np.round(np.arange(round(ny_pol, 2)-step, ny_bnd+step*1.5, step), 2)

    cs = ax.contour(X, Y, NU, levels=levels, colors='red', linewidths=0.8)
    ax.clabel(cs, fmt='%.2f', fontsize=7)

    # Country boundary
    ax.plot(bx, by, 'b-', lw=1.5)
    ax.set_aspect('equal')
    plt.show()


vals = run_stereo(SHP_PATH)
plot(SHP_PATH, *vals)
