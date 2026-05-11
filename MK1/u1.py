from math import *
from uvtosd import *
import requests
from io import BytesIO
from PIL import Image
import matplotlib.pyplot as plt

def WGSToJTSK(phi_WGS, la_WGS):
    #WGS84 parameters
    a_WGS = 6378137.00
    b_WGS = 6356752.3142
    
    e2_WGS = (a_WGS*a_WGS - b_WGS*b_WGS)/(a_WGS*a_WGS)
    W_WGS = sqrt(1-e2_WGS*(sin(phi_WGS))**2)
    N_WGS = a_WGS/W_WGS
    
    #XYZ coordinates, WGS 84
    X_WGS = N_WGS * (cos(phi_WGS) * cos(la_WGS))
    Y_WGS = N_WGS * (cos(phi_WGS) * sin(la_WGS))
    Z_WGS = N_WGS * (1 - e2_WGS) * sin(phi_WGS)
    
    #Helmert transformation, parameters
    om_x = 4.9984 / 3600 * pi / 180
    om_y = 1.5867 / 3600 * pi / 180
    om_z = 5.2611 / 3600 * pi / 180
    m = 1 - 3.5623e-06
    dlt_x = -570.8285
    dlt_y = -85.6769
    dlt_z = -462.8420
    
    #Helmert transformation, Bessel ellipsoid
    X_Bes = m * (X_WGS + Y_WGS * om_z - om_y * Z_WGS) + dlt_x
    Y_Bes = m * (-om_z * X_WGS + Y_WGS + om_x * Z_WGS) + dlt_y
    Z_Bes = m * (om_y * X_WGS - om_x * Y_WGS + Z_WGS) + dlt_z
    
    #Bessel parameters
    a_Bes = 6377397.155
    b_Bes = 6356078.963
    e2_Bes = (a_Bes*a_Bes - b_Bes*b_Bes)/(a_Bes*a_Bes)
    
    #Phi, lam, Bessel
    la_Bes = atan2(Y_Bes, X_Bes)
    tan_phi_Bes = Z_Bes / ((1 - e2_Bes) * sqrt(X_Bes**2 + Y_Bes**2))
    phi_Bes = atan(tan_phi_Bes)
    
    #Shift to Ferro
    la_Ferro = la_Bes + (17 + 2/3) * pi / 180
    
    #Gauss conformal projection, parameters
    phi0 = 49.5 * pi / 180
    alpha = sqrt(1 + e2_Bes * (cos(phi0))**4 / (1 - e2_Bes))
    u0 = asin(sin(phi0) / alpha)
    
    kn = (tan(phi0/2+pi/4)**alpha * ((1-sqrt(e2_Bes)*sin(phi0))/(1+sqrt(e2_Bes)*sin(phi0)))**(alpha*sqrt(e2_Bes)/2))
    kd = tan(u0/2+pi/4)
    k = kn / kd
    
    R = (a_Bes*sqrt(1-e2_Bes)) / (1-e2_Bes*(sin(phi0)**2))
    
    #Gauss conformal projection
    u = 2*(atan(1/k*(tan(phi_Bes/2+pi/4)*((1-sqrt(e2_Bes)*sin(phi_Bes))/(1+sqrt(e2_Bes)*sin(phi_Bes)))**(sqrt(e2_Bes)/2))**alpha))-pi/2
    v = alpha * la_Ferro
    
    #Cartographic pole
    uk = (59+(42/60)+(42.6969/3600)) * (pi/180)
    vk = (42+(31/60)+(31.41725/3600)) * (pi/180)
    
    #Conversion (u, v) -> (s, d)
    s, d = uvTosd(u, v, uk, vk)
    
    #LCC
    s0 = 78.5 * pi/180
    rho0 = R * 1/tan(s0) * 0.9999
    c = sin(s0)
    
    rho = rho0 * ((tan(s0/2+pi/4)) / (tan(s/2+pi/4)))**c
    eps = c * d
    
    # (rho, eps) -> (y, x)
    y_jtsk = rho * sin(eps)
    x_jtsk = rho * cos(eps)
    
    # Local linear scale
    m_scale = (c * rho) / (R * cos(s))
    
    print(f"Y_JTSK = {y_jtsk:.3f} m")
    print(f"X_JTSK = {x_jtsk:.3f} m")
    print(f"Length distortion m-1 = {m_scale - 1:.3e}")
    
    return y_jtsk, x_jtsk, m_scale


#Disclaimer: Visualization and OSM tiles related code was created with the help of Claude Opus 4.6

#Geographic coordinates (deg) -> OSM tile index (x, y)
def _deg2tile(lat_deg, lon_deg, zoom):
    lat_r = radians(lat_deg)
    n = 2 ** zoom
    x = int((lon_deg + 180.0) / 360.0 * n)
    y = int((1.0 - log(tan(lat_r) + 1.0/cos(lat_r)) / pi) / 2.0 * n)
    return x, y

#OSM tile index
def _tile2deg(x, y, zoom):
    n = 2 ** zoom
    lon = x / n * 360.0 - 180.0
    lat = degrees(atan(sinh(pi * (1 - 2*y/n))))
    return lat, lon

#Download and adjust OSM tiles
def fetch_osm_tiles(minlat, minlon, maxlat, maxlon, zoom=17):
    headers = {"User-Agent": "Mozilla/5.0 (geodesy student project)"}
    #Tile index range covering the bounding box
    tx_min, ty_max = _deg2tile(minlat, minlon, zoom)
    tx_max, ty_min = _deg2tile(maxlat, maxlon, zoom)
    tile_size = 256
    cols = tx_max - tx_min + 1
    rows = ty_max - ty_min + 1
    #Empty canvas for the mosaic
    mosaic = Image.new("RGB", (cols * tile_size, rows * tile_size))
    for row, ty in enumerate(range(ty_min, ty_max + 1)):
        for col, tx in enumerate(range(tx_min, tx_max + 1)):
            url = f"https://tile.openstreetmap.org/{zoom}/{tx}/{ty}.png"
            try:
                r = requests.get(url, headers=headers, timeout=15)
                r.raise_for_status()
                tile = Image.open(BytesIO(r.content)).convert("RGB")
            except Exception as e:
                print(f"  Tile {tx}/{ty} failed: {e}")
                tile = Image.new("RGB", (tile_size, tile_size), (220, 220, 220))
            mosaic.paste(tile, (col * tile_size, row * tile_size))
    #Geographic extent of the assembled mosaic
    top_lat, left_lon  = _tile2deg(tx_min,     ty_min,     zoom)
    bot_lat, right_lon = _tile2deg(tx_max + 1, ty_max + 1, zoom)
    return mosaic, (left_lon, right_lon, bot_lat, top_lat)


def visualize(phi_deg, la_deg, y_jtsk, x_jtsk, m_scale):
    #Map extent around the point
    MARGIN = 0.0008
    minlat = phi_deg - MARGIN
    maxlat = phi_deg + MARGIN
    minlon = la_deg  - MARGIN * 1.5
    maxlon = la_deg  + MARGIN * 1.5

    fig, ax = plt.subplots(figsize=(7, 7))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("#e8e0d8")

    #OSM background
    print("Loading OSM tiles...")
    img, ext = fetch_osm_tiles(minlat, minlon, maxlat, maxlon, zoom=19)
    ax.imshow(img, extent=ext, origin="upper", aspect="auto", zorder=0)

    #WGS84 input point
    ax.plot(la_deg, phi_deg, "o", ms=11, color="white",
            mec="black", mew=1.5, zorder=6, label="WGS84")
    #S-JTSK output point – blue triangle (same physical location, drawn on top)
    ax.plot(la_deg, phi_deg, "^", ms=8, color="#1B47C2",
            mec="white", mew=1.0, zorder=7, label="S-JTSK")

    #WGS84 label
    ax.annotate(f"WGS84\nphi = {phi_deg:.6f} deg\nla  = {la_deg:.6f} deg",
                xy=(la_deg, phi_deg), xytext=(-90, 30),
                textcoords="offset points", fontsize=8.5, ha="center",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                          edgecolor="#333", lw=1.2, alpha=0.95),
                arrowprops=dict(arrowstyle="->", color="#333", lw=1.0), zorder=10)

    #S-JTSK label
    ax.annotate(f"S-JTSK\nY = {y_jtsk:.3f} m\nX = {x_jtsk:.3f} m\nm-1 = {m_scale-1:.3e}",
                xy=(la_deg, phi_deg), xytext=(90, -45),
                textcoords="offset points", fontsize=8.5, ha="center",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="#DDEEFF",
                          edgecolor="#1565C0", lw=1.4, alpha=0.95),
                arrowprops=dict(arrowstyle="->", color="#1565C0", lw=1.0), zorder=10)

    #Axes, legend, grid
    ax.legend(loc="upper right", fontsize=8.5, framealpha=0.9,
              edgecolor="#aaa", facecolor="white")
    ax.set_xlim(minlon, maxlon)
    ax.set_ylim(minlat, maxlat)
    ax.set_xlabel("Longitude [deg]", fontsize=9)
    ax.set_ylabel("Latitude [deg]", fontsize=9)
    ax.set_title("WGS-84 to S-JTSK", fontsize=10, pad=8)
    ax.tick_params(labelsize=7.5)
    ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.4f'))
    ax.xaxis.set_major_formatter(plt.FormatStrFormatter('%.4f'))
    ax.grid(True, color="#bbb", lw=0.4, ls="--", alpha=0.6)

    plt.show()

#P1
#phi_WGS = 50.064270 * pi/180
#la_WGS = 14.419012 * pi/180

#P2
phi_WGS = 50.050572 * pi/180
la_WGS = 14.384470 * pi/180

#Compute S-JTSK coordinates and length distortion
y_jtsk, x_jtsk, m_scale = WGSToJTSK(phi_WGS, la_WGS)

#Convert input to degrees
phi_deg = phi_WGS * 180 / pi
la_deg  = la_WGS  * 180 / pi

visualize(phi_deg, la_deg, y_jtsk, x_jtsk, m_scale)