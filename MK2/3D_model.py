import numpy as np
import pyvista as pv
import os


R = 1.0

# Path to the folder
CONT_PATH = r"C:\Users\kadle\OneDrive - CUNI\A-Dokumenty\.Uni\_4.MGR\Matkarto\2_uloha\Kontinenty"

# Continent files
CONT_FILES = ["afr_mad.txt", "afr.txt", "amer.txt", "antar.txt", "austr.txt", "eur.txt", "greenl.txt", "newzel1.txt", "newzel2.txt", "tasm.txt"]


# Dodecahedron vertex latitudes
us1 = np.radians(52.6226)
us2 = np.radians(10.8123)
uj2 = np.radians(-10.8123)
uj1 = np.radians(-52.6226)

# Helper functions for conversion and object creation
def latlon_to_xyz(lat, lon, r=R):
    # Convert geographic coordinates to Cartesian
    x = r * np.cos(lat) * np.cos(lon)
    y = r * np.cos(lat) * np.sin(lon)
    z = r * np.sin(lat)
    if hasattr(lat, '__len__'):
        return np.column_stack([x, y, z])
    return np.array([x, y, z])


def load_continent(filepath):
    """Loads continent points (lat lon in degrees)."""
    try:
        data = np.loadtxt(filepath)
        return np.radians(data[:, 0]), np.radians(data[:, 1])
    except Exception as e:
        print(f"  Failed to load {filepath}: {e}")
        return None, None


def make_polyline(points):
    """Creates a PyVista line from points."""
    n = len(points)
    if n < 2:
        return None
    lines = np.zeros(n + 1, dtype=int)
    lines[0] = n
    lines[1:] = np.arange(n)
    return pv.PolyData(points, lines=lines)

# Dodecahedron definition
vertex_names = ['A', 'B', 'C', 'D', 'E',
                'F', 'H', 'J', 'L', 'N',
                'G', 'I', 'K', 'M', 'O',
                'P', 'Q', 'R', 'S', 'T']

vertex_latlon = {
    'A': (uj1, 0),     'B': (uj1, 72),    'C': (uj1, 144),
    'D': (uj1, 216),   'E': (uj1, 288),
    'F': (uj2, 0),     'H': (uj2, 72),    'J': (uj2, 144),
    'L': (uj2, 216),   'N': (uj2, 288),
    'G': (us2, 36),    'I': (us2, 108),   'K': (us2, 180),
    'M': (us2, 252),   'O': (us2, 324),
    'P': (us1, 36),    'Q': (us1, 108),   'R': (us1, 180),
    'S': (us1, 252),   'T': (us1, 324),
}

vertex_xyz = {}
vertex_array = []
name_to_idx = {}
for i, name in enumerate(vertex_names):
    lat, lon_deg = vertex_latlon[name]
    xyz = latlon_to_xyz(lat, np.radians(lon_deg))
    vertex_xyz[name] = xyz
    vertex_array.append(xyz)
    name_to_idx[name] = i
vertex_array = np.array(vertex_array)

faces_names = [
    ['A', 'F', 'G', 'H', 'B'],   # Face 1  (K6)
    ['F', 'O', 'T', 'P', 'G'],   # Face 2  (K1)
    ['G', 'P', 'Q', 'I', 'H'],   # Face 3  (K2)
    ['I', 'Q', 'R', 'K', 'J'],   # Face 4  (K3)
    ['K', 'R', 'S', 'M', 'L'],   # Face 5  (K4)
    ['M', 'S', 'T', 'O', 'N'],   # Face 6  (K5)
    ['B', 'H', 'I', 'J', 'C'],   # Face 7  (K7)
    ['C', 'J', 'K', 'L', 'D'],   # Face 8  (K8)
    ['D', 'L', 'M', 'N', 'E'],   # Face 9  (K9)
    ['E', 'N', 'O', 'F', 'A'],   # Face 10 (K10)
    ['P', 'Q', 'R', 'S', 'T'],   # Face 11 (K11 - N pole)
    ['A', 'B', 'C', 'D', 'E'],   # Face 12 (K12 - S pole)
]

faces_flat = []
for face in faces_names:
    indices = [name_to_idx[v] for v in face]
    faces_flat.append(len(indices))
    faces_flat.extend(indices)
faces_flat = np.array(faces_flat)

# Dodecahedron faces
def compute_face_planes():
   
    face_planes = []
    for face in faces_names:
        verts = np.array([vertex_xyz[v] for v in face])
        centroid = np.mean(verts, axis=0)
        normal = centroid / np.linalg.norm(centroid)
        dist = np.dot(normal, verts[0])  # plane distance from the center
        face_planes.append({
            'normal': normal,
            'dist': dist,
            'centroid': centroid,
            'vertices': verts,
        })
    return face_planes


def find_faces(xyz_sphere, face_planes):
    
    normals = np.array([fp['normal'] for fp in face_planes])  # (12, 3)
    dots = xyz_sphere @ normals.T  # (N, 12)
    return np.argmax(dots, axis=1)


def project_points_to_dodeca(xyz_sphere, face_planes):
   
    face_ids = find_faces(xyz_sphere, face_planes)
    
    normals = np.array([fp['normal'] for fp in face_planes])[face_ids]
    dists = np.array([fp['dist'] for fp in face_planes])[face_ids]

    p_hat = xyz_sphere / np.linalg.norm(xyz_sphere, axis=1, keepdims=True)
    denom = np.sum(normals * p_hat, axis=1)
    
    # Depth buffer
    t = (dists * 1.002) / denom
    
    projected = p_hat * t[:, np.newaxis]

    return projected, face_ids


def split_into_segments(projected, face_ids, threshold=0.15):
   
    segments = []
    start = 0
    for i in range(1, len(projected)):
        do_break = False
        # Face change
        if face_ids[i] != face_ids[i - 1]:
            do_break = True
        # Large jump
        elif np.linalg.norm(projected[i] - projected[i - 1]) > threshold:
            do_break = True

        if do_break:
            if i - start >= 2:
                segments.append(projected[start:i])
            start = i

    if len(projected) - start >= 2:
        segments.append(projected[start:])

    return segments


# Creation of objects for visualization
def create_dodecahedron_mesh():
    """Dodecahedron mesh (opaque faces)."""
    return pv.PolyData(vertex_array, faces=faces_flat)


def create_dodecahedron_edges():
    """Dodecahedron edges as lines."""
    edges_set = set()
    points = []
    lines = []
    offset = 0

    for face in faces_names:
        n = len(face)
        for i in range(n):
            v1, v2 = face[i], face[(i + 1) % n]
            edge = tuple(sorted([v1, v2]))
            if edge not in edges_set:
                edges_set.add(edge)
                points.append(vertex_xyz[v1])
                points.append(vertex_xyz[v2])
                lines.extend([2, offset, offset + 1])
                offset += 2

    return pv.PolyData(np.array(points), lines=np.array(lines))


def create_graticule_on_dodeca(face_planes, du_deg=10, dv_deg=10, step_deg=1):
    """Geographic graticule projected onto the faces of the dodecahedron."""
    step = np.radians(step_deg)
    all_polylines = []

    # Meridians
    for v_deg in range(-180, 180, dv_deg):
        v = np.radians(v_deg)
        u = np.arange(-np.pi / 2, np.pi / 2 + step, step)
        xyz_sphere = latlon_to_xyz(u, np.full_like(u, v))
        projected, face_ids = project_points_to_dodeca(xyz_sphere, face_planes)
        segs = split_into_segments(projected, face_ids, threshold=0.2)
        for seg in segs:
            poly = make_polyline(seg)
            if poly is not None:
                all_polylines.append(poly)

    # Parallels
    for u_deg in range(-80, 90, du_deg):
        u = np.radians(u_deg)
        v = np.arange(-np.pi, np.pi + step, step)
        xyz_sphere = latlon_to_xyz(np.full_like(v, u), v)
        projected, face_ids = project_points_to_dodeca(xyz_sphere, face_planes)
        segs = split_into_segments(projected, face_ids, threshold=0.2)
        for seg in segs:
            poly = make_polyline(seg)
            if poly is not None:
                all_polylines.append(poly)

    if all_polylines:
        return all_polylines[0].merge(all_polylines[1:])
    return None


def create_continents_on_dodeca(cont_path, cont_files, face_planes):
    """Continent outlines projected onto the faces of the dodecahedron."""
    all_polylines = []

    for filename in cont_files:
        filepath = os.path.join(cont_path, filename)
        if not os.path.exists(filepath):
            print(f"  File not found: {filepath}")
            continue

        lat, lon = load_continent(filepath)
        if lat is None:
            continue

        xyz_sphere = latlon_to_xyz(lat, lon)
        projected, face_ids = project_points_to_dodeca(xyz_sphere, face_planes)
        segs = split_into_segments(projected, face_ids, threshold=0.15)

        for seg in segs:
            poly = make_polyline(seg)
            if poly is not None:
                all_polylines.append(poly)

        name = os.path.splitext(filename)[0]
        print(f"  {name}: {len(lat)} points -> {len(segs)} segments")

    if all_polylines:
        return all_polylines[0].merge(all_polylines[1:])
    return None


# Main script
if __name__ == '__main__':

    print("=" * 60)
    print("3D polyhedral globe - dodecahedron (PyVista)")
    print("Continents and graticule on dodecahedron faces")
    print("=" * 60)

    # Calculation of face planes
    face_planes = compute_face_planes()

    print("\n1. Creating dodecahedron...")
    dodeca_mesh = create_dodecahedron_mesh()
    dodeca_edges = create_dodecahedron_edges()

    print("2. Projecting geographic graticule onto faces...")
    graticule = create_graticule_on_dodeca(face_planes, du_deg=10, dv_deg=10)

    print("3. Projecting continents onto faces...")
    continents = create_continents_on_dodeca(CONT_PATH, CONT_FILES, face_planes)

    # Visualization
    print("\n4. Launching visualization...")

    plotter = pv.Plotter(window_size=[1400, 1000])
    plotter.set_background('white')

    # Dodecahedron faces
    plotter.add_mesh(dodeca_mesh, color='white', opacity=1.0,
                     show_edges=False, smooth_shading=False)

    # Dodecahedron edges (red)
    plotter.add_mesh(dodeca_edges, color='red', line_width=3)


    # Geographic graticule on faces (black)
    if graticule is not None:
        plotter.add_mesh(graticule, color='black', line_width=0.8,
                         opacity=0.6)

    # Continents on faces (blue)
    if continents is not None:
        plotter.add_mesh(continents, color='blue', line_width=1.8)


    plotter.add_title("Polyhedral globe - dodecahedron", font_size=14)
    plotter.camera_position = 'xz'

    # Export
    plotter.show(screenshot='dvanactisten_3d.png')
    dodeca_mesh.save('dvanactisten.stl')
    dodeca_mesh.save('dvanactisten_mesh.vtk')
    print("\n5. Exported:")
    print("   - dvanactisten_3d.png")
    print("   - dvanactisten_mesh.vtk")
    print("\nDone!")
