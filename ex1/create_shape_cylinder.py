import numpy as np


def generate_cylinder(radius, height, n_points):
    angles = np.linspace(0.0, 2.0 * np.pi, n_points, endpoint=False)

    x = radius * np.cos(angles)
    y = radius * np.sin(angles)

    z_top = +0.5 * height
    z_bottom = -0.5 * height

    top = np.column_stack((x, y, np.full(n_points, z_top)))
    bottom = np.column_stack((x, y, np.full(n_points, z_bottom)))

    vertices = np.vstack((top, bottom))
    return vertices


def generate_faces(n_points):
    faces = []

    # Faces latérales (quadrilatères)
    for i in range(n_points):
        j = (i + 1) % n_points
        faces.append([i, j, j + n_points, i + n_points])

    # Face supérieure
    faces.append(list(range(n_points)))

    # Face inférieure
    faces.append(list(range(n_points, 2 * n_points)))

    return faces


def generate_edges(faces):
    edges = set()
    for face in faces:
        for i in range(len(face)):
            a = face[i]
            b = face[(i + 1) % len(face)]
            edges.add(tuple(sorted((a, b))))
    return sorted(edges)


def write_shape_file(
    filename,
    vertices,
    edges,
    faces,
    name="Cylinder",
    radius=0.1,
):
    nv = len(vertices)
    ne = len(edges)
    nf = len(faces)

    # OBB (axis-aligned)
    min_corner = vertices.min(axis=0)
    max_corner = vertices.max(axis=0)
    extent = 0.5 * (max_corner - min_corner)
    center = 0.5 * (max_corner + min_corner)

    with open(filename, "w") as f:
        f.write("<\n")
        f.write(f"name {name}\n")
        f.write(f"radius {radius}\n")
        f.write("preCompDone y\n")

        # Sommets
        f.write(f"nv {nv}\n")
        for v in vertices:
            f.write(f"{v[0]} {v[1]} {v[2]}\n")

        # Edges
        f.write(f"ne {ne}\n")
        for e in edges:
            f.write(f"{e[0]} {e[1]}\n")

        # Faces
        f.write(f"nf {nf}\n")
        for face in faces:
            f.write(f"{len(face)} " + " ".join(str(i) for i in face) + "\n")

        # OBB
        f.write(
            f"obb.extent {extent[0]} {extent[1]} {extent[2]}\n"
        )
        f.write("obb.e1 1.0 0.0 0.0\n")
        f.write("obb.e2 0.0 1.0 0.0\n")
        f.write("obb.e3 0.0 0.0 1.0\n")
        f.write(
            f"obb.center {center[0]} {center[1]} {center[2]}\n"
        )

        # Pose
        f.write("position 0.0 0.0 0.0\n")
        f.write("orientation 1.0 0.0 0.0 0.0\n")

        # Volume analytique du cylindre
        volume = np.pi * radius**2 * (max_corner[2] - min_corner[2])
        f.write(f"volume {volume}\n")

        # Inertie massique (approx analytique, axe z)
        h = max_corner[2] - min_corner[2]
        i_x = (1.0 / 12.0) * (3 * radius**2 + h**2)
        i_z = 0.5 * radius**2
        f.write(f"I/m {i_x} {i_x} {i_z}\n")

        f.write(">\n")


# ----------------------
# Paramètres
# ----------------------
radius = 0.5
height = 10
n_points = 32

# ----------------------
# Génération
# ----------------------
vertices = generate_cylinder(radius, height, n_points)
faces = generate_faces(n_points)
edges = generate_edges(faces)

# ----------------------
# Écriture fichier shape
# ----------------------
write_shape_file(
    filename="cylinder.shp",
    vertices=vertices,
    edges=edges,
    faces=faces,
    name="Cylinder",
    radius=radius,
)

