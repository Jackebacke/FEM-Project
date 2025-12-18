import os
import matplotlib.pyplot as plt
import numpy as np
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader

FIG_FOLDER = "C/Figures"

def data_from_vtu(filename):
    # Read the first VTU file
    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()

    mesh = reader.GetOutput()
    print(f"Loaded mesh with {mesh.GetNumberOfPoints()} points and {mesh.GetNumberOfCells()} cells")

    # Get the vector field
    vector_array = mesh.GetPointData().GetVectors()
    if vector_array is None:
        # Try by name
        vector_array = mesh.GetPointData().GetArray("f_12")

    if vector_array is None:
        print("Error: Could not find vector data")
        exit(1)

    print(f"Vector array: {vector_array.GetName()} with {vector_array.GetNumberOfComponents()} components")

    # Extract components into numpy arrays
    n_points = mesh.GetNumberOfPoints()
    components = {}

    for comp_idx in range(vector_array.GetNumberOfComponents()):
        comp_data = []
        for point_idx in range(n_points):
            value = vector_array.GetValue(point_idx * vector_array.GetNumberOfComponents() + comp_idx)
            comp_data.append(value)
        components[comp_idx] = np.array(comp_data)

    print(f"Component 0 (u): min={components[0].min():.3e}, max={components[0].max():.3e}")
    print(f"Component 1 (v): min={components[1].min():.3e}, max={components[1].max():.3e}")
    print(f"Component 2 (w): min={components[2].min():.3e}, max={components[2].max():.3e}")

    # Get coordinates
    coords = []
    for i in range(n_points):
        coords.append(mesh.GetPoint(i))
    coords = np.array(coords)

    print(f"Coordinates shape: {coords.shape}")
    
    # Extract triangulation from mesh cells
    triangles = []
    for cell_id in range(mesh.GetNumberOfCells()):
        cell = mesh.GetCell(cell_id)
        if cell.GetNumberOfPoints() == 3:
            point_ids = [cell.GetPointId(i) for i in range(3)]
            triangles.append(point_ids)
    triangles = np.array(triangles)
    print(f"Number of triangles: {len(triangles)}")
    
    return coords, components[0], components[1], components[2], triangles

def plot_solution(coords, u, v, w, triangles, output_filename, title="Population densities at t=0"): 
    # Extract x, y coordinates
    x = coords[:, 0]
    y = coords[:, 1]
    
    # Create 3 subplots with colored mesh (ParaView 2D style)
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    titles = ["u(x,t)", "v(x,t)", "w(x,t)"]
    components_list = [u, v, w]

    for idx in range(3):
        ax = axes[idx]
        component = components_list[idx]
        
        # Average component values at triangle centers for coloring
        triangle_values = component[triangles].mean(axis=1)
        
        # Create colored mesh with triangle edges
        tripcolor = ax.tripcolor(x, y, triangles, facecolors=triangle_values, cmap='cool', 
                                edgecolors='black', linewidths=0.1, shading='flat')
        
        ax.set_title(titles[idx], fontsize=12, fontweight='bold')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_aspect('equal')
        
        # Add colorbar
        cbar = plt.colorbar(tripcolor, ax=ax)
        cbar.set_label('Value', fontsize=10)
    plt.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    # Save figure
    os.makedirs(FIG_FOLDER, exist_ok=True)
    output_path = os.path.join(FIG_FOLDER, output_filename)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved figure to {output_path}")
    plt.close()
 
# Setup
data="sverige"
id = "000008"


if data == "sverige":
    data_filename = f"{data}_solution{id}.vtu"
else:
    data_filename = f"{data}{id}.vtu"
    
data_folder = f"C/Solutions/{data}"
vtu_path = os.path.join(data_folder, data_filename)

coords, u, v, w, triangles = data_from_vtu(vtu_path)
plot_solution(coords, u, v, w, triangles, output_filename=f"solution.png", title=f"Population Densities for {data} at t=0")
