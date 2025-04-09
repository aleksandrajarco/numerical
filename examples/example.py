# example_usage.py

from densityCube import DensityCube
import os
# Specify the path to your example cube file
cube_file_path = os.path.join(os.path.dirname(__file__), 'density.cube')


# Initialize the DensityCube object with the cube file and the desired range
cube = DensityCube(cube_file_path, d_range=(0.0, 10.0))

# Print basic information about the cube
print("Grid Shape:", cube.grid_shape)
print("Voxel Dimensions (x, y, z):", cube.xVoxel, cube.yVoxel, cube.zVoxel)
print("Voxel Widths (x, y, z):", cube.xVoxelWidth, cube.yVoxelWidth, cube.zVoxelWidth)

# Print atom info
print("Atom Count:", cube.atom_count)
print("Atom Position (x, y, z):", cube.atom_x_origin, cube.atom_y_origin, cube.atom_z_origin)

# Compute and print the density integral within the specified range
integrated_density = cube.integrate_density()
print(f"Integrated Density: {integrated_density:.15e}")

# Calculate the derivative in the x-direction at point (1, 1, 1)
dfdx = cube.dfdx(1, 1, 1)
print(f"df/dx at (1, 1, 1): {dfdx}")

# Calculate the second derivative in the x-direction at point (1, 1, 1)
d2fdx2 = cube.d2f_dx2(1, 1, 1)
print(f"d^2f/dx^2 at (1, 1, 1): {d2fdx2}")

# Calculate and print the Hessian matrix at all grid points
hessians = cube.calculate_hessian_for_all_points()
print("Shape of Hessian matrix at all grid points:", hessians.shape)

