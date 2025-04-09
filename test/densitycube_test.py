import unittest
import tempfile
from numpy.ma.testutils import assert_equal
import os

from densityCube import DensityCube


class TestDensityCube(unittest.TestCase):
    def setUp(self):
        # Create a simple fake cube file for testing
        self.cube_content = ''' Title Line1
        Title Line2
          1  0.0    0.0  0.0
          3    1.0    0.0  0.0
          3    0.0    1.0  0.0
          3    0.0    0.0  1.0
          1    1.0    0.0  0.0 
         8.48011E-21  1.74768E-20  3.54520E-20  7.08371E-20  1.39548E-19  2.71342E-19
         5.21439E-19  9.91787E-19  1.86997E-18  3.50034E-18  6.51351E-18  1.20601E-17
         2.22265E-17  4.07605E-17  7.43042E-17  1.34425E-16  2.40842E-16  4.26328E-16
         7.43817E-16  1.27610E-15  2.14822E-15  3.54184E-15  5.71005E-15  8.98933E-15
         1.38043E-14  2.06590E-14  3.01094E-14'''

        # Create a temporary cube file
        self.temp_file = tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.cube')
        self.temp_file.write(self.cube_content)
        self.temp_file.close()

        # Create the cube object here so it's available to all tests
        self.cube = DensityCube(self.temp_file.name, d_range=(0.0, 10.0))

    def testGridShape(self):
        actual_grid_shape = self.cube.grid_shape
        desired_grid_shape = [3, 3, 3]
        assert_equal(actual_grid_shape, desired_grid_shape)

    def test_voxel_dimensions(self):
        assert_equal(self.cube.xVoxel, 3)
        assert_equal(self.cube.yVoxel, 3)
        assert_equal(self.cube.zVoxel, 3)

    def test_atom_info(self):
        assert_equal(self.cube.atom_count, 1)
        assert_equal(self.cube.atom_x_origin, 0.0)
        assert_equal(self.cube.atom_y_origin, 0.0)
        assert_equal(self.cube.atom_z_origin, 0.0)

    def test_voxel_info(self):
        assert_equal(self.cube.xVoxel, 3)
        assert_equal(self.cube.xVoxelWidth, 1)
        assert_equal(self.cube.yVoxel, 3)
        assert_equal(self.cube.yVoxelWidth, 1)
        assert_equal(self.cube.zVoxel, 3)
        assert_equal(self.cube.zVoxelWidth, 1)

    def test_dfdx(self):
        # Test derivative at point (1,1,1)
        result = self.cube.dfdx(1, 1, 1)
        expected = self.cube.grid_data[1, 1, 2] - self.cube.grid_data[1, 1, 1]
        self.assertEqual(result, expected)

    def test_d2f_dx2(self):
        # At the center, it should compute correctly
        result = self.cube.d2f_dx2(1, 1, 1)
        expected = (self.cube.grid_data[1, 1, 2] - 2 * self.cube.grid_data[1, 1, 1] +
                    self.cube.grid_data[1, 1, 0])
        self.assertEqual(result, expected)

    def test_integrate_density(self):
        self.cube.d_range = [0, 1000]
        # For a known cube of values 0..26, total = sum * voxel_volume (1)
        calculated_density = self.cube.integrate_density()

        # Example: expected integrated density value
        expected_density = 8.79e-14

        # Assert that the values are almost equal with a tolerance of 1e-15
        self.assertAlmostEqual(calculated_density, expected_density, places=15)

    def test_hessian_output_shape(self):
        hessians = self.cube.calculate_hessian_for_all_points()
        self.assertEqual(hessians.shape[0], 27)  # 3*3*3 grid points
        self.assertEqual(hessians.shape[1:], (3, 3))

    def tearDown(self):
        os.unlink(self.temp_file.name)  # Delete the temp file after tests


if __name__ == '__main__':
    unittest.main()
