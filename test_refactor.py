import unittest

from refactor import *

class TestBox(unittest.TestCase):
    
    def setUp(self):
        self.box = BoxWithCartesianBoundary() # default setup

    def test_is_outside_on_boundary(self):
        pt = self.box.x_limits[1], self.box.y_limits[0]
        self.assertTrue( np.allclose( [self.box.is_outside(pt)], [0.0] )  )

    def test_is_boundary_point(self):
        pt1 = [-0.01, 1.019]
        pt2 = [0.005, 1.02]

        check1 = self.box.is_boundary_point(pt1) 
        check2 = self.box.is_boundary_point(pt2)

        self.assertTrue( check1 and check2 )

    def test_get_grid_block_vertices_interior(self):
        pt = [0.006, 1.006]
        nbrs = self.box.get_grid_block_vertices(pt)

        expected_nbrs = [
            [0.00, 1.00],
            [0.01, 1.00],
            [0.00, 1.01],
            [0.01, 1.01]
        ]

        check = True 
        for expected in expected_nbrs:
            flag = False
            for nbr in nbrs:
                if np.allclose(expected, nbr):
                    flag = True 
            if flag is False:
                check = False 
                break 
        
        self.assertTrue(check)

    def test_get_grid_block_vertices_boundary(self):
        # pt = [-0.09, 1.02] --> behavior on the boundary is currently not specified.
        pass 


    def test_nearest_boundary_point(self):
        pt = [-0.01, 1.0196] # don't test on a boundary point
        bdry_pt = self.box.hit_boundary(pt)

        expected = [-0.01, 1.020]
        self.assertTrue( np.allclose(expected, bdry_pt) )

    def test_nearest_boundary_point_exterior(self):
        pt = [-0.011, 1.0196] # don't test on a boundary point
        bdry_pt = self.box.hit_boundary(pt)

        expected = None
        self.assertTrue( expected == bdry_pt )

    def test_nearest_boundary_point_snap_y(self):
        pt = np.array([-0.009, 1.0186]) # don't test on a boundary point
        bdry_pt = self.box.nearest_boundary_point(pt)

        expected = np.array([-0.01, 1.019])

        #print(np.linalg.norm(expected - pt))
        #print(np.linalg.norm(bdry_pt - pt))
        
        self.assertTrue( np.allclose(expected, bdry_pt) )

    def test_nearest_boundary_point_not_in_bdry(self):
        pt = np.array([-0.011, 1.019]) # an exterior point; should snap to x axis
        expected = np.array([-0.01, 1.019] )
        bdry_pt = self.box.nearest_boundary_point(pt)
        self.assertTrue( np.allclose(expected, bdry_pt) )

    

if __name__ == '__main__':
    unittest.main()