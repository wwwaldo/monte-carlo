import numpy as np
import matplotlib.pyplot as plt

class BoxWithCartesianBoundary():
    '''
    Params for this test box.
    '''
    def __init__(self):
        self.dh = 0.01
        self.dboundary = 0.001

        # The x, y bounds on the box.
        self.x_limits = [-0.015, 0.015]
        self.y_limits = [0.975, 1.025]

        # TODO: make generate_test_case depend on x,y limits
        self.bdry_X , self.bdry_Y = BoxWithCartesianBoundary.generate_test_case(
            self.dh, self.dboundary
        )
    
    '''
    "indexing" functions.
    '''
    def get_grid_block_vertices(self, pt):
        x, y = pt
        h = self.dh
        basex, basey = h * np.floor( x / h ), h * np.floor( y / h )
        return list(map(np.array, [ [basex, basey], [basex + h, basey],
            [basex, basey + h], [basex + h, basey + h]]))

    def get_index(self, pt):
        X, Y = self.bdry_X, self.bdry_Y
        temp_zeros = np.zeros(len(X))
        
        # Throws an obscure error if a non-bdry point is passed in
        index = np.intersect1d( np.where(np.isclose(X - pt[0], temp_zeros)) , 
                            np.where(np.isclose( Y - pt[1], temp_zeros)) 
                            )[0]
        return index

    '''
    "Nearest" functions.
    Gets the nearest mesh grid point, or the nearest boundary grid point, if it exists,
    or the nearest point on *any* 1d mesh given a mesh spacing.
    '''
    # TODO: refactor other nearest functions in terms of this one
    def nearest1d(self, x, h):

        basex = h * np.floor( x / h )
        cx = h * np.round( (x / h) % 1 ) + basex
        return cx

    def nearest_grid_point(self, pt):
        x, y = pt
        h = self.dh

        basex, basey = h * np.floor( x / h ), h * np.floor( y / h )
        cx = h * np.round( (x / h) % 1 ) + basex
        cy = h * np.round( (y / h) % 1 ) + basey
        return cx, cy

    def nearest_boundary_point(self, pt):
        h = self.dh 
        dboundary = self.dboundary


        x, y = pt
        cx, cy = self.nearest_grid_point(pt)

        deltax = x - cx # flipped this sign
        deltay = y - cy 
        
        nbrs:[[int], [int]] = [None, None] # every coarse grid point has 2 neighbours
            
        # Get the two adjacent coarse-grid mesh points to this point
        if deltax <= 0:
            nbrs[0] = [cx - h, cy]
        else:
            nbrs[0] = [cx + h, cy]
        if deltay <= 0:
            nbrs[1] = [cx, cy - h]
        else:
            nbrs[1] = [cx, cy + h]
        
        # closest boundary point
        bdry_pt = None
        delta = np.inf

        # TODO: refactor this.
        if nbrs[0][0] in self.bdry_X and nbrs[0][1] in self.bdry_Y:
            bdry_pt = [self.nearest1d(x, dboundary) , cy ]
            delta = np.linalg.norm( bdry_pt - pt)
        
        if nbrs[1][0] in self.bdry_X and nbrs[1][1] in self.bdry_Y:
            temp = [cx, self.nearest1d(y, dboundary)]
            if bdry_pt is None or np.linalg.norm( temp - pt ) < delta:
                bdry_pt = temp
        
        if bdry_pt is not None:
            bdry_pt = np.round(bdry_pt,
                int(np.ceil(np.abs(np.log10(dboundary)))))
        return bdry_pt


    '''
    Hit functions.
    Check if we've moved outside the domain, 
    or if we've hit the cartesian boundary.
    '''

    def is_boundary_point(self, pt):
        X, Y = self.bdry_X, self.bdry_Y
        
        temp_zeros = np.zeros(len(X))
        # The same check as get_index.
        index = np.intersect1d( np.where(np.isclose(X - pt[0], temp_zeros)) , 
                    np.where(np.isclose( Y - pt[1], temp_zeros)) 
                    )

        assert len(index) == 1 or len(index) == 0
        return True if len(index) == 1 else False

    def is_exterior_point(self, pt):
        return True if self.is_outside(pt) > 0 else False

    def hit_boundary(self, pt):
        dboundary, h = self.dboundary, self.dh
        X, Y = self.bdry_X, self.bdry_Y
        x, y = pt
        

        # If the point is in a 'boundary block', it hasn't hit the boundary.        
        nbrs = self.get_grid_block_vertices(pt)        

        if any( [ self.is_exterior_point(nbr) for nbr in nbrs ]): # the current point lives inside a boundary square #TODO: fix this.
            return None

        # is the nearest grid point also a boundary point?
        grid_point = self.nearest_grid_point(pt)
        if grid_point[0] in self.bdry_X and grid_point[1] in self.bdry_Y:
            bdry_pt = self.nearest_boundary_point(pt)

            # is the boundary point close enough to the original point?
            if bdry_pt is not None:
                dist = np.linalg.norm( bdry_pt - np.array(pt) )
                if not dist < 2. * dboundary: # 95 % rule for stdev
                    bdry_pt = None
        return bdry_pt

    def is_outside(self, pt):
        '''
        Return: negative float if point is inside box domain.
                positive float if point is outside box domain.
                zero if point is exactly on box domain 

                # TODO: The zero case is actually mishandled.
                Right now remove the zero case; come back to it later.
        '''
        x, y = pt
        xmin, xmax = self.x_limits
        ymin, ymax = self.y_limits
        return_value = 1. 
        if x > xmin and x < xmax: # todo: flip this cond (better perf?)
            if y > ymin and y < ymax: 
                return_value = -1. 
        
        if x == xmin or x == xmax and y == ymin or y == ymax:
            return_value = 0.
        return return_value

    '''
    Function which makes the only instance of this class.
    '''
    @classmethod
    def generate_test_case(cls, dh, dboundary):
        '''
        Inputs: dx, the size of the spatial mesh.
                dboundary, the size of the boundary mesh.
                Boundary should be less than the spatial mesh.
        '''

        interior_x = 0.00 # only 3 interior points
        interior_y = [0.99, 1.00, 1.01]

        boundary_x = -0.01, 0.01
        boundary_y = (0.98 + np.linspace(0.0, 0.01, 11)[:-1])
        boundary_y = np.append( boundary_y, (0.99 + np.linspace(0.0, 0.01, 11)[:-1]) )
        boundary_y = np.append( boundary_y, (1. + np.linspace(0.0, 0.01, 11)[:-1]) )
        boundary_y = np.append( boundary_y, (1.01 + np.linspace(0.0, 0.01, 11)) )

        X, Y = np.meshgrid( boundary_x , boundary_y )
        X, Y = X.reshape(-1,), Y.reshape(-1,)

        boundary_y = 0.98, 1.02
        boundary_x = (-0.01 + np.linspace(0.0, 0.01, 11)[1:-1]) # truncate leftmost point so we don't duplicate
        boundary_x = np.append(boundary_x, 0 + np.linspace(0.0, 0.01, 11)[:-1])

        Xt, Yt = np.meshgrid( boundary_x, boundary_y )
        Xt, Yt = Xt.reshape(-1,), Yt.reshape(-1,)

        X, Y = np.append(X, Xt), np.append(Y, Yt)


        #Convoluted Rounding
        X = np.round(X ,int(np.ceil(np.abs(np.log10(dboundary)))))
        Y = np.round(Y ,int(np.ceil(np.abs(np.log10(dboundary)))))

        # convoluted check for duplicate points
        assert( len(np.unique( np.array( 
                            list(map( np.array, zip(X, Y) ))), axis=0)) == len(X))
        assert np.allclose( [dboundary], [boundary_x[-1] - boundary_x[-2]])

        return X, Y # The set of boundary points.


class MonteCarloSimulator():
    
    def __init__(self, g=None):
        self.domain = BoxWithCartesianBoundary()

        N = len(self.domain.bdry_X) # number of boundary points
        self.frequencies = np.zeros( (N, N) )
        self.total_boundary_values = np.zeros(N)

        # internal info counters
        self.ofreq, self.bfreq = 0, 0

        self.nsamples = 100
        self.maxsteps = 100 # tweakable

        # so 95% of points take a step of size 0.5 dboundary or less
        self.dt = 0.5 * (self.domain.dboundary ** 2) 

        if g: 
            self.g = g 
        else:
            self.g = lambda x, y: 0. # laplacian

    def simulate_point(self, pt):
        np.random.seed(2)

        dt = self.dt
        
        pt_index = self.domain.get_index(pt)

        for i in range(self.nsamples): #self.nsamples
            x, y = pt # reset the point location
            
            for j in range(self.maxsteps): # nsteps of random walk
                x, y = x + np.random.randn() * np.sqrt(dt), y + np.random.randn() * np.sqrt(dt)
                # Did I hit the domain exterior?

                # careful here.
                if self.domain.is_outside( (x, y) ) >= 0:  # Cost: one function evaluation                                        
                    self.ofreq += 1

                    #print(f"{j}, g(j) is {self.g(x,y)}")
                    self.total_boundary_values[pt_index] += self.g(x, y)
                    break

                # Did I hit the grid boundary?
                bdry_pt = self.domain.hit_boundary( (x, y) )
                if bdry_pt is not None:
                    self.bfreq += 1
                    index = self.domain.get_index(bdry_pt)
                    self.frequencies[pt_index, index] += 1
                    break

                continue

        print(f" Coupled: {self.bfreq}, Outside: {self.ofreq}")

    
    
    def simulate(self):
        #points = zip(self.domain.bdry_X, self.domain.bdry_Y)

        points = [[0.01, 1.]]
        for pt in points:
            self.simulate_point(pt)

    def solve_coupling(self):
        N = self.nsamples 
        X, Y = self.domain.bdry_X, self.domain.bdry_Y

        # Maybe divide by N later if the boundary values overflow
        A = ( N * np.eye( len(X), len(X) ) - self.frequencies )
        b = -self.total_boundary_values

        print(b)

        assert np.linalg.matrix_rank(A) == len(X)

        u = np.linalg.solve(A , b)
        gg = np.vectorize(self.g)

        soln = gg(X, Y)
        
        print( (u - soln)) # My god this convergence is slow as shit
        print(soln)
        print(u)

        # The L2 norm error:
        # error = (dboundary * (len(X) - 1) * np.sum(np.abs(u - soln))) - 0.5 * (np.abs((u -soln))[0] + np.abs((u - soln))[-1])

if __name__ == '__main__':
    
    print("nothing here")

    simulator = MonteCarloSimulator(lambda x, y: 
        x ** 3 + y ** 3 - 3 * x ** 2 * y -3 * y ** 2 * x + 5.
    )
    simulator.simulate()
    simulator.solve_coupling()
    


