import numpy as np
import matplotlib.pyplot as plt
import traceback
import itertools
from fd import *

class SetWithCartesianBoundary():
    '''
    Params for this test box.
    '''
    def __init__(self, boundary_dict, grid_to_point, is_outside_fn, h):
        self.dh = h
        self.dboundary = h

        bdry = []
        self.bdry_idx = []
        self.bdry_dict = {}
        i = 0
        for p in boundary_dict:
            bdry.append(grid_to_point(p))
            self.bdry_idx.append(p)
            self.bdry_dict[p] = i
            i += 1
        bdry = np.array(bdry)
        self.bdry_X , self.bdry_Y = bdry[:,0], bdry[:,1]

        self.is_outside = is_outside_fn

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
        
        # Throws an obscure error if a non-bdry point is passed in
        #try:
        index = np.intersect1d( np.where(np.isclose(X, pt[0])) , 
                        np.where(np.isclose(Y, pt[1])) 
                        )[0]
        #except Exception as e:
        #    print(f'Tried to index {pt} as boundary point:')
        #    raise e

        return index

    '''
    "Nearest" functions.
    Gets the nearest mesh grid point, or the nearest boundary grid point, if it exists,
    or the nearest point on *any* 1d mesh given a mesh spacing.
    '''
    # TODO: refactor other nearest functions in terms of this one
    def nearest1d(self, x, h):
        return h * np.round(x / h)

    def nearest_grid_point(self, pt):
        x, y = pt
        h = self.dh

        cx = h * np.round(x / h)
        cy = h * np.round(y / h)
        return cx, cy

    def nearest_boundary_point(self, pt):
        '''
        Bug: The point returned may not be a real boundary point.
        Should be fixed now.
        '''
        h = self.dh 
        dboundary = self.dboundary


        x, y = pt
        # is the nearest grid point also a boundary point?
        cx, cy = self.nearest_grid_point(pt)
        if not self.is_boundary_point((cx, cy)):
            return None

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
        if self.is_boundary_point(nbrs[0]):
            bdry_pt = [self.nearest1d(x, dboundary) , cy ]
            delta = np.linalg.norm( np.array(bdry_pt) - np.array(pt))
        
        if self.is_boundary_point(nbrs[1]):
            temp = [cx, self.nearest1d(y, dboundary)]
            if bdry_pt is None or np.linalg.norm( np.array(temp) - np.array(pt) ) < delta:
                bdry_pt = temp
        
        if bdry_pt is not None:
            bdry_pt = dboundary * np.round(np.array(bdry_pt) / dboundary)
        return bdry_pt


    '''
    Hit functions.
    Check if we've moved outside the domain, 
    or if we've hit the cartesian boundary.
    '''

    def is_boundary_point(self, pt):
        X, Y = self.bdry_X, self.bdry_Y
        
        # The same check as get_index.
        index = np.intersect1d( np.where(np.isclose(X, pt[0])) , 
                    np.where(np.isclose(Y, pt[1])) 
                    )

        assert len(index) == 1 or len(index) == 0
        return (len(index) == 1)

    def is_exterior_point(self, pt):
        return (self.is_outside(pt) > 0)

    def hit_boundary(self, pt):
        dboundary, h = self.dboundary, self.dh
        X, Y = self.bdry_X, self.bdry_Y
        x, y = pt
        

        # If the point is in a 'boundary block', it hasn't hit the boundary.        
        nbrs = self.get_grid_block_vertices(pt)        

        if any( [ self.is_exterior_point(nbr) for nbr in nbrs ]): # the current point lives inside a boundary square #TODO: fix this.
            return None

        bdry_pt = self.nearest_boundary_point(pt)

        # is the boundary point close enough to the original point?
        if bdry_pt is not None:
            dist = np.linalg.norm( bdry_pt - np.array(pt) )
            if not dist < 2. * dboundary: # 95 % rule for stdev
                bdry_pt = None
        return bdry_pt

    def visualize_boundary(self):
        fig, ax = plt.subplots()
        ax.scatter(self.bdry_X, self.bdry_Y)

        X, Y = np.linspace(-1.5, 1.5), np.linspace(-1.5, 1.5) # hard-coded
        X, Y = np.meshgrid(X, Y)
        outside_v = squiggly_domain(0.45, True)

        ax.contour(X, Y, outside_v(X, Y), levels=[0.0])
        plt.show()

    @classmethod
    def generate_test_case_2(cls):
        dboundary = 0.0001
        '''
        Inputs: dx, the size of the spatial mesh.
                dboundary, the size of the boundary mesh.
                Boundary should be less than the spatial mesh.
        '''

        interior_x = 0.00 # only 3 interior points
        interior_y = [0.99, 1.00, 1.01]

        boundary_x = -0.01, 0.01
        boundary_y = (0.98 + np.linspace(0.0, 0.01, 11)[:-1])
        boundary_y = np.append( boundary_y, (0.99 + np.linspace(0.0, 0.01, 101)[:-1]) )
        boundary_y = np.append( boundary_y, (1. + np.linspace(0.0, 0.01, 101)[:-1]) )
        boundary_y = np.append( boundary_y, (1.01 + np.linspace(0.0, 0.01, 101)) )

        X, Y = np.meshgrid( boundary_x , boundary_y )
        X, Y = X.reshape(-1,), Y.reshape(-1,)

        boundary_y = 0.98, 1.02
        boundary_x = (-0.01 + np.linspace(0.0, 0.01, 11)[1:-1]) # truncate leftmost point so we don't duplicate
        boundary_x = np.append(boundary_x, 0 + np.linspace(0.0, 0.01, 101)[:-1])

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
    
    def __init__(self, domain, nsamples, g=None):
        self.domain = domain

        N = len(self.domain.bdry_X) # number of boundary points
        self.frequencies = np.zeros( (N, N) )
        self.total_boundary_values = np.zeros(N)

        # internal info counters
        self.ofreq, self.bfreq = 0, 0

        self.nsamples = nsamples
        self.maxsteps = 100 # tweakable

        # so 95% of points take a step of size 0.5 dboundary or less
        self.dt = 0.25 * (self.domain.dboundary ** 2) 

        if g: 
            self.g = g 
        else:
            self.g = lambda x, y: 0. # laplacian


    def simulate_point(self, pt):
        dt = self.dt
        
        pt_index = self.domain.get_index(pt)

        for i in range(self.nsamples): #self.nsamples
            x, y = pt # reset the point location
            
            for j in itertools.count():
            #for j in range(self.maxsteps): # nsteps of random walk
                x, y = x + np.random.randn() * np.sqrt(dt), y + np.random.randn() * np.sqrt(dt)
                # print(x, y)
                # Did I hit the domain exterior?

                # careful here.
                if self.domain.is_outside( (x, y) ) >= 0:  # Cost: one function evaluation                                        
                    self.ofreq += 1

                    #print(f"{j}, g(j) is {self.g(x,y)}")
                    self.total_boundary_values[pt_index] += self.g(x, y)
                    break

                """# Did I hit the grid boundary?
                bdry_pt = self.domain.hit_boundary( (x, y) )
                if bdry_pt is not None: # and not np.allclose(bdry_pt,pt):
                    self.bfreq += 1
                    index = self.domain.get_index(bdry_pt)
                    self.frequencies[pt_index, index] += 1
                    break"""
                
                continue

        print(f" Coupled: {self.bfreq}, Outside: {self.ofreq}")

    
    
    def simulate(self):
        points = zip(self.domain.bdry_X, self.domain.bdry_Y)
        #points = list(points)[0:1]

        #points = [[0.01, 1.]]
        for pt in points:
            self.simulate_point(pt)

    def solve_coupling(self):
        N = self.nsamples 
        X, Y = self.domain.bdry_X, self.domain.bdry_Y

        # Maybe divide by N later if the boundary values overflow
        A = ( np.eye( len(X), len(X) ) - self.frequencies / N )
        b = self.total_boundary_values / N
        # A = ( np.eye( len(X), len(X) ))
        # b = self.total_boundary_values / N + self.frequencies.dot(self.total_boundary_values) / N / N

        print(b)

        assert np.linalg.matrix_rank(A) == len(X)

        u = np.linalg.solve(A , b)
        gg = np.vectorize(self.g)

        soln = gg(X, Y)

        print(self.frequencies)
        
        print( (u - soln))
        print(soln)
        print(u)

        # The L2 norm error:
        error = ((self.domain.dboundary * np.sum(np.abs(u - soln))) - 0.5 * self.domain.dboundary * (np.abs((u -soln))[0] + np.abs((u - soln))[-1]))/len(X)
        print(error)
        rerror = ((self.domain.dboundary * np.sum(np.abs(u - soln)/soln)) - 0.5 * self.domain.dboundary * (np.abs((u -soln)/soln)[0] + np.abs((u - soln)/soln)[-1]))/len(X)
        print(rerror)
        merror = np.max(np.abs(u-soln)/soln)
        print(merror)
        return u


def squiggly_domain(spikiness, meshgrid=False):
    if meshgrid: # for compatibility
        def squiggly(X, Y):
            Z = (X ** 2 + Y ** 2 ) 
            theta = np.arctan2( X, Y )
            for i in range(3, 7):
                Z -= (1. / ((X ** 2 + Y ** 2) * (0.5 * i)) ) * np.abs(np.sin(spikiness * np.pi * i * theta))
            return Z
    else:
        def squiggly(pt):
            X, Y = pt
            Z = (X ** 2 + Y ** 2 ) 

            if np.allclose(pt, 0):
                return -100

            theta = np.arctan2( X, Y )
            

            for i in range(3, 7):
                Z -= (1. / ((X ** 2 + Y ** 2) * (0.5 * i)) ) * np.abs(np.sin(spikiness * np.pi * i * theta))
            return Z
    return squiggly


if __name__ == '__main__':
    import shelve 
    data = None
    
    with shelve.open('data.dat') as shelf:
        try:
            data = shelf['squiggly_N100_nsamples_400']
        except Exception as e:
            print("Need to regenerate data")
    
    print("nothing here")

    boundary_func = lambda x: x[0]*x[0] + x[1]*x[1] - 1
    boundary_func = squiggly_domain(0.45)

    N = 150
    nsamples = 400
    grid_to_vec, vec_to_grid, boundary, grid_to_point = grid_gen(boundary_func, np.array([-1.5,-1.5]), np.array([1.5,1.5]), np.array([N+1,N+1]))
    h = 2/N

    domain = SetWithCartesianBoundary(boundary, grid_to_point, boundary_func, h)

    domain.visualize_boundary()
    


    f = lambda x: x[0]**3 + x[1]**3 - 3*x[0]*x[0]*x[1] - 3*x[0]*x[1]*x[1] + 1
    g = lambda x: 0
    # f = lambda x: x[0]**2 + x[1]**2
    # g = lambda x: -4

    simulator = MonteCarloSimulator(domain, nsamples,
        lambda x, y: x ** 3 + y ** 3 - 3 * x ** 2 * y -3 * y ** 2 * x + 1
    )

    if data is None:
        simulator.simulate()
        u = simulator.solve_coupling()

        # shelve the data 
        import shelve 
        with shelve.open('data.dat') as shelf:
            shelf['squiggly_N100_nsamples_400'] = u
        data = u

    u = data

    def find_val(p):
        return u[domain.bdry_dict[p]]

    L, b = gen_L(find_val, g, grid_to_vec, vec_to_grid, boundary, grid_to_point, h)

    x = sppl.spsolve(L, b)
    x_true = np.vectorize(lambda x: f(grid_to_point(x)), signature='(n)->()')(np.array(vec_to_grid))
    print(np.max(np.abs(x-x_true)))
    print(np.max(np.abs(x-x_true)/x_true))


    find_val_exact = lambda x: f(grid_to_point(x))
    for k in boundary:
        boundary[k] = None

    L, b = gen_L(find_val_exact, g, grid_to_vec, vec_to_grid, boundary, grid_to_point, h)

    x = sppl.spsolve(L, b)
    x_true = np.vectorize(lambda x: f(grid_to_point(x)), signature='(n)->()')(np.array(vec_to_grid))
    print(np.max(np.abs(x-x_true)))
    print(np.max(np.abs(x-x_true)/x_true))
