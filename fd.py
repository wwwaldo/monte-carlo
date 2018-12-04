import numpy as np
import scipy as sp
import scipy.sparse as spp
import scipy.sparse.linalg as sppl

def grid_gen(boundary_func, min_xs, max_xs, Ns):
    grid_to_vec = {}
    vec_to_grid = []
    boundary = {}
    k = 0
    hs = (max_xs - min_xs)/(Ns-1)
    grid_to_point = lambda p: np.array([p[0]*hs[0], p[1]*hs[1]]) + min_xs
    for i in range(Ns[0]):
        for j in range(Ns[1]):
            p = np.array([i*hs[0], j*hs[1]]) + min_xs
            p = grid_to_point((i,j))
            if boundary_func(p) < 0:
                b1 = boundary_func(grid_to_point((i,j+1)))
                b2 = boundary_func(grid_to_point((i,j-1)))
                b3 = boundary_func(grid_to_point((i+1,j)))
                b4 = boundary_func(grid_to_point((i-1,j)))
                if b1 < 0 and b2 < 0 and b3 < 0 and b4 < 0:
                    grid_to_vec[(i,j)] = k
                    vec_to_grid.append((i,j))
                    k += 1
                else:
                    boundary[(i,j)] = None
    return grid_to_vec, vec_to_grid, boundary, grid_to_point


def gen_L(find_val, poisson_f, grid_to_vec, vec_to_grid, boundary, grid_to_point, h):
    # Generates Lu = b
    # This tries to solve -Delta*u = f, with u=g on the boundary.
    # f is given by poisson_f, and g is found via find_val.
    N = len(vec_to_grid)
    L = -4*spp.identity(N).todok()
    b = (spp.eye(N,1) - spp.eye(N,1)).todok()
    for i in range(N):
        p = vec_to_grid[i]
        b[i] = -poisson_f(grid_to_point(p))*h*h
        for (di,dj) in [(1,0), (-1,0), (0,1), (0,-1)]:
            new_p = (p[0] + di, p[1] + dj)
            if new_p in grid_to_vec:
                j = grid_to_vec[new_p]
                L[i,j] = 1
            else:
                if boundary[new_p] is None:
                    boundary[new_p] = find_val(new_p)
                b[i] += -boundary[new_p]
    return L.tocsr(), b.tocsr()

if __name__ == "__main__":
    boundary_func = lambda x: x[0]*x[0] + x[1]*x[1] - 1
    N = 100
    grid_to_vec, vec_to_grid, boundary, grid_to_point = grid_gen(boundary_func, np.array([-1,-1]), np.array([1,1]), np.array([N+1,N+1]))

    f = lambda x: (lambda x: x[0]**3 + x[1]**3 - 3*x[0]*x[0]*x[1] - 3*x[0]*x[1]*x[1])(grid_to_point(x))
    g = lambda x: 0
    # f = lambda x: x[0]**2 + x[1]**2
    # g = lambda x: -4
    L, b = gen_L(f, g, grid_to_vec, vec_to_grid, boundary, grid_to_point, 2/N)

    x = sppl.spsolve(L, b)
    x_true = np.vectorize(lambda x: f(grid_to_point(x)), signature='(n)->()')(np.array(vec_to_grid))
    print(np.max(np.abs(x-x_true)))
