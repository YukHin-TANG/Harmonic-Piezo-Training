"""
This file contains the piezo-solid truss class, which is used to calculate the statics of a piezo-solid truss structure
"""
import numpy as np
import time
import tqdm as tqdm


class TransientPiezoTruss:
    """
    static_piezo_truss class
    ========================
        This class is used to calculate the temporal response of a piezo-solid truss structure
    """
    edges, pts = [], []
    node_state = []
    Kij, Cij = [], []
    extForce, extDisp, extVel, extAccel = [], [], [], []
    m, C, G, freq = 0, 0, 0, None
    T, calDisp, calQ = None, None, None
    Vt = None

    def __init__(self, size=(3, 2), m=1E-2, eta=0.04, kij=2.2E4, G=2.72E-3, C=1.44E-8, freq=[0], load=False, path=None):
        """
        transient_piezo_truss constructor
        ===================
        :param
        """
        #
        self.G, self.C = G, C
        self.freq = np.array(freq)
        self.edges = []
        self.pts = []
        if load:
            lat = np.load(path)
            self.pts = lat['nodes']
            self.edges = lat['edges']
        else:
            # size of the grid
            M, N = size
            X, Y = np.meshgrid(np.arange(M), np.arange(N))
            # connectivity from the grid of points
            triangles = np.array([X.flatten(), Y.flatten()]).T
            #
            self.edges = []
            self.pts = []
            for i, j in triangles:
                self.pts += [[i + j * np.sin(np.pi / 6), j * np.cos(np.pi / 6)]]
                # Equilateral triangle edges indices
                if i < M - 1:
                    self.edges += [[i + j * M, i + j * M + 1]]
                if j < N - 1:
                    self.edges += [[i + j * M, i + (j + 1) * M]]
                if i > 0 and j < N - 1:
                    self.edges += [[i + j * M, i + (j + 1) * M - 1]]
            #
            self.pts = np.array(self.pts)
            self.edges = np.array(self.edges)
            #
            self.triangle = triangles
        #
        cij = 2 * eta * np.sqrt(kij * m)
        #
        # length of the spring, cos(theta), sin(theta)
        self.Kij = np.zeros((len(self.pts) * 2, len(self.pts) * 2))
        self.Cij = np.zeros((len(self.pts) * 2, len(self.pts) * 2))
        self.m = m
        # Mij = np.zeros((len(pts)*2, len(pts)*2))
        for p, q in self.edges:
            _l = np.linalg.norm(self.pts[p] - self.pts[q])
            c = (self.pts[p, 0] - self.pts[q, 0]) / _l
            s = (self.pts[p, 1] - self.pts[q, 1]) / _l
            # print(f"l_{n} = {l:.2f}, c_{n} = {c:.2f}, s_{n} = {s:.2f}")
            A = np.array([[c ** 2, c * s],
                          [c * s, s ** 2]])
            #
            self.Kij[2 * p:(2 * p + 2), 2 * p:(2 * p + 2)] += +A * kij
            self.Kij[2 * q:(2 * q + 2), 2 * q:(2 * q + 2)] += +A * kij
            self.Kij[2 * p:(2 * p + 2), 2 * q:(2 * q + 2)] += -A * kij
            self.Kij[2 * q:(2 * q + 2), 2 * p:(2 * p + 2)] += -A * kij
            #
            self.Cij[2 * p:(2 * p + 2), 2 * p:(2 * p + 2)] += +A * cij
            self.Cij[2 * q:(2 * q + 2), 2 * q:(2 * q + 2)] += +A * cij
            self.Cij[2 * p:(2 * p + 2), 2 * q:(2 * q + 2)] += -A * cij
            self.Cij[2 * q:(2 * q + 2), 2 * p:(2 * p + 2)] += -A * cij
        #
        self.extForce = np.zeros((len(self.pts), 2))
        self.extDisp = np.zeros((len(self.pts), 2))
        self.extVel = np.zeros((len(self.pts), 2))
        self.extAccel = np.zeros((len(self.pts), 2))

    def extU(self, _t, nodes, u, reset=True):
        self.node_state = np.array([True] * len(self.pts) * 2)
        self.node_state[nodes] = False
        #
        if reset:
            self.extDisp = np.zeros((len(self.pts), 2))
            self.extVel = np.zeros((len(self.pts), 2))
            self.extAccel = np.zeros((len(self.pts), 2))
        #
        __u = np.array(u)
        #
        __Ut = np.zeros((len(self.pts) * 2,))
        __Ut[nodes] = __u[:, :, 0] @ np.sin(2 * np.pi * self.freq * _t)
        __Ut[nodes] += __u[:, :, 1] @ np.cos(2 * np.pi * self.freq * _t)
        #
        __Vt = np.zeros((len(self.pts) * 2,))
        __Vt[nodes] = __u[:, :, 0] * 2 * np.pi * self.freq @ np.cos(2 * np.pi * self.freq * _t)
        __Vt[nodes] += -__u[:, :, 1] * 2 * np.pi * self.freq @ np.sin(2 * np.pi * self.freq * _t)
        #
        __At = np.zeros((len(self.pts) * 2,))
        __At[nodes] = -__u[:, :, 0] * (2 * np.pi * self.freq) ** 2 @ np.sin(2 * np.pi * self.freq * _t)
        __At[nodes] += -__u[:, :, 1] * (2 * np.pi * self.freq) ** 2 @ np.cos(2 * np.pi * self.freq * _t)
        #
        self.extDisp += __Ut.reshape((-1, 2))
        self.extVel += __Vt.reshape((-1, 2))
        self.extAccel += __At.reshape((-1, 2))

    def extF(self, _t, nodes, f, reset=True):
        #
        if reset:
            self.extForce = np.zeros((len(self.pts), 2))
        #
        __f = np.array(f)
        #
        __Ft = np.zeros((len(self.pts) * 2,))
        __Ft[nodes] = __f[:, :, 0] @ np.sin(2 * np.pi * self.freq * _t)
        __Ft[nodes] += __f[:, :, 1] @ np.cos(2 * np.pi * self.freq * _t)
        #
        self.extForce += __Ft.reshape((-1, 2))

    def extVolForce(self, _t, V, reset=True):
        # Calculate the direction of the edges
        vec = self.pts[self.edges[:, 1], :] - self.pts[self.edges[:, 0], :]
        norm_vec = np.linalg.norm(vec, axis=1)[:, np.newaxis]
        vec /= norm_vec
        # Create an empty array to store the voltage force
        voltage_force = np.empty((self.edges.shape[0], 2))
        # Assign the voltage force to the array
        #
        __V = np.array(V)
        Vt = __V[:, :, 0] @ np.sin(2 * np.pi * self.freq * _t)
        Vt += __V[:, :, 1] @ np.cos(2 * np.pi * self.freq * _t)
        self.Vt = Vt
        #
        voltage_force[:, 0] = vec[:, 0] * (Vt * self.G)
        voltage_force[:, 1] = vec[:, 1] * (Vt * self.G)
        # Add the voltage force to the external force
        if reset:
            self.extForce = np.zeros((len(self.pts), 2))
        for n, e in enumerate(self.edges):
            self.extForce[e[0], :] -= voltage_force[n, :]
            self.extForce[e[1], :] += voltage_force[n, :]

    def derivation(self, t, y, UNodes, u0, FNodes=None, f=None, V=None):
        #
        self.extU(t, UNodes, u0)
        #
        if FNodes is not None:
            self.extF(t, FNodes, f)
        #
        if V is not None:
            self.extVolForce(t, V)
        #
        _Ft = self.extForce.flatten()
        _Ft -= np.dot(self.Kij, self.extDisp.flatten()) + \
               np.dot(self.Cij, self.extVel.flatten()) + \
               self.extAccel.flatten() * self.m
        #
        Ftff = _Ft.flatten()[self.node_state]
        Kff = self.Kij[self.node_state][:, self.node_state]
        Cff = self.Cij[self.node_state][:, self.node_state]
        #
        u = y.reshape(2, -1)[0]
        du = y.reshape(2, -1)[1]
        #
        ddu = (Ftff - np.dot(Cff, du) - np.dot(Kff, u)) / self.m

        return np.concatenate([du, ddu])

    def solve(self, tspan, UNodes, u, FNodes=None, f=None, V=None):
        """
        solve method
        ============
        :param tspan:
        :param UNodes:
        :param u:
        :param FNodes:
        :param f:
        :param V:
        :return:
        """
        from scipy.integrate import solve_ivp
        disp = np.zeros((len(self.pts) * 2, len(tspan)))
        if len(self.pts) * 2 != len(UNodes):
            y0 = np.zeros(((len(self.pts) * 2 - len(UNodes)) * 2))
            sol = solve_ivp(self.derivation, [tspan[0], tspan[-1]], y0, method='RK45', t_eval=tspan,
                            args=(UNodes, u, FNodes, f, V))
            self.T = sol.t
            #
            disp[self.node_state] = sol.y[:len(sol.y) // 2]
        #
        V_series = np.zeros((len(self.edges), len(tspan)))
        for n in range(len(tspan)):
            self.extU(tspan[n], UNodes, u)
            self.extVolForce(tspan[n], V)
            disp[~self.node_state, n] = self.extDisp.flatten()[~self.node_state]
            V_series[:, n] = self.Vt
        self.calDisp = disp.reshape((-1, 2, len(tspan)))
        #
        self.calQ = -self.calc_bond_stretch() * self.G + V_series * self.C

    def calc_bond_stretch(self):
        """
        Calculate the stress of each element.
        :return: bond stretch, bond inner force
        """
        __idx = self.edges[:, 0]
        __idy = self.edges[:, 1]
        __disps_diff = self.calDisp[__idy, :] - self.calDisp[__idx, :]  # displacement difference
        # print("_______________")
        # print(__disps_diff)
        #
        __bond_vector = self.pts[__idy, :] - self.pts[__idx, :]  # bond vector
        __bond_length = np.linalg.norm(__bond_vector, axis=1)  # bond length
        __bond_direction = __bond_vector / __bond_length[:, np.newaxis]  # bond direction
        # bond stretch  # bond inner force
        __bond_stretch = np.sum(__disps_diff * __bond_direction[:, :, np.newaxis], axis=1)
        #
        return __bond_stretch

    def plotStructure(self):
        """
        plotStructure method
        ====================
        :return:
        """
        import matplotlib.pyplot as plt
        plt.axis('equal')
        plt.plot(np.array(self.pts)[:, 0], np.array(self.pts)[:, 1], 'o', markersize=10)
        for i, j in self.edges:
            plt.plot([self.pts[i][0] , self.pts[j][0]], [self.pts[i][1], self.pts[j][1]], 'k:')
        # label the nodes
        for n, (x, y) in enumerate(self.pts):
            plt.text(x, y, str(n), fontsize=20)
        plt.show()


def sinBoundary(A, B, phi=0, U0=1E-6):
    nodes = []
    u = []
    for a in range(A):
        nodes += [2 * a + 1]
        u += [[-np.sin(np.pi * a / (A - 1)) * np.sin(phi), -np.sin(np.pi * a / (A - 1)) * np.cos(phi)]]
        #
        nodes += [2 * a + 1 + 2 * (B - 1) * A]
        u += [[np.sin(2 * np.pi * a / (A - 1)) * np.sin(phi), np.sin(2 * np.pi * a / (A - 1)) * np.cos(phi)]]
    ##
    for b in range(B):
        nodes += [2 * A - 2 + 2 * b * A]
        u += [[np.sin(0.5 * np.pi * b / (B - 1)) * np.sin(phi), np.sin(0.5 * np.pi * b / (B - 1)) * np.cos(phi)]]
        nodes += [2 * b * A]
        u += [[np.sin(0.5 * np.pi * b / (B - 1)) * np.sin(phi), np.sin(0.5 * np.pi * b / (B - 1)) * np.cos(phi)]]
    return np.array(nodes), np.array(u) * U0


def test_transient_piezo_truss():
    """
    Test the static_piezo_truss class
    ===================
    :return:
    """
    # Create a static_piezo_truss object
    M, N = 3, 2
    freq = [5]
    transient_piezo_truss = TransientPiezoTruss(size=(M, N), freq=freq)
    # Plot the structure
    UNodes_free = [0, 1, 2 * M - 1]
    u_free = [[[0, 0]], [[0, 0]], [[0, 0]]]
    #
    UNodes_clamped = [0, 1, 4, 2 * M - 1, ]
    u_clamped = [[[0, 0]], [[0, 0]], [[1E-7, 0]], [[0, 0]]]
    #
    T = np.linspace(0, 1, 1000)
    N = len(T) // 2
    Vol = np.ones((len(transient_piezo_truss.edges), len(freq), 2)) * 10
    Vol[:, :, 1] = 0
    #
    Epoch = 25
    lr = 1E9
    # Vol = np.ones((len(transient_piezo_truss.edges), len(freq), 2)) * 10
    for i in tqdm.tqdm(range(Epoch)):
        #
        transient_piezo_truss.solve(T, UNodes_free, u_free, V=Vol)
        Q_free = transient_piezo_truss.calQ.copy()
        #
        transient_piezo_truss.solve(T, UNodes_clamped, u_clamped, V=Vol)
        Q_clamped = transient_piezo_truss.calQ.copy()
        #
        for n, fq in enumerate(freq):
            Q_free_sin = 2 * Q_free[:, N:] @ np.sin(fq * 2 * np.pi * T[N:]) / len(T[N:])
            Q_free_cos = 2 * Q_free[:, N:] @ np.cos(fq * 2 * np.pi * T[N:]) / len(T[N:])
            #
            Q_clamped_sin = 2 * Q_clamped[:, N:] @ np.sin(fq * 2 * np.pi * T[N:]) / len(T[N:])
            Q_clamped_cos = 2 * Q_clamped[:, N:] @ np.cos(fq * 2 * np.pi * T[N:]) / len(T[N:])
            #
            Vol[:, n, 0] += lr * (Q_clamped_sin - Q_free_sin)
            Vol[:, n, 1] += lr * (Q_clamped_cos - Q_free_cos)
    #
    print(Vol)
    # Plot the displacement
    # import matplotlib.pyplot as plt
    # plt.plot(T, transient_piezo_truss.calc_bond_stretch()[0], label='Q')
    # plt.legend()
    # plt.show()
    # print(f"Time elapsed: {time.time() - timer:.2f} s")


if __name__ == '__main__':
    test_transient_piezo_truss()
