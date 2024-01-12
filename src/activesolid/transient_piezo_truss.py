"""
This file contains the piezo-solid truss class, which is used to calculate the statics of a piezo-solid truss structure
"""
import numpy as np

# Import standard library modules
from truss import *
import matplotlib.pyplot as plt
import scipy.sparse as ss  # sparse matrix
import scipy


class TransientPiezoTruss(Truss):
    """
    static_piezo_truss class
    ========================
        This class is used to calculate the temporal response of a piezo-solid truss structure
    """
    Steps, eta = None, None
    G, C, extVoltage = None, None, None
    M, H = None, None
    extVel, extAccel = None, None
    calVel, calAccel = None, None

    def __init__(self, node_num_4_each_edge: int = 3, path=None, tspan=None, dt: float = 1.0e-3,
                 gamma: float = 1.0e-4, capacitance: float = 1.0e-8, k: float = 1.0e4,
                 xi: float = 1.2, m: float = 1.0E-1, freq: list = 100.0):
        """
        transient_piezo_truss constructor
        ===================
        :param node_num_4_each_edge: int
            Number of nodes for each edge of the lattice.
        """
        super().__init__(node_num_4_each_edge, path, )  # call the constructor of the parent class
        #
        if tspan is None:
            tspan = [0, 1]
        self.k, self.m, self.dt = k, m, dt
        self.omega = np.array(freq) * 2 * np.pi
        #
        self.Steps = np.arange(tspan[0], tspan[1], dt)
        #
        self.eta = 2 * xi * (self.m * self.weights) ** 0.5
        self.G = np.ones((self.elements.shape[0],)) * gamma
        self.C = np.ones((self.elements.shape[0],)) * capacitance
        #
        self.extVoltage = np.vectorize(lambda t: np.zeros((self.nodes.shape[0], 2, 1)))
        self.extForces = np.vectorize(lambda t: np.zeros((self.nodes.shape[0], 2, 1)))
        self.extDisps = np.vectorize(lambda t: np.zeros((self.nodes.shape[0], 2, 1))*np.nan)  # nan means free
        self.extVel = np.vectorize(lambda t: np.zeros((self.nodes.shape[0], 2, 1)))
        self.extAccel = np.vectorize(lambda t: np.zeros((self.nodes.shape[0], 2, 1)))

    def add_volForce(self, voltage: list = None, reset: bool = False, seed: int = 0):
        """
        Add voltage to the piezo-solid truss structure
        ===================
        :param voltage: float
            Voltage to be added to the piezo-solid truss structure
        :param reset: bool
            Whether to reset the external force
        :param seed: int
            Random seed
                0 means identical voltage for each node otherwise uniform distribution.
        :return:
        """
        # print(self.extVoltage.shape)
        if reset:
            self.extForces *= 0
        if voltage is not None:
            __Vi = np.array(voltage)
            if seed != 0:
                np.random.seed(seed)
                __rnd = np.random.uniform(0, 2, (self.elements.shape[0]))
            else:
                __rnd = np.ones((self.elements.shape[0],))
            #
            self.extVoltage = __rnd[:, np.newaxis, np.newaxis] * __Vi
        __Voltage = self.extVoltage[:, 0, :] @ np.sin(self.omega[:, np.newaxis] * self.Steps)
        __Voltage += self.extVoltage[:, 1, :] @ np.cos(self.omega[:, np.newaxis] * self.Steps)
        # Calculate the direction of the edges
        vec = self.nodes[self.elements[:, 1], :] - self.nodes[self.elements[:, 0], :]
        # print(vec.shape)
        norm_vec = np.linalg.norm(vec, axis=1)[:, np.newaxis]
        vec /= norm_vec
        # Create an empty array to store the voltage force
        voltage_force = np.empty((self.elements.shape[0], 2, len(self.Steps)))
        # Assign the voltage force to the array
        voltage_force[:, 0, :] = (vec[:, 0] * self.G)[:, np.newaxis] * __Voltage
        voltage_force[:, 1, :] = (vec[:, 1] * self.G)[:, np.newaxis] * __Voltage
        # Add the voltage force to the external force
        for n, e in enumerate(self.elements):
            self.extForces[e[0], :] -= voltage_force[n, :]
            self.extForces[e[1], :] += voltage_force[n, :]

    def add_extTemporalForces(self, nodes: np.ndarray = None, forces: np.ndarray = None, reset: bool = False):
        """
        Add external forces to the truss.
        Zero forces are added to the nodes that are not specified.
        ==================
        :param nodes: shape = (n, )
        :param forces: shape = (*, 2, 2, m) or (n, 2)
        :param reset: reset external forces to zero
        :return:
        """
        if reset:
            self.extForces = np.vectorize(lambda t: np.zeros((self.nodes.shape[0], 2, 1)))
        if nodes is not None:
            def __extForces(t):
                __calForces = np.zeros((self.nodes.shape[0], 2, 1))
                __calForces[nodes] = forces[:, :, 0, :] @ np.sin(self.omega[:, np.newaxis] * t)
                __calForces[nodes] += forces[:, :, 1, :] @ np.cos(self.omega[:, np.newaxis] * t)
                return __calForces

            self.extForces = np.vectorize(__extForces)

    def add_extTemporalDisps(self, nodes: np.ndarray = None, disps: np.ndarray = None, reset: bool = False):
        """
        Add external displacements to the truss.
        NaN displacements are added to the nodes that are not specified.
        NaN displacements are free nodes.
        ==================
        :param nodes:
        :param disps:
        :param reset: reset external displacements to NaN
        :return:
        """
        if reset:
            self.extDisps *= np.nan
        if nodes is not None:
            #
            __disps = np.zeros((nodes.shape[0], 2, len(self.Steps)))
            __disps = disps[:, :, 0, :] @ np.sin(self.omega[:, np.newaxis] * self.Steps)
            __disps += disps[:, :, 1, :] @ np.cos(self.omega[:, np.newaxis] * self.Steps)
            #
            __vel = np.zeros((nodes.shape[0], 2, len(self.Steps)))
            __vel = (disps[:, :, 0, :] * self.omega) @ np.cos(self.omega[:, np.newaxis] * self.Steps)
            __vel += -(disps[:, :, 1, :] * self.omega) @ np.sin(self.omega[:, np.newaxis] * self.Steps)
            #
            __accel = np.zeros((nodes.shape[0], 2, len(self.Steps)))
            __accel = -(disps[:, :, 0, :] * self.omega ** 2) @ np.sin(self.omega[:, np.newaxis] * self.Steps)
            __accel += -(disps[:, :, 1, :] * self.omega ** 2) @ np.cos(self.omega[:, np.newaxis] * self.Steps)
            #
            self.extDisps[nodes] = __disps
            self.extVel[nodes] = __vel
            self.extAccel[nodes] = __accel

    def calc_sparse_dynamical_matrix(self):
        """
        Calculate the global stiffness matrix of the truss.
        The stiffness matrix is a special form of the structure.
        ==================
        :return:
        """
        __row = np.array(self.Aidx)[:, 0]
        __col = np.array(self.Aidx)[:, 1]
        #
        __Kval = self.Aval * np.repeat(self.weights, 16)
        self.K = ss.coo_matrix((__Kval, (__row, __col)),
                               shape=(self.nodes.shape[0] * 2, self.nodes.shape[0] * 2)).tocsr()
        #
        __Hval = self.Aval * np.repeat(self.eta, 16)
        self.H = ss.coo_matrix((__Hval, (__row, __col)),
                               shape=(self.nodes.shape[0] * 2, self.nodes.shape[0] * 2)).tocsr()
        #
        # __M = np.repeat(np.array([[self.m, self.m, 0, 0,
        #                            self.m, self.m, 0, 0,
        #                            0, 0, self.m, self.m,
        #                            0, 0, self.m, self.m]]), self.elements.shape[0], axis=0).flatten()
        # __Mval = self.Aval * __M
        # self.M = ss.coo_matrix((__Mval, (__row, __col)),
        #                        shape=(self.nodes.shape[0] * 2, self.nodes.shape[0] * 2)).tocsr()
        # Identy matrix
        self.M = ss.eye(self.nodes.shape[0] * 2, format='csr') * self.m

    def dynamicSolve(self):
        __meanDisps = self.extDisps.mean(axis=2)
        self.calDisps = self.extDisps.copy().reshape((-1, len(self.Steps)))
        self.calVel = self.extVel.copy().reshape((-1, len(self.Steps)))
        self.calAccel = self.extAccel.copy().reshape((-1, len(self.Steps)))
        self.calForces = self.extForces.copy().reshape((-1, len(self.Steps)))
        # Find free and fixed nodes using boolean indexing
        free_nodes = np.isnan(__meanDisps.flatten())
        # print(self.extForces[:,:,10])
        # print(free_nodes)
        fixed_nodes = ~free_nodes
        print(free_nodes.sum(), fixed_nodes.sum())
        print(fixed_nodes)
        #
        K_free = self.K[free_nodes][:, free_nodes]
        H_free = self.H[free_nodes][:, free_nodes]
        # M_free = self.M[free_nodes][:, free_nodes]
        #
        F_free = (self.calForces
                  - self.K[:, fixed_nodes] @ self.calDisps[fixed_nodes]
                  - self.H[:, fixed_nodes] @ self.calVel[fixed_nodes]
                  - self.M[:, fixed_nodes] @ self.calAccel[fixed_nodes])[free_nodes]

        # Ode system solver
        plt.plot(self.Steps, F_free[-2, :], 'r')
        plt.show()
        print(self.Steps)

        def ode_system(t, _y):
            # print(F_free[:, int(t / self.dt)] - H_free @ y[-free_nodes.sum():] - K_free @ y[:free_nodes.sum()])
            # print(1/M_free.toarray())
            return np.array([_y[-free_nodes.sum():],
                             (F_free[:, int(t / self.dt)] -
                              H_free @ _y[-free_nodes.sum():] -
                              K_free @ _y[:free_nodes.sum()]) / self.m]).flatten()

        y0 = np.zeros((2 * free_nodes.sum(),))
        # y0[-1] = 1E-3
        sol = scipy.integrate.solve_ivp(ode_system, [self.Steps[0], self.Steps[-1]], y0,
                                        t_eval=self.Steps)
        #
        self.calDisps[free_nodes] = sol.y[:free_nodes.sum()]
        self.calVel[free_nodes] = sol.y[free_nodes.sum():]
        #
        plt.plot(self.Steps, self.calDisps[-2, :], label='x', color='orange')
        plt.show()
        print(self.calDisps[:, -1])

    #
    # def calc_static_Q(self):
    #     """
    #     Calculate the statical charge of the piezo-solid truss structure
    #     ===================
    #     :return: numpy array
    #         Statical charge of the piezo-solid truss structure
    #     """
    #     __bond_stretch, __bond_inner_force = self.calc_bond_stretch()
    #     return self.C * self.extVoltage + self.G * __bond_stretch


def test_transient_piezo_truss():
    """
    Test the static_piezo_truss class
    ===================
    :return:
    """
    # Create a static_piezo_truss object
    static_piezo_truss = TransientPiezoTruss(node_num_4_each_edge=2, path='./Transient_Piezo_Truss',
                                             freq=[10], tspan=[0, 10])
    # # Add voltage to the static_piezo_truss structure
    # static_piezo_truss.add_volForce(voltage=[[1, 1, 0], [0, 0, 0]], seed=1, reset=True)
    # print("static_piezo_truss.extVoltage = \n", static_piezo_truss.extVoltage)
    # print("static_piezo_truss.G = \n", static_piezo_truss.G)
    # print("static_piezo_truss.extForces = \n", static_piezo_truss.extForces)
    # # Add displacement to the static_piezo_truss structure
    static_piezo_truss.add_extTemporalDisps(nodes=np.array([0, 1]),
                                            disps=np.array([
                                                [[[0, ], [0, ]], [[0, ], [0, ]]],
                                                [[[np.nan, ], [np.nan, ]], [[0, ], [0, ]]],
                                                # [[[0, ], [1E-3, ]], [[0, ], [0, ]]],
                                            ]),
                                            reset=True)
    static_piezo_truss.add_extTemporalForces(nodes=np.array([-1]),
                                             forces=np.array([
                                                 [[[0, ], [1, ]], [[0, ], [0, ]]],
                                             ]),
                                             reset=True)
    static_piezo_truss.calc_sparse_dynamical_matrix()
    static_piezo_truss.dynamicSolve()
    # print("static_piezo_truss.extDisps = \n", static_piezo_truss.extDisps)
    # # Solve the static_piezo_truss structure
    # static_piezo_truss.linearSolve()
    # print("static_piezo_truss.calDisps = \n", static_piezo_truss.calDisps)
    # print("static_piezo_truss.calForces = \n", static_piezo_truss.calForces)
    # # Calculate the statical charge of the static_piezo_truss structure
    # static_Q = static_piezo_truss.calc_static_Q()
    # print("static_Q = \n", static_Q)


if __name__ == '__main__':
    test_transient_piezo_truss()
