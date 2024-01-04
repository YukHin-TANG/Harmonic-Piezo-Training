"""
This file contains the piezo-solid truss class, which is used to calculate the harmonic response of a piezo-solid truss structure
"""
import numpy as np

# Import standard library modules
from truss import *
import pymunk
import matplotlib.pyplot as plt


class HarmonicPiezoTruss(Truss):
    """
    HarmonicPiezoTruss class
    """
    B, C, Frq, extVoltage = None, None, None, None
    space = None

    def __init__(self, node_num_4_each_edge: int = 2, path=None,
                 gamma: float = 1.0e-4, capacitance: float = 1.0e-8, omega: float = 1.0e3,
                 xi: float = 1.0e-2, m: float = 1.0):
        """
        HarmonicPiezoTruss constructor

        ===================
        :param node_num_4_each_edge: int
            Number of nodes for each edge of the lattice.
        :param path: str
            Path to the lattice file
        :param gamma: float
            Piezoelectric coupling coefficient
        :param capacitance: float
            Capacitance
        :param omega: float
            Angular frequency
        :param xi: float
            Damping ratio
        :param m: float
            Mass
        """
        super().__init__(node_num_4_each_edge, path, )
        #
        self.B = np.ones((self.elements.shape[0],)) * gamma
        self.C = np.ones((self.elements.shape[0],)) * capacitance
        self.xi, self.m, self.Omega = xi, m, omega
        self.weights = np.ones((self.elements.shape[0],)) * 10000
        #
        self.extVoltage = np.zeros((self.nodes.shape[0],))
        self.add_dynamic_model()
        # self.calc_sparse_stiffness_matrix()

    def add_dynamic_model(self):
        """
        Add dynamic model to the piezo-solid truss structure with pymunk
        :return:
        """
        print('Add dynamic model to the piezo-solid truss structure')
        # Calculate the stiffness matrix
        # Initialize space
        self.space = pymunk.Space()
        self.space.gravity = (0, 0)  # zero gravity
        # Create bodies
        __r = 1e-3
        for pos in self.nodes:
            moment = pymunk.moment_for_circle(self.m, 0, __r)  # moment of inertia for the mass
            bt = pymunk.Body.DYNAMIC
            if pos[0] == pos[1] == -0.5:
                bt = pymunk.Body.STATIC
            body = pymunk.Body(self.m, moment, body_type=bt)
            body.position = pos.tolist()
            shape = pymunk.Circle(body, __r)
            self.space.add(body, shape)

        # Create links
        for n, (i, j) in enumerate(self.elements):
            __length = np.linalg.norm(self.nodes[i, :] - self.nodes[j, :])
            __beta = 2 * self.xi * (self.weights[n] / self.m) ** 0.5
            link = pymunk.DampedSpring(self.space.bodies[i], self.space.bodies[j],
                                       (0, 0), (0, 0),
                                       __length, self.weights[n], __beta)
            self.space.add(link)

    def add_nodeForce(self):
        pass

    def dynamic_sim(self, dt: float = 1e-3, steps: int = 1000):
        """
        Dynamic simulation
        :param dt: float
            Time step
        :param steps: int
            Number of steps
        :return:
        """
        force = (1, 0)
        positions = []
        # Add a y-direction displacement constraint to the second node
        body = self.space.bodies[1]
        pivot = pymunk.PivotJoint(body, self.space.static_body, (0, 0), body.position)
        self.space.add(pivot)

        for i in range(steps):
            self.space.bodies[-1].apply_force_at_local_point(force)
            self.space.step(dt)
            # print(self.space.bodies[1].position)
            positions += [self.space.bodies[-1].position]

        plt.plot(np.arange(steps) * dt, np.array(positions)[:, 0], ':', label='x')
        plt.plot(np.arange(steps) * dt, np.array(positions)[:, 1], '-.', label='y')
        plt.legend()
        plt.show()


def test():
    """
    Test function
    ===================
    :return:
    """
    # Create a piezo-solid truss structure
    piezo_truss = HarmonicPiezoTruss(node_num_4_each_edge=3)
    piezo_truss.dynamic_sim(steps=30000, dt=1e-3)


if __name__ == '__main__':
    test()
