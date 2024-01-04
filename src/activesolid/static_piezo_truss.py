"""
This file contains the piezo-solid truss class, which is used to calculate the statics of a piezo-solid truss structure
"""

# Import standard library modules
from truss import *


class StaticPiezoTruss(Truss):
    """
    static_piezo_truss class
    ========================
        This class is used to calculate the statics of a piezo-solid truss structure
    """
    G, C, extVoltage = None, None, None

    def __init__(self, node_num_4_each_edge: int = 3, path=None, gamma: float = 1.0e-4, capacitance: float = 1.0e-8):
        """
        static_piezo_truss constructor
        ===================
        :param node_num_4_each_edge: int
            Number of nodes for each edge of the lattice.
        """
        super().__init__(node_num_4_each_edge, path, )  # call the constructor of the parent class
        #
        self.G = np.ones((self.elements.shape[0],)) * gamma
        self.C = np.ones((self.elements.shape[0],)) * capacitance
        self.extVoltage = np.zeros((self.nodes.shape[0],))

    def add_volForce(self, voltage: float = None, reset: bool = False, seed: int = 0):
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
        if reset:
            self.extForces *= 0
        if voltage is not None:
            if seed != 0:
                np.random.seed(seed)
                self.extVoltage = np.random.uniform(0, 2, (self.elements.shape[0],)) * voltage
            else:
                self.extVoltage = np.ones((self.elements.shape[0],)) * voltage
        # Calculate the direction of the edges
        vec = self.nodes[self.elements[:, 1], :] - self.nodes[self.elements[:, 0], :]
        norm_vec = np.linalg.norm(vec, axis=1)[:, np.newaxis]
        vec /= norm_vec
        # Create an empty array to store the voltage force
        voltage_force = np.empty((self.elements.shape[0], 2))
        # Assign the voltage force to the array
        voltage_force[:, 0] = vec[:, 0] * (self.extVoltage * self.G)
        voltage_force[:, 1] = vec[:, 1] * (self.extVoltage * self.G)
        # Add the voltage force to the external force
        for n, e in enumerate(self.elements):
            self.extForces[e[0], :] -= voltage_force[n, :]
            self.extForces[e[1], :] += voltage_force[n, :]

    def calc_static_Q(self):
        """
        Calculate the statical charge of the piezo-solid truss structure
        ===================
        :return: numpy array
            Statical charge of the piezo-solid truss structure
        """
        __bond_stretch, __bond_inner_force = self.calc_bond_stretch()
        return self.C * self.extVoltage + self.G * __bond_stretch


def test_static_piezo_truss():
    """
    Test the static_piezo_truss class
    ===================
    :return:
    """
    # Create a static_piezo_truss object
    static_piezo_truss = StaticPiezoTruss(node_num_4_each_edge=3, path='./Static_Piezo_Truss')
    # Add voltage to the static_piezo_truss structure
    static_piezo_truss.add_volForce(voltage=1.0, seed=0, reset=True)
    print("static_piezo_truss.extVoltage = \n", static_piezo_truss.extVoltage)
    print("static_piezo_truss.G = \n", static_piezo_truss.G)
    print("static_piezo_truss.extForces = \n", static_piezo_truss.extForces)
    # Add displacement to the static_piezo_truss structure
    static_piezo_truss.add_extDisps(nodes=[0, 1], disps=[[0, 0], [np.nan, 0]], reset=True)
    print("static_piezo_truss.extDisps = \n", static_piezo_truss.extDisps)
    # Solve the static_piezo_truss structure
    static_piezo_truss.linearSolve()
    print("static_piezo_truss.calDisps = \n", static_piezo_truss.calDisps)
    print("static_piezo_truss.calForces = \n", static_piezo_truss.calForces)
    # Calculate the statical charge of the static_piezo_truss structure
    static_Q = static_piezo_truss.calc_static_Q()
    print("static_Q = \n", static_Q)


if __name__ == '__main__':
    test_static_piezo_truss()
