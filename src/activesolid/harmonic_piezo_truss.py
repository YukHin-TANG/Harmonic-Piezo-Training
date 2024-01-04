"""
This file contains the piezo-solid truss class, which is used to calculate the harmonic response of a piezo-solid truss structure


"""
import numpy as np

# Import standard library modules
from truss import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import os


class HarmonicPiezoTruss(Truss):
    """
    HarmonicPiezoTruss class
    """
    G, C, L, Frq, extVoltage = None, None, None, None, None
    loss = None

    def __init__(self, node_num_4_each_edge: int = 3, path=None,
                 gamma: float = 1.0e-4, capacitance: float = 1.0e-8, k: float = 1.0e3,
                 xi: float = 0.9, m: float = 1.0, omega: float = 100.0):
        """
        HarmonicPiezoTruss constructor
        ===================
        :param node_num_4_each_edge: int
            Number of nodes for each edge of the lattice.
        """
        super().__init__(node_num_4_each_edge, path, )
        #
        self.nn4ee = node_num_4_each_edge
        self.m, self.w, self.k = m, omega, k
        self.extForces = np.zeros((self.nodes.shape[0], 2), dtype=complex)
        self.extDisps = np.zeros((self.nodes.shape[0], 2), dtype=complex) * np.nan  # nan means free
        self.weights = self.weights.astype(complex) * k
        self.G = np.ones((self.elements.shape[0],), dtype=complex) * gamma
        self.C = np.ones((self.elements.shape[0],), dtype=complex) * capacitance
        self.eta = 2 * xi * (self.m * self.weights) ** 0.5
        #
        self.extVoltage = np.zeros((self.nodes.shape[0],), dtype=complex)
        # self.calc_sparse_stiffness_matrix()
        self.calc_sparse_harm_stiffness_matrix()

    def add_harmVolForce(self, V: complex = None, reset: bool = False, seed: int = 0):
        """
        Add voltage to the piezo-solid truss structure
        ===================
        :param V: complex
            => v = V * exp(jwt) = [V1 + j*V2] * [cos(wt) + j*sin(wt)]
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
        if V is not None:
            self.extVoltage = np.ones((self.elements.shape[0],)) * V
            if seed != 0:
                np.random.seed(seed)
                self.extVoltage.real *= np.random.uniform(0, 2, (self.elements.shape[0],))
                np.random.seed(seed + 1)
                self.extVoltage.imag *= np.random.uniform(0, 2, (self.elements.shape[0],))
        # Calculate the direction of the edges
        vec = self.nodes[self.elements[:, 1], :] - self.nodes[self.elements[:, 0], :]
        norm_vec = np.linalg.norm(vec, axis=1)[:, np.newaxis]
        vec /= norm_vec
        # Create an empty array to store the voltage force
        voltage_force = np.empty((self.elements.shape[0], 2), dtype=complex)
        # Assign the voltage force to the array
        voltage_force[:, 0] = vec[:, 0] * (self.extVoltage * self.G)
        voltage_force[:, 1] = vec[:, 1] * (self.extVoltage * self.G)
        # Add the voltage force to the external force
        for n, e in enumerate(self.elements):
            self.extForces[e[0], :] -= voltage_force[n, :]
            self.extForces[e[1], :] += voltage_force[n, :]
        # print(self.extForces)

    def calc_sparse_harm_stiffness_matrix(self):
        """
        Calculate the global stiffness matrix of the truss.
        The stiffness matrix is a special form of the structure.
        ==================
        :return:
        """
        # print(self.Aval)
        __W2 = self.w * self.eta * 1j + self.weights
        __W1 = -self.m * self.w ** 2 + __W2
        __W = []
        for w1, w2 in zip(__W1, __W2):
            __W += [
                w1, w1, w2, w2,
                w1, w1, w2, w2,
                w2, w2, w1, w1,
                w2, w2, w1, w1,
            ]

        __Kval = self.Aval * np.array(__W, dtype=complex)
        __row = np.array(self.Aidx)[:, 0]
        __col = np.array(self.Aidx)[:, 1]
        self.K = ss.coo_matrix((__Kval, (__row, __col)),
                               shape=(self.nodes.shape[0] * 2, self.nodes.shape[0] * 2)).tocsr()
        # print(self.K)

    def calc_bond_stretch(self):
        """
        Calculate the stress of each element.
        :return: bond stretch, bond inner force
        """
        __idx = self.elements[:, 0]
        __idy = self.elements[:, 1]
        __disps_diff = self.calDisps[__idy, :] - self.calDisps[__idx, :]  # displacement difference
        # print("_______________")
        # print(__disps_diff)
        #
        __bond_vector = self.nodes[__idy, :] - self.nodes[__idx, :]  # bond vector
        __bond_length = np.linalg.norm(__bond_vector, axis=1)  # bond length
        __bond_direction = __bond_vector / __bond_length[:, np.newaxis]  # bond direction
        #
        __bond_stretch = np.sum(__disps_diff * __bond_direction, axis=1)  # bond stretch
        #
        __bond_inner_force = __bond_stretch * self.weights  # bond inner force
        #
        return __bond_stretch, __bond_inner_force

    def calc_dynamical_Q(self):
        """
        Calculate the statical charge of the piezo-solid truss structure
        ===================
        :return: numpy array
            Statical charge of the piezo-solid truss structure
        """
        __bond_stretch, __bond_inner_force = self.calc_bond_stretch()
        return self.C * self.extVoltage + self.G * __bond_stretch

    def add_sinXYDeform(self, A: float = 1e-6, angle: float = 0.0, reset: bool = False):
        """
        Add a sinX deformation to the mechanical meta-material.
        :return:
        """
        if reset:
            self.extDisps *= np.nan # nan means free
        __element_num = 4 * (self.nn4ee - 1)
        self.extDisps[:__element_num] = 0.0
        # __element_num = self.nodes.shape[0]
        for n in range(__element_num):
            x, y = self.nodes[n, 0], self.nodes[n, 1]
            if abs(x) == 0.5:
                self.extDisps[n] += np.array([
                    np.sin(0.5 * np.pi * (y + 0.5)),
                    0
                ])
            if y == 0.5:
                self.extDisps[n] += np.array([
                    0,
                    np.sin(2 * np.pi * (x + 0.5))
                ])

            if y == -0.5:
                self.extDisps[n] += np.array([
                    0,
                    -np.sin(np.pi * (x + 0.5))
                ])
        #
        self.extDisps *= A * np.exp(1j * angle / 180 * np.pi)

    def plot_extDeform(self, Disps, angle):
        Disps = Disps * np.exp(-1j * angle / 180 * np.pi)
        s = np.nanmax(np.linalg.norm(np.abs(Disps), axis=1))
        #
        plt.figure(figsize=(10 / 2.54, 10 / 2.54))
        plt.gca().set_aspect("equal")
        plt.xlim(-0.7, 0.7)
        plt.ylim(-0.7, 0.7)
        __weights = np.abs(self.weights)
        for (i, j), __k in zip(self.elements, __weights):
            plt.plot(self.nodes[[i, j], 0], self.nodes[[i, j], 1],
                     color="green", linewidth=1.5 * __k / self.k, zorder=1)
        plt.plot(self.nodes[:, 0], self.nodes[:, 1], "o", markersize=4, zorder=1)
        plt.quiver(self.nodes[:, 0], self.nodes[:, 1], Disps[:, 0].real, Disps[:, 1].real,
                   color='orange', scale=8 * s, zorder=2, linewidth=0.5, edgecolor='orange')
        plt.title("External Deformation")
        circ = plt.Circle((-0.57, -0.6), 0.07, color='grey', fill=False)
        plt.gcf().gca().add_artist(circ)
        plt.quiver(-0.57, -0.6, 0.07 * np.cos(angle), 0.07 * np.sin(angle), color='red', zorder=2, linewidth=1, )
        plt.plot([-0.57, -0.57 + 0.07], [-0.6, -0.6], 'k--', linewidth=1)

        if self.path is not None:
            if not os.path.exists(self.path):
                os.makedirs(self.path)
            plt.savefig(self.path + "extDeform.svg", dpi=300)

    def plot_harmVar(self, var: np.ndarray, title: str = None, cm: str = 'RdBu_r'):
        """
        Plot the voltage of the piezo-solid truss structure
        ===================
        :param var: numpy array
            Variable to be plotted
        :param title: str
            Title of the plot
        :param cm: str
            Colormap of the plot
        :return:
        """
        plt.figure(figsize=(10 / 2.54, 10 / 2.54))
        plt.title(title)
        plt.gca().set_aspect("equal")
        plt.xlim(-0.7, 0.7)
        plt.ylim(-0.7, 0.7)
        # plot lattice and colorbar for external thermal
        vmin = np.min(var)
        vmax = np.max(var)
        #  Center the colorbar on 0
        if vmin > 0:
            norm = mpl.colors.TwoSlopeNorm(vmin=-vmin, vcenter=0, vmax=vmax)
        elif vmax < 0:
            norm = mpl.colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=-vmax)
        else:
            norm = mpl.colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
        cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.colormaps[cm])
        cmap.set_array([])
        for (i, j), w in zip(self.elements, var):
            plt.plot(self.nodes[[i, j], 0], self.nodes[[i, j], 1], linewidth=2, color=cmap.to_rgba(w))
        # Add colorbar
        plt.colorbar(cmap, ax=plt.gca())
        if not os.path.exists(self.path):
            os.makedirs(self.path)
        plt.savefig(self.path + title + ".svg", dpi=300)

    def plot_harmVoltage(self):
        """
        Plot the voltage magnitude and phase angle of the piezo-solid truss structure
        ===================
        :return:
        """

        self.plot_harmVar(var=np.abs(self.extVoltage), title="Voltage Amplitude (V)")
        self.plot_harmVar(var=np.angle(self.extVoltage, deg=True), title=r"Voltage Phase Angle ( $^o$ )")

    def plot_loss(self):
        """
        Plot the loss of the piezo-solid truss structure
        ===================
        :return:
        """
        self.plot_harmVar(var=self.loss, title="Loss")

    def save_data(self):
        """
        Save the data of the piezo-solid truss structure
        ===================
        :return:
        """
        print("Saving data...")
        np.savez(self.path + "data.npz", nodes=self.nodes,
                 elements=self.elements,
                 weights=self.weights,
                 extForces=self.extForces,
                 extDisps=self.extDisps,
                 extVoltage=self.extVoltage,
                 G=self.G,
                 C=self.C,
                 L=self.L,
                 Frq=self.Frq,
                 k=self.k,
                 m=self.m,
                 w=self.w,
                 nn4ee=self.nn4ee,
                 Aidx=self.Aidx,
                 Aval=self.Aval,
                 K=self.K,
                 A=self.A)
        print("Save loss data...")
        np.savez(self.path + "loss.npz", loss=self.loss)
        print("Data saved.")


def test():
    """
    Test function
    ===================
    :return:
    """
    # Create a piezo-solid truss structure
    piezo_truss = HarmonicPiezoTruss(node_num_4_each_edge=3, omega=10, xi=1.2)
    # Add voltage to the piezo-solid truss structure
    piezo_truss.add_harmVolForce(V=1 + 1j, reset=True, seed=0)
    print("piezo_truss.extVoltage = \n", piezo_truss.extVoltage)
    print("piezo_truss.G = \n", piezo_truss.G)
    print("piezo_truss.extForces = \n", piezo_truss.extForces)
    # Add displacement to the piezo-solid truss structure
    piezo_truss.add_extDisps(nodes=[0, 1], disps=[[0, 0], [np.nan + 1j * np.nan, 0]], reset=True)
    print("piezo_truss.extDisps = \n", piezo_truss.extDisps)
    # Solve the piezo-solid truss structure
    piezo_truss.linearSolve()
    print("piezo_truss.calDisps = \n", piezo_truss.calDisps)
    print("piezo_truss.calForces = \n", piezo_truss.calForces)
    # Calculate the statical charge of the piezo-solid truss structure
    dynamical_Q = piezo_truss.calc_dynamical_Q()
    print("dynamical_Q = \n", dynamical_Q)


if __name__ == '__main__':
    test()
