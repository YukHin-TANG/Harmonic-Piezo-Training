from harmonic_piezo_truss import *
import tqdm
import argparse as arg
import time


def train():
    #
    parser = arg.ArgumentParser(description='ElasticNet')

    parser.add_argument('--path', type=str, default='./log/Harmonic_Piezo_Truss', help='Path to save data.')
    parser.add_argument('--Ne', type=int, default=8, help='Number of edges nodes.')
    parser.add_argument('--m', type=float, default=1.0E-3, help='Mass of the node.')
    parser.add_argument('--k', type=float, default=1.0E+4, help='Weight of the edge.')
    parser.add_argument('--omega', type=float, default=10, help='Frequency of the harmonic excitation.')
    parser.add_argument('--xi', type=float, default=0.2, help='Damping ratio.')
    parser.add_argument('--gamma', type=float, default=1.0E-4, help='Piezoelectric coefficient.')
    parser.add_argument('--C', type=float, default=1.0e-7, help='Capacitance.')
    parser.add_argument('--W_seed', type=int, default=0, help='Network weight random seed.')
    parser.add_argument('--G_seed', type=int, default=0, help='Network gamma random seed.')
    parser.add_argument('--C_seed', type=int, default=0, help='Network capacitance random seed.')
    parser.add_argument('--V_seed', type=int, default=0, help='Network voltage random seed.')
    parser.add_argument('--lr', type=float, default=1E9, help='Learning rate.')
    parser.add_argument('--Vr', type=float, default=10, help='Voltage real part.')
    parser.add_argument('--Vi', type=float, default=1, help='Voltage imaginary part.')
    parser.add_argument('--Ep', type=int, default=20000, help='Number of epochs.')
    parser.add_argument('--Vm', type=float, default=0.0, help='Maximum voltage.')
    parser.add_argument('--A', type=float, default=1E-6, help='Amplitude of the deformation.')
    parser.add_argument('--angle', type=float, default=0, help='Deformation angle.')
    # Convert to dictionary.
    __params = vars(parser.parse_args())
    #
    __path = __params['path']
    for k, v in __params.items():
        if k != 'path':
            __path += k + ':' + str(v) + '/'
    __path += time.strftime("%Y%m%d-%H%M%S") + '/'
    print(__path)
    if os.path.exists(__path) is False:
        os.makedirs(__path)
    # Create a harmonic piezo-solid truss structure with 3 nodes on each edge,
    # and the frequency of the harmonic excitation is 10 rad*Hz, and the damping ratio is 0.8
    harm_piezo_truss = HarmonicPiezoTruss(node_num_4_each_edge=__params['Ne'], m=__params['m'], k=__params['k'],
                                          omega=__params['omega'], xi=__params['xi'], gamma=__params['gamma'],
                                          capacitance=__params['C'], path=__path)

    # Add voltage to the piezo-solid truss structure
    harm_piezo_truss.add_harmVolForce(V=__params['Vr'] + 1j * __params['Vi'], reset=True, seed=__params['V_seed'])
    #
    harm_piezo_truss.add_sinXYDeform(A=__params['A'], angle=__params['angle'], reset=True)
    desired = harm_piezo_truss.extDisps.copy()
    num_slice = 4 * harm_piezo_truss.nn4ee - 4
    #
    harm_piezo_truss.plot_deform(desired, angle=__params['angle'], c='orange', title='Desired Deformation')
    #
    harm_piezo_truss.loss = []
    #
    for _ in tqdm.tqdm(range(__params['Ep'])):
        # Free vibration
        harm_piezo_truss.add_extDisps(np.array([0, 1]), np.array([[0, 0], [np.nan + 1j * np.nan, 0]]), reset=True)
        harm_piezo_truss.linearSolve()
        Q_free = harm_piezo_truss.calc_dynamical_Q()
        harm_piezo_truss.loss += [
            np.linalg.norm(harm_piezo_truss.calDisps[:num_slice] - desired[:num_slice])
            /
            np.linalg.norm(desired[:num_slice])
        ]
        # Forced vibration
        harm_piezo_truss.extDisps = desired.copy()
        harm_piezo_truss.linearSolve()
        Q_forced = harm_piezo_truss.calc_dynamical_Q()

        # Update the voltage
        harm_piezo_truss.extVoltage += __params['lr'] * (Q_forced - Q_free)
        harm_piezo_truss.add_harmVolForce(reset=True)
    #
    harm_piezo_truss.plot_deform(harm_piezo_truss.calDisps,  c='blue', title='Calculated Deformation')
    harm_piezo_truss.plot_loss()
    harm_piezo_truss.save_data()


if __name__ == '__main__':
    train()
