## -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 17:55:46 2022

@author: cosmi
"""

import numpy as np

from qutip import Qobj, rand_dm, fidelity, displace, qdiags, qeye, expect
from qutip.states import coherent, coherent_dm, thermal_dm, fock_dm
from qutip.visualization import plot_wigner, hinton, matrix_histogram
from qutip.wigner import qfunc
import qutip

import matplotlib.pyplot as plt

import matplotlib.animation as ani
import matplotlib.colors as clr
import types

import IPython.display as ipyd


# some pretty printing and animation stuff
from IPython.display import clear_output

"""
Iterative Maximum Likelihood estimation based on photon number counting.
"""


def measure_population(beta, rho):
    """
    Measures the photon number statistics for state rho when displaced
    by angle alpha.
    
    Parameters
    ----------    
    alpha: np.complex
        A complex displacement.

    rho:
        The density matrix as a QuTiP Qobj (`qutip.Qobj`)

    Returns
    -------
    population: ndarray
        A 1D array for the probabilities for populations.
    """
    hilbertsize = rho.shape[0]
    # Apply a displacement to the state and then measure the diagonals.

    D = displace(hilbertsize, beta)
    rho_disp = D*rho*D.dag()
    populations = np.real(np.diagonal(rho_disp.full()))
    return populations

def roperator(beta, rho, measured):
    """
    Calculates the iterative ratio operator for measured probability for photons
    (j) to the analytical prediction for some rho.

    Parameters
    ----------
    beta: list_like
        A list of the displacements that were applied to the state before
        measurement.

    rho: `qutip.Qobj`
        The current estimate of the density matrix.

    measured: list_like
        A list of the measurement statistics (diagonal terms) for each beta.

    Returns
    -------
    R: `qutip.Qobj`
        The iterative operator which we are going to apply for state
        reconstruction.
    """
    hilbert_size = rho.shape[0]

    # initialize an empty operator and build it

    R = 0*qeye(hilbert_size)
    calculated_measurements = measure_population(beta, rho)

    for n in range(hilbert_size):
        op = fock_dm(hilbert_size, n)
        D = displace(hilbert_size, beta)
        displaced_D = D.dag()*op*D
        ratio = measured[n]/(calculated_measurements[n] + 1e-6)
        displaced_D = ratio*displaced_D
        R += displaced_D

    return R


hilbert_size = 4
alpha_range = 1.9

alphas = np.array([alpha_range, -alpha_range - 1j*alpha_range,
                   -alpha_range + 1j*alpha_range])

rho_true = sum([coherent_dm(hilbert_size, a) for a in alphas])/3


betas = [1.7, -2, 2.2j, -2.1 - 2.4j, -2 + 2j]
measured_populations = [measure_population(b, rho_true) for b in betas]
width = 1

random_rho = rand_dm(hilbert_size)
hinton(random_rho)
plt.show()


# visualize H

lbls_list = [[str(d) for d in range(hilbert_size)], ["u", "d"]]

xlabels = []

fig, ax = matrix_histogram(random_rho, xlabels, xlabels, limits=[-1,1])

ax.view_init(azim=-55, elev=45)

plt.show()




# Prepare figure
fname = "rabi_oscill.gif"

fig, ax = plt.subplots(1, 3, figsize=(15, 5))
indices = np.arange(hilbert_size)

plot_wigner(random_rho, fig, ax[0])
ax[0].scatter(np.real(betas), np.imag(betas), marker="x")
ax[0].set_title("Random inital state wigner function")
for i in range(len(betas)):
    ax[1].bar(indices, measured_populations[i],
              label = r"$beta = {}$".format(i), width=(i+1)/12)

ax[1].set_title("Population measurement statistics")
ax[1].set_xlabel("n")
ax[1].set_ylabel("Photon number probability")


plot_wigner(rho_true, fig, ax[2])
ax[2].scatter(np.real(betas), np.imag(betas), marker="x")
ax[2].set_title("Target state wigner function")
plt.show()




Nfock = 10

psi0 = qutip.coherent(Nfock, 0)
Tsim = np.arange(0, 2.5, 1e-3)
Ntraj = 25


SYS_num = measured_populations; SYS_num

Nfock = 10
#SYS_num.dimension = Nfock

psi0 = qutip.coherent(Nfock, 0)
Tsim = np.arange(0, 2.5, 1e-3)
Ntraj = 25

H_num, L_num = lbls_list
qmc = qutip.mcsolve(H_num, psi0, Tsim, L_num, [], ntraj=Ntraj)

N_num = qutip.num(Nfock)

plt.figure(figsize=(12,5))
plt.xlabel("Time", fontsize=14); plt.ylabel("Photon Number", fontsize=14)
plt.tick_params(labelsize=14)
plt.plot(qmc.times, np.ones(qmc.times.size), "--k")

for traj in qmc.states:
    plt.plot(qmc.times, qutip.expect(N_num, traj), "b", alpha=0.2)
    
    
Tsim = np.arange(0, 6, 1e-2)

qme = qutip.mesolve(H_num, psi0, Tsim, L_num, [])
rho_ss = qutip.steadystate(H_num, L_num)


plt.figure(figsize=(12,5))
plt.xlabel("Time", fontsize=14); plt.ylabel("Photon Number", fontsize=14)
plt.tick_params(labelsize=14)

plt.plot(qme.times, qutip.expect(rho_ss, N_num)*np.ones(qme.times.size), "--k")
plt.plot(qme.times, qutip.expect(qme.states, N_num)) # Compute steady-state of the system


# Set Wigner function scale
xvec = np.linspace(-7, 7, 150)

# Parallelized Wigner computation
def compute_wigner(rho):
    return qutip.wigner(rho, xvec, xvec)
W_list = qutip.parfor(compute_wigner, qme.states)



# Prepare figure
fname = "rabi_oscill.gif"
fig, ax = plt.subplots(1,2, figsize=(12,5))

# Plot steady state
ax[1].contourf(xvec, xvec, qutip.wigner(rho_ss, xvec, xvec), 100, norm=clr.Normalize(-0.25,0.25), cmap=plt.get_cmap("RdBu"))
ax[1].set_aspect("equal"); ax[1].set_title("Steady State", fontsize=14); ax[1].tick_params(labelsize=14)
ax[1].set_xlabel(r"$x$", fontsize=14); ax[1].set_ylabel(r"$p$", fontsize=14)

# Animate evolution
def animate(n):
    ax[0].cla(); ax[0].set_aspect("equal"); ax[0].tick_params(labelsize=14)
    ax[0].set_title("Time: %.2f"%(qme.times[n]), fontsize=14);
    ax[0].set_xlabel(r"$x$", fontsize=14); ax[0].set_ylabel(r"$p$", fontsize=14)
    im = ax[0].contourf(xvec, xvec, W_list[n], 100, norm=clr.Normalize(-0.25,0.25), cmap=plt.get_cmap("RdBu"))
    def setvisible(self, vis): # Work-around for visibility bug in contourf
        for c in self.collections: c.set_visible(vis)
    im.set_visible = types.MethodType(setvisible, im)
anim = ani.FuncAnimation(fig, animate, frames=len(qme.times))
anim.save(fname, writer="imagemagick", fps=20)
plt.close()


ipyd.Image(url=fname)

