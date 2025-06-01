import numpy as np
import glob
import os
import matplotlib.pyplot as plt
from cal_order_para import read_op_series


def PD_gamma_T(L, rho0, Dt, sigma=0, h=0.05, seed=3001):
    folder = f"../data/PD/L{L}"
    if sigma == 0:
        gamma_arr = np.array([0, 0.002, 0.004, 0.006, 0.008, 0.01, 0.012, 0.014, 0.016, 0.018, 0.02,
                            0.022, 0.024, 0.026, 0.028, 0.03])
        T_arr = np.array([0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2,
                        0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4,
                        0.42, 0.44, 0.46, 0.48, 0.5])
    elif sigma == 0.05:
        gamma_arr = np.array([0, 0.002, 0.004, 0.006, 0.008, 0.01, 0.012, 0.014, 0.016, 0.018, 0.02,
                            0.022, 0.024, 0.026, 0.028, 0.03])
        T_arr = np.array([0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2,
                        0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4])
    phi_arr = np.zeros((T_arr.size, gamma_arr.size))

    for j, T in enumerate(T_arr):
        for i, gamma in enumerate(gamma_arr):
            t, phi = read_op_series(L, rho0, Dt, T, gamma, sigma, h, seed, folder)
            ncut = 10000
            phi_arr[j, i] = np.mean(phi[ncut:])

    fig, ax = plt.subplots(1, 1, constrained_layout=True)

    dx = gamma_arr[1] - gamma_arr[0]
    dy = T_arr[1] - T_arr[0]
    extent = [gamma_arr[0]-dx/2, gamma_arr[-1]+dx/2, T_arr[0]-dy/2, T_arr[-1]+dy/2]
    im = ax.imshow(phi_arr, origin="lower", extent=extent, aspect="auto")
    ax.set_xlabel(r"$\gamma$", fontsize="x-large")
    ax.set_ylabel(r"$T$", fontsize="x-large")
    cb = fig.colorbar(im)
    cb.set_label(r"Global Polarity", fontsize="x-large")
    # fig.suptitle(r"$L=%g, D = 0, \sigma=0$", fontsize="x-large")
    plt.show()
    plt.close()

if __name__ == "__main__":
    L = 128
    rho0 = 1.5
    Dt = 0.01
    sigma = 0.05
    PD_gamma_T(L, rho0, Dt, sigma=sigma)