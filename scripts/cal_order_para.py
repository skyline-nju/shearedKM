import numpy as np
import glob
import os
import matplotlib.pyplot as plt
from add_line import add_line


def read_op_series(L, rho0, Dt, T, gamma, sigma, h, seed, folder, ret_theta=False):
    pat = f"{folder}/L{L:d}_{L:d}_r{rho0:g}_Dt{Dt:g}_T{T:g}_g{gamma:g}_s{sigma:g}_h{h:g}_S{seed:d}_t*.dat"
    files = glob.glob(pat)
    lines_dict = {}
    n_lines = 0
    for file in files:
        t_beg = int((os.path.basename(file).rstrip(".dat").split("_")[-1]).lstrip("t"))
        with open(file, "r") as fin:
            lines = fin.readlines()
            lines_dict[t_beg] = lines
            n_lines += len(lines)
    phi_arr = np.zeros(n_lines)
    print("find", n_lines, "lines")
    if ret_theta:
        theta_arr = np.zeros_like(phi_arr)
    i = 0
    for t_beg in sorted(lines_dict.keys()):
        print(t_beg)
        print(file)
        lines = lines_dict[t_beg]
        for line in lines:
            s = line.rstrip("\n").split("\t")
            try:
                phi_arr[i] = float(s[1])
                if ret_theta:
                    theta_arr[i] = float(s[2])
            except IndexError:
                phi_arr[i] = phi_arr[i-1]
                if ret_theta:
                    theta_arr[i] = theta_arr[i-1]
            except ValueError:
                phi_arr[i] = phi_arr[i-1]
                if ret_theta:
                    theta_arr[i] = theta_arr[i-1]

            i += 1
    t_arr = (np.arange(phi_arr.size) + 1) * 100 * h
    if not ret_theta:
        return t_arr, phi_arr
    else:
        return t_arr, phi_arr, theta_arr


def plot_time_series(L, rho0, Dt, T, gamma, sigma, h, seed):
    folder = f"D:/code/shearedKM/data/PD/L{L}"

    t_arr, phi_arr, theta_arr = read_op_series(L, rho0, Dt, T, gamma, sigma, h, seed, folder, True)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 4), constrained_layout=True)
    ax1.plot(t_arr, phi_arr)
    ax2.plot(t_arr, theta_arr)
    plt.show()
    plt.close()


def phi_vs_T_varied_gamma(L, gamma_arr, rho0, Dt, sigma=0, h=0.05, seed=3001):
    folder = f"../data/PD/L{L}"
    T_arr = np.array([0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2,
                      0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4,
                      0.42, 0.44, 0.46, 0.48, 0.5])
    fig, ax = plt.subplots(1, 1, constrained_layout=True)

    for gamma in gamma_arr:
        phi_arr = np.zeros(T_arr.size)

        for j, T in enumerate(T_arr):
            t, phi = read_op_series(L, rho0, Dt, T, gamma, sigma, h, seed, folder)
            if T > 0.4 or gamma > 0.02:
                ncut = 1000
            else:
                ncut = 10000
            phi_arr[j] = np.mean(phi[ncut:])

        ax.plot(T_arr, phi_arr, "-o", label=r"$\gamma=%g$" % gamma)
    ax.set_xlabel(r"$T$", fontsize="x-large")
    ax.set_ylabel(r"$\langle \phi\rangle$", fontsize="x-large")
    fig.legend()
    plt.show()
    plt.close()


def phi_vs_T_varied_L(L_arr, gamma, rho0, Dt, sigma=0, h=0.05, seed=3001):
    fig, ax = plt.subplots(1, 1, constrained_layout=True)

    T_arr = np.array([0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2,
                    0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4,
                    0.42, 0.44, 0.46, 0.48, 0.5])

    for L in L_arr:
        folder = f"../data/PD/L{L}"
        phi_arr = np.zeros(T_arr.size)
        for j, T in enumerate(T_arr):
            t, phi = read_op_series(L, rho0, Dt, T, gamma, sigma, h, seed, folder)
            if phi.size > 10000:
                ncut = 10000
            else:
                ncut = 1000
            phi_arr[j] = np.mean(phi[ncut:])
        ax.plot(T_arr, phi_arr, "-o", label=r"$L=%d$" % L)
    ax.set_xlabel(r"$T$", fontsize="x-large")
    ax.set_ylabel(r"$\langle \phi\rangle$", fontsize="x-large")
    fig.legend()
    plt.show()
    plt.close()


if __name__ == "__main__":
    # L = 128
    # gamma_arr = [0, 0.002, 0.01]
    # rho0 = 1.5
    # Dt = 0.01
    # phi_vs_T_varied_gamma(L, gamma_arr, rho0, Dt)

    # L_arr = [32, 64, 128, 256, 512]
    # gamma = 0.
    # rho0 = 1.5
    # Dt = 0.01
    # phi_vs_T_varied_L(L_arr, gamma, rho0, Dt, sigma=0)


    L = 128
    gamma = 0.002
    plot_time_series(L, 1.5, 0.01, 0.3, gamma=gamma, sigma=0, h=0.05, seed=3001)