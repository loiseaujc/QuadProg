import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

    # Data.
    k, p, v, u = np.loadtxt("mpc_response.dat", unpack=True)

    # Figure.
    fig, ax = plt.subplots(2, 1, figsize=(9, 3), sharex=True)

    ax[0].plot(k, p, label="p")
    ax[0].plot(k, v, label="v")
    ax[0].axhline(0, c="k", ls="--")
    ax[0].set(ylim=(-abs(p).max(), abs(p).max()))
    ax[0].legend(loc="upper right", ncols=2)

    ax[1].plot(k, u, label="u")
    ax[1].axhline(0, c="k", ls="--")
    ax[1].legend(loc="upper right")
    ax[1].set(xlim=(0, k.max()), xlabel="Time")
    ax[1].set(ylim=(-np.ceil(1.1*abs(u).max()), np.ceil(1.1*abs(u).max())))

    plt.show()
