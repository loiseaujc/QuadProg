import numpy as np
import matplotlib.pyplot as plt

def equicorr(r):
    cov = np.zeros((3, 3)) ; cov[:, :] = r
    np.fill_diagonal(cov, 1)
    return cov

def plot_efficient_frontier(expected_return, expected_std, portfolios, mu, cov):

    fig, ax = plt.subplots(1, 1, figsize=(6, 4))

    # Efficient frontier.
    ax.plot(expected_std, expected_return, lw=2, c="k")
    ax.set(xlabel="Standard deviation", ylabel="Expected return")
    ax.set(xlim=(0.5, 1), ylim=(1, 2))

    # Compute expected returns and standard deviations of the random portfolios.
    r = portfolios @ mu
    std = np.array([np.sqrt(x.T @ cov @ x) for x in portfolios])

    # Compute Sharpe ratio.
    sharpe = r / std

    # Plot random portfolios.
    h = ax.scatter(std, r, c=sharpe, s=0.2)
    fig.colorbar(h, label="Sharpe ratio")

if __name__ == "__main__":

    mu = np.ones(3); mu[0] = 2.0
    
    #---------------------------------------
    #-----     Uncorrelated assets     -----
    #---------------------------------------

    # Load data.
    data = np.loadtxt("uncorrelated_efficient_frontier.dat")

    # Covariance matrix.
    cov = equicorr(0.0)

    # Generate random portfolios.
    portfolios = np.random.rand(10000, 3) ; portfolios /= portfolios.sum(axis=1).reshape(-1, 1)
    plot_efficient_frontier(data[:, 0], data[:, 1], portfolios, mu, cov)

    plt.savefig("Pareto_font_uncorrelated_assets.png", bbox_inches="tight", dpi=1200)

    #-----------------------------------------
    #-----     0.5-correlated assets     -----
    #-----------------------------------------

    # Load data.
    data = np.loadtxt("mild_correlation_efficient_frontier.dat")

    # Covariance matrix.
    cov = equicorr(0.5)

    # Generate random portfolios.
    portfolios = np.random.rand(10000, 3) ; portfolios /= portfolios.sum(axis=1).reshape(-1, 1)
    plot_efficient_frontier(data[:, 0], data[:, 1], portfolios, mu, cov)

    plt.savefig("Pareto_font_mild_correlation_assets.png", bbox_inches="tight", dpi=1200)

    #-----------------------------------------
    #-----     0.8-correlated assets     -----
    #-----------------------------------------

    # Load data.
    data = np.loadtxt("strong_correlation_efficient_frontier.dat")

    # Covariance matrix.
    cov = equicorr(0.8)

    # Generate random portfolios.
    portfolios = np.random.rand(10000, 3) ; portfolios /= portfolios.sum(axis=1).reshape(-1, 1)
    plot_efficient_frontier(data[:, 0], data[:, 1], portfolios, mu, cov)

    plt.savefig("Pareto_font_strong_correlation_assets.png", bbox_inches="tight", dpi=1200)
    
    plt.show()
