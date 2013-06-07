#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

# This is an example of sampling of f(x) from a prior and posterior Gaussian Process

f_x_samples = 15

def k_squared_exp(x0, x1, v0, l):
    return v0*np.exp(-1.0/(2.0 * l**2) * (x0-x1)**2.0)

def k_rational_quadratic(x0, x1, alpha, l):
    return (1.0 + 1.0/(2.0 * alpha * l**2) * (x0-x1)**2.0 ) ** (-alpha)

class GaussianProcess:
    """ A class which supports sampling from a Gaussian Process
        prior and posterior.
    """
    def __init__(self, m, k, l = 1.0, v0 = 1.0, alpha = 0.1, sigma = 1.0):
        """
            Parameters
            ----------

            m : the mean function. Callable with 1 argument.

            k : the covariance function. "Squared_exponential",
               "rational_quadratic", or your own callable with 2 arguments
               that can be used as a covariance function (positive-semidefinite)

            l : length scale of the covariance function (if one of the presets
                is used, otherwise ignored)

            v0 : controls the pointwise prior variability

            alpha : the negative exponent in the rational quadratic cov. function
                
            sigma : standard deviation of the likelihood (it is squared inside
                    the computation)
            
        """

        # Covariance function hyperparameters
        self.l = l
        self.v0 = v0
        self.alpha = alpha
        
        self.m = m

        # Accepted values: squared_exponential, rational_quadratic, or callable
        if k is "squared_exponential":
            self.k = lambda x0, x1: k_squared_exp(x0, x1, self.v0, self.l)
        elif k is "rational_quadratic":
            self.k = lambda x0, x1: k_rational_quadratic(x0, x1, self.alpha, self.l)
        else:
            self.k = k

        self.sigma = sigma

        # Used for storing fitted results when sampling from posterior
        self._M = None
        self._K = None

        # This is the updated mean and variance function.
        self.m_hat = self.m
        self.k_hat = self.k

    def fit(self, X, Y):
        """
           Fits the posterior m_hat, k_hat functions. Discards all data seen so far.

           Parameters
           ----------

           X : a vector of points at which measurements were taken

           Y : a vector of measurements

           Returns
           -------

           self : with the re-computed functions m_hat and k_hat
        """

        # If we are using fit and not partial_fit, re-write posterior
        self.m_hat = self.m
        self.k_hat = self.k

        self.Y = np.array([])
        self.X = np.array([])

        return self.partial_fit(X,Y)

    def partial_fit(self, X, Y):
        """
           Fits the posterior m_hat, k_hat functions.

           Parameters
           ----------

           X : a vector of points at which measurements were taken

           Y : a vector of measurements

           Returns
           -------

           self : with the re-computed functions m_hat and k_hat
        """

        # To compute:
        #
        #  self.M_x = m(X)
        #  self.delta_M_x = Y - m(X)
        #  self.K_xx = (K(X,X) + I*(sigma**2))**-1
        #  self.K_xx_inverse = self.K_xx**-1


        self.Y = np.append(self.Y, Y)
        self.X = np.append(self.X, X)

        self.M_x = np.zeros(len(X))
        self.delta_M_x = np.zeros(len(X)) # stores: y - m(x)

        for i, x_i in enumerate(X):
            self.M_x[i] = self.m(x_i)
            self.delta_M_x[i] = Y[i] - self.M_x[i]

        self.K_XX = np.zeros((len(X),len(X)))

        for i, x_i in enumerate(X):
            self.K_XX[i,i] += self.sigma**2
            for j, x_j in enumerate(X):
                self.K_XX[i,j] += self.k(x_i, x_j)

        self.K_XX_inverse = np.linalg.inv(self.K_XX)

        self.m_hat = self._update_m()
        self.k_hat = self._update_k()

        return self

    def _update_m(self):
        """Updates m to m_hat according to the posterior equations."""
        m_hat_1 = self.m
        K_xX = lambda x: np.array([self.k(x,x_i) for x_i in self.X])       
        m_hat_2 = lambda x: np.dot(np.dot(K_xX(x), self.K_XX_inverse), self.delta_M_x.T)

        return lambda x: m_hat_1(x) + m_hat_2(x)

    def _update_k(self):
        """Updates k to k_hat according to the posterior equations."""
        K_xX = lambda x: np.array([self.k(x,x_i) for x_i in self.X])

        return lambda x_1, x_2: self.k(x_1, x_2) - np.dot(np.dot(K_xX(x_1),self.K_XX_inverse),K_xX(x_2).T)

    def draw_posterior(self, x):
        """
           Parameters
           ----------

           x : a vector of points at which to sample

           Returns
           -------

           y : a vector of sampled values
        """
        M = np.zeros((len(x)))

        for i, x_i in enumerate(x):
            M[i] = self.m_hat(x_i)

        K = np.zeros((len(x),len(x)))

        for i, x_i in enumerate(x):
            for j, x_j in enumerate(x):
                K[i,j] += self.k_hat(x_i,x_j)

        y = np.random.multivariate_normal(M,K)

        return y

    def draw_prior(self, x):
        """
           Parameters
           ----------

           x : a vector of points at which to sample.

           Returns
           -------

           y : a vector of sampled values
        """
        M = np.zeros((len(x)))

        for i, x_i in enumerate(x):
            M[i] = self.m(x_i)

        K = np.zeros((len(x),len(x)))

        for i, x_i in enumerate(x):
            for j, x_j in enumerate(x):
                K[i,j] = self.k(x_i, x_j)

        y = np.random.multivariate_normal(M,K)

        return y

##################################################
def sample_and_plot_prior(GP, num_samples, x, fig_name = 'prior_sample'):
    fig = plt.figure(fig_name)
    ax = fig.add_subplot('111')
    ax.set_title('Prior samples from GP')
    for i in xrange(num_samples):
        y = GP.draw_prior(x)
        ax.plot(x,y)
    plt.xlim(min(x),max(x))

    return fig

def sample_and_plot_posterior(GP, X, Y, num_samples, x):
    GP.fit(X,Y)
    
    fig = plt.figure()
    ax = fig.add_subplot('111')

    ax.plot(X,Y,'o')

    max_x = max(max(X),max(x))
    min_x = min(min(X),min(x))
    max_y = max(Y)
    min_y = min(Y)

    for i in xrange(num_samples):
        y = GP.draw_posterior(x)
        
        max_y = max(max_y,max(y))
        min_y = min(min_y,min(y))

        ax.plot(x,y)

    plt.xlim(min_x,max_x)
    plt.ylim(min_y - 0.5, max_y + 0.5)

    return fig

def posterior_with_length(l, sampling_grid, num_samples, X, Y):
    fig = plt.figure()

    ax1 = fig.add_subplot('221')
    ax1.set_title('SE Prior samples, l = ' + str(l))
    ax1.plot(X, Y, 'o')
    GP1 = GaussianProcess(lambda x: 0, k = 'squared_exponential', l = l)
    for i in xrange(num_samples):
        prior_sample = GP1.draw_prior(sampling_grid)
        ax1.plot(x, prior_sample)

    ax2 = fig.add_subplot('222')
    ax2.set_title('SE Posterior samples, l = ' + str(l))
    ax2.plot(X, Y, 'o')
    GP1.fit(X,Y)
    for i in xrange(num_samples):
        posterior_sample = GP1.draw_posterior(sampling_grid)
        ax2.plot(x, posterior_sample)

    ax3 = fig.add_subplot('223')
    ax3.set_title('RQ Prior samples, l = ' + str(l))
    ax3.plot(X, Y, 'o')
    GP2 = GaussianProcess(lambda x: 0, k = 'squared_exponential', l = l)
    for i in xrange(num_samples):
        prior_sample = GP2.draw_prior(sampling_grid)
        ax3.plot(x, prior_sample)

    ax4 = fig.add_subplot('224')
    ax4.set_title('RQ Posterior samples, l = ' + str(l))
    ax4.plot(X, Y, 'o')
    GP2.fit(X,Y)
    for i in xrange(num_samples):
        posterior_sample = GP2.draw_posterior(sampling_grid)
        ax4.plot(x, posterior_sample)
        
    plt.show()

def squared_exponential_pointwise_variance(v0, sampling_grid, num_samples, X, Y):
    fig = plt.figure()

    ax1 = fig.add_subplot('211')
    ax1.set_title('SE Prior samples, v0 = ' + str(v0))
    ax1.plot(X, Y, 'o')
    GP1 = GaussianProcess(lambda x: 0, k = 'squared_exponential', l = 0.3, v0 = v0)
    for i in xrange(num_samples):
        prior_sample = GP1.draw_prior(sampling_grid)
        ax1.plot(x, prior_sample)

    ax2 = fig.add_subplot('212')
    ax2.set_title('SE Posterior samples, v0 = ' + str(v0))
    ax2.plot(X, Y, 'o')
    GP1.fit(X,Y)
    for i in xrange(num_samples):
        posterior_sample = GP1.draw_posterior(sampling_grid)
        ax2.plot(x, posterior_sample)

    plt.show()

if __name__ == '__main__':

    # The function we are 'recovering'
    sigma = 2.0
    l_range = np.array([0.01, 0.05, 0.1, 0.5, 1.0, 1.5, 5.0, 100.0])
    v0_range = np.array([0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0])
    #f = lambda x: x**5 + 0.7*(x**4) + 3*(x**3) - 0.9*x + 0.3

    f = lambda x: x**3 + 1.2*(x**2) - x + 0.3
    
    x = np.arange(-2.0, 2.0, 0.1)
    
    X = np.array([-1.8, -1.1, -1.0, -0.95 -0.9, -0.4, 1.7, 1.9])
    Y = np.array([np.random.normal(f(x_i), sigma**2) for x_i in X])

    GP_1 = GaussianProcess(lambda x: 0, "squared_exponential", l = 0.01, sigma = sigma)
    GP_2 = GaussianProcess(lambda x: 0, "rational_quadratic", l = 0.01, sigma = sigma)

    #    print "Sampling GP 1 prior..."
    #    plt.figure()
    #    for n in range(f_x_samples):
    #        y = GP_1.draw_prior(x)
    #        plt.plot(x, y)
    #    plt.show()

    plt.figure()
    plt.plot(x, [f(x_i) for x_i in x])
    plt.show()

    for l in l_range:
        print "Sampling prior and posterior for length scale", l
        posterior_with_length(l, x, 5, X, Y)

    for v0 in v0_range:
        print "Sampling prior and posterior for pointwise variance", v0
        squared_exponential_pointwise_variance(v0, x, 5, X, Y)

    #print "Sampling GP 1 posterior..."
    #fig = plt.figure()

    #ax1 = fig.add_subplot('211')
    #GP_1.fit(X,Y)
    #for n in range(f_x_samples):
    #    y = GP_1.draw_posterior(x)
    #    ax1.plot(x,y)

    #    print "Sampling GP 2 prior..."
    #    plt.figure()
    #    for n in range(f_x_samples):
    #        y = GP_2.draw_prior(x)
    #        plt.plot(x, y)
    #    plt.show()

    #print "Sampling GP 2 posterior..."
    #ax2 = fig.add_subplot('212')
    #GP_2.fit(X,Y)
    #for n in range(f_x_samples):
    #    y = GP_2.draw_posterior(x)
    #    ax2.plot(x,y)

    #plt.show()

