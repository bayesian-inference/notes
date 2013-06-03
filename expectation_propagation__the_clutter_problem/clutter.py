#!/usr/bin/env python
# vim: set fileencoding=utf-8
#
# This code implements the clutter problem and its solution using EP, as
# described in
#
# @inproceedings{Minka:2001:EPA:647235.720257,
#  author = {Minka, Thomas P.},
#  title = {Expectation Propagation for approximate Bayesian inference},
#  booktitle = {Proceedings of the 17th Conference in Uncertainty in Artificial
#               Intelligence},
#  series = {UAI '01},
#  year = {2001},
#  isbn = {1-55860-800-1},
#  pages = {362--369},
#  numpages = {8},
#  url = {http://dl.acm.org/citation.cfm?id=647235.720257},
#  acmid = {720257},
#  publisher = {Morgan Kaufmann Publishers Inc.},
#  address = {San Francisco, CA, USA},
# }
#
# Written in partial fulfillment of the subject requirements of MFF UK's
# NPFL108 2012/13.
#
# author: Matěj Korvas

__author__ = u"Matěj Korvas"
__date__ = u"2013-05"

import copy
from itertools import izip
from math import sqrt
import numpy as np
from scipy.stats import bernoulli, norm


DEBUG = True
# DEBUG = False
INTERACTIVE = True
if DEBUG:
    from alex.utils import pdbonerror
if INTERACTIVE:
    import matplotlib.pyplot as plt
    plt.ion()
was_interactive = INTERACTIVE


DEFAULT_CFG = {
    'n_dimensions': 1,
    'n_observations': 41,
    'clutter_ratio': 0.5,   # the higher, the more clutter
    'clutter_mean': 0.,
    'clutter_var': 10.,     # variance-covariance matrix of clutter will be
                            # <this number> * I
    'prior_var': 100.,      # variance-covariance matrix of x will be
                            # <this number> * I
    'max_n_iterations': 20,    # maximum number of iterations
    'tol': 10 ** -4,        # tolerance in change of parameters to determine
                            # when convergence has been reached
}


config = DEFAULT_CFG  # To be overwritten.


def psphere_norm(point, mu, variances):
    """Computes the multivariate normal PDF.  Multiplication could cause
    underflow, so using log-scale might be needed."""
    return np.prod(map(lambda pt_uni, mu_uni, sd: norm.pdf(pt_uni, mu_uni, sd),
                       point, mu, np.sqrt(variances)))


def gen_data():
    """Generates a sample of data conformant to options specified in the
    config."""
    # Introduce short variable names.
    N = config['n_observations']
    d = config['n_dimensions']
    cvar = config['clutter_var']
    w = config['clutter_ratio']

    # Define relevant distributions.
    ris_clutter = bernoulli(w).rvs
    clutter_uni = norm(0, sqrt(cvar))
    rclutter = lambda: [clutter_uni.rvs() for _ in xrange(d)]
    prior = norm(0., sqrt(config['prior_var']))

    # Sample from the distributions.
    x = np.array([prior.rvs() for _ in xrange(d)])
    x_rvss = [norm(x_i).rvs for x_i in x]
    rx = lambda: [xrvs() for xrvs in x_rvss]
    Y = np.vstack([(rclutter() if ris_clutter() else rx())
                   for _ in xrange(N)])

    # Return.
    return x, Y


def param_dist(new, old):
    """Measures distance between two vectors of parameters."""
    # if    any(np.isinf(num) for num in new.flat) or (
          # any(np.isinf(num) for num in old.flat)):
        # return np.infty
    # return np.max(np.abs(new - old))
    new_finfo = np.finfo(new.dtype)
    new_clipped = np.clip(new, new_finfo.min, new_finfo.max)
    old_finfo = np.finfo(old.dtype)
    old_clipped = np.clip(old, old_finfo.min, old_finfo.max)
    return max(np.abs(new_clipped - old_clipped))


def true_distro(true_x):
    w = config['clutter_ratio']
    d = config['n_dimensions']
    d_zeros = np.zeros(d)
    d_cvars = np.repeat(config['clutter_var'], d)
    pclutter = lambda y: psphere_norm([y], d_zeros, d_cvars)

    prob_y = lambda y: ((1 - w) * psphere_norm([y], true_x, np.repeat(1., d)) +
                        w * pclutter(y))

    return prob_y


def update_plot(Y, m_x, v_x, cavity_mx, cavity_vx, ms, vs, idxs_used,
                idxs_skipped, x_true, pclutter):
    global INTERACTIVE
    w = config['clutter_ratio']
    min_y = np.min(Y, axis=0)[0]
    max_y = np.max(Y, axis=0)[0]
    plt_ys = np.linspace(max_y + 1.1 * (min_y - max_y),
                         min_y + 1.5 * (max_y - min_y), 201)

    plt.cla()
    # Draw the posterior.
    post_pdf = norm(m_x, sqrt(v_x)).pdf
    plt.plot(plt_ys, map(post_pdf, plt_ys), label="x posterior")
    cluttered_x_post = [(1 - w) * post_pdf(y) + w * pclutter([y])
                        for y in plt_ys]  # sorry for confusing
                                            # xs and ys...
    plt.plot(plt_ys, cluttered_x_post, label="cluttered x posterior")
    # Draw the true posterior.
    plt.plot(plt_ys, map(true_distro(x_true), plt_ys),
                label="original distribution")
    # Draw the points used.
    plt.plot((x_true, x_true), (0, 0.5), label="true x")
    plt.scatter(Y[idxs_used], np.zeros(len(idxs_used)), color='b',
                label="used points")
    plt.scatter(Y[idxs_skipped], np.zeros(len(idxs_skipped)), color='r',
                label="skipped points")
    # Draw the cavity distribution.
    plt.plot(plt_ys, map(norm(cavity_mx, sqrt(cavity_vx)).pdf, plt_ys),
             label="cavity")
    # Draw the last point's distribution.
    last_y = Y[idxs_used[-1]][0]
    plt.scatter([last_y], [0.], color='g', label="last point")
    last_post_pdf = norm(ms[idxs_used[-1]], sqrt(vs[idxs_used[-1]])).pdf
    plt.plot(plt_ys, map(last_post_pdf, plt_ys),
             label="last point's posterior")
    plt.legend()
    plt.show()
    try:
        raw_input("Press enter...")
    except EOFError:
        INTERACTIVE = False


def plot_factors(Y, x_true, m_0, v_0, ms, vs):
    min_y = np.min(Y, axis=0)[0]
    max_y = np.max(Y, axis=0)[0]
    plt_ys = np.linspace(max_y + 1.1 * (min_y - max_y),
                         min_y + 1.5 * (max_y - min_y), 201)

    plt.cla()
    plt.plot((x_true, x_true), (0, 0.5), label="true x")
    plt.scatter(Y, np.zeros(Y.shape[0]))
    prior = norm(m_0, sqrt(v_0)).pdf
    plt.plot(plt_ys, map(prior, plt_ys), label="prior")
    # import ipdb; ipdb.set_trace()
    for obs_idx, (mi, vi) in enumerate(izip(ms, vs)):
        y_post = norm(mi, sqrt(vi)).pdf
        plt.plot(plt_ys, map(y_post, plt_ys))
    plt.legend()
    plt.show()
    try:
        raw_input("Press enter...")
    except EOFError:
        INTERACTIVE = False


def main():
    # Introduce short variable names.
    d = config['n_dimensions']
    N = config['n_observations']
    w = config['clutter_ratio']
    cvar = config['clutter_var']
    tol = config['tol']

    x_true, Y = gen_data()
    print "Data generated:"
    print "  x  = {x}".format(x=x_true)
    print "  Y  = {y}".format(y=Y)

    # Initialise the prior (not to change it ever more).
    m_0 = np.zeros(d)  # The prior of x shall be non-informative.
    v_0 = config['prior_var']
    s_0 = (2 * np.pi * v_0) ** (-.5 * d)
    # Initialise data terms.
    vs = np.repeat(np.infty, N)
    ms = np.zeros((N, d))
    ss = np.ones(N)
    old_vs = copy.copy(vs)
    old_ms = copy.copy(ms)
    old_ss = copy.copy(ss)
    # Initialise the estimated distribution for x to the prior.
    m_x = m_0
    v_x = v_0

    # Some precomputations.
    d_zeros = np.zeros(d)
    d_cvars = np.repeat(cvar, d)
    pclutter = lambda y: psphere_norm(y, d_zeros, d_cvars)

    # Iterate EP.
    # if DEBUG:
        # import ipdb; ipdb.set_trace()
    for iteration in xrange(config['max_n_iterations']):
        iteration_h = iteration + 1  # human-readable
        print "Iteration {i}:".format(i=iteration_h)

        idxs_used = list()
        idxs_skipped = list()

        for obs_idx in xrange(N):
            # Find the cavity distribution parameters.
            # (I.e., compute $\frac{q(x)}{\tilde{t}_i(x)}$ where $i$ is our
            # `obs_idx'.)
            # These are the well-known formulas for fraction of Gaussians.
            inv_vi = 1. / vs[obs_idx]
            inv_vx = 1. / v_x
            # Handle division by zero.
            if inv_vi == inv_vx:
                cavity_vx = np.infty
            else:
                cavity_vx = 1. / (inv_vx - inv_vi)
            # Minka:
            # cavity_mx = m_x + cavity_vx * inv_vi * (m_x - ms[obs_idx])
            # José:
            cavity_mx = cavity_vx * (m_x * inv_vx - ms[obs_idx] * inv_vi)
            # ...probably equivalent (checked on a few concrete cases) but not
            # trivially reducible one to the other, at least I can't see it.

            # Recompute m_x, v_x, Z_i from cavity_mx, cavity_vx.
            yi = Y[obs_idx]
            Zi = ((1 - w) * psphere_norm(yi, cavity_mx,
                                         np.repeat(cavity_vx + 1., d)) +
                  w * pclutter(yi))
            r = 1. - w * pclutter(yi) / Zi
            # Minka:
            new_m_x = (cavity_mx +
                       r * cavity_vx / (1. + cavity_vx) * (yi - cavity_mx))
            # José:
            # new_m_x = (cavity_mx +
            #            r * cavity_vx / (1. + cavity_vx) * (yi - m_0))
            # ...this yields bad results.
            new_v_x = (cavity_vx - r * (cavity_vx ** 2) / (cavity_vx + 1.) +
                       r * (1 - r) * cavity_vx ** 2
                       * np.sum((yi - cavity_mx) ** 2)
                       / (d * (cavity_vx + 1) ** 2))

            # Update $\tilde{t}$.
            inv_cavity_vx = 1. / cavity_vx
            inv_new_vx = 1. / new_v_x
            # Handle division by zero.
            if inv_cavity_vx == inv_new_vx:
                ## This would perhaps be the most principled choice.
                vi = np.infty
                ## The other option below, now commented, leads to very low Z at
                ## the expense of bad approximation of the true x.
                # print "Skipping {n}".format(n=obs_idx)
                # continue
                ## Another option is setting vi to a very high number.
                # vi = 2 ** 512
            else:
                vi = 1. / (inv_new_vx - inv_cavity_vx)
            # Eschew negative variance.
            if vi + cavity_vx < 0:
                print "Skipping {n}".format(n=obs_idx)
                idxs_skipped.append(obs_idx)
                continue
            else:
                m_x = new_m_x
                v_x = new_v_x
                vs[obs_idx] = vi
                idxs_used.append(obs_idx)

            # Eschew evaluating 0 * ∞.
            if vi == np.infty and (
                    inv_cavity_vx == 0 or any(m_x == cavity_mx)):
                ms[obs_idx] = cavity_mx
            else:
                ms[obs_idx] = cavity_mx + (vi + cavity_vx) * inv_cavity_vx * (
                    m_x - cavity_mx)
            m_prob = psphere_norm(ms[obs_idx], cavity_mx,
                                  np.repeat(vi + cavity_vx, d))
            ss[obs_idx] = Zi / m_prob if m_prob > 0. else np.infty
            if INTERACTIVE:
                update_plot(Y, m_x, v_x, cavity_mx, cavity_vx, ms,
                            vs, idxs_used, idxs_skipped, x_true, pclutter)

        # Check for convergence, update old_anything.
        dist = max(param_dist(new, old) for (new, old) in
                   ((vs, old_vs), (ms, old_ms), (ss, old_ss)))
        print "m_x = {mx}; v_x = {vx}".format(mx=m_x, vx=v_x)
        print "Maximum distance from last parameter values: {d}".format(
            d=dist)
        if dist <= tol:
            print "Done."
            break
        old_vs = copy.copy(vs)
        old_ms = copy.copy(ms)
        old_ss = copy.copy(ss)

    # Compute the normalising constant.
    # Ignore points (skipped in the last iteration) with ss = ∞.
    valid_idxs = np.isfinite(ss)
    B = m_x.dot(m_x) / v_x - sum(m.dot(m) / v for (m, v) in
                                 izip(ms[valid_idxs], vs[valid_idxs]))
    Z = ((2 * np.pi * v_x) ** (.5 * d) * np.exp(.5 * B)
         * s_0 * np.prod(ss[valid_idxs]))

    print "Results:"
    print "  vs = {vs}".format(vs=vs)
    print "  ms = {ms}".format(ms=ms)
    print "  ss = {ss}".format(ss=ss)
    print " m_x = {mx}".format(mx=m_x)
    print " v_x = {vx}".format(vx=v_x)
    print "  B  = {B}".format(B=B)
    print "  Z  = {Z}".format(Z=Z)

    if was_interactive:
        update_plot(Y, m_x, v_x, cavity_mx, cavity_vx, ms, vs,
                    idxs_used, idxs_skipped, x_true, pclutter)
        plot_factors(Y, x_true, m_0, v_0, ms, vs)


if __name__ == "__main__":
    # TODO Allow for changing the configuration from CL.
    main()
