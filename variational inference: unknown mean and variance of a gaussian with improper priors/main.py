import numpy as np

from scipy.stats import gamma, norm


def gen_data(n, mu_low, mu_high, tau_low, tau_high):
    data = []
    print "-----------------------------------------"
    print "prior for mean:      [%.2f, %.2f]" % (mu_low, mu_high)
    print "prior for precision: [%.2f, %.2f]" % (tau_low, tau_high)
    print "-----------------------------------------"
    print "%-9s %-9s %-10s" % ("mean", "variance", "x ~ N(mean, variance)")
    print "-----------------------------------------"

    mu = np.random.uniform(mu_low, mu_high)
    tau = np.random.uniform(tau_low, tau_high)
    for i in range(n):
        data.append(np.random.normal(mu, 1/tau))
        print "%f %9f %9f" % (mu, 1/tau, data[-1])
    print "-----------------------------------------"
    return np.array(data)

def mean_field(x, tau_len):
    N = len(x)
    a = N / 2.
    l = 1/(1/12. * tau_len**2)
    S = np.sum((x - np.mean(x))**2)
    #compute_b = lambda l: (1./2) * ((N/l) + S)
    sx = np.sum(x)
    sx2 = np.sum(x**2)
    mx = np.mean(x)
    compute_b = lambda l: (1./2) * (sx2 - 2*sx*mx + N * ((1/l) + mx**2))
    b = compute_b(l)

    for i in range(10):
        l = N*a/b
        b = compute_b(l)

    return np.mean(x), l, a, b


print "Generated data:"
n=100000
mu_low=0
mu_high=0.01
tau_low=0.5
tau_high=5
x = gen_data(n, mu_low, mu_high, tau_low, tau_high)

mu, precision, a, b = mean_field(x, tau_high - tau_low)

print
print "Mean parameters:"
print "\tMean:      %9f" % mu
print "\tPrecision: %9f" % (precision)
print
print "Precision parameteres:"
print "\tShape:     %9f" % a
print "\tRate:      %9f" % b
print
print "Data: "
print "\tMean:      %9f" % np.mean(x)
print "\tPrecision: %9f" % (1/np.std(x))

