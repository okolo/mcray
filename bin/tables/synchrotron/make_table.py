import numpy as np
from sys import stderr
from scipy.special import kv
from scipy.integrate import quad

def main():
# printing dimensionless synchrotron power spectrum
# see e.g. formula	(5.66) from https://www.cv.nrao.edu/~sransom/web/Ch5.html
    def integrand(x):
        return kv(5/3, x)
    max_rel_error = 1e-6
    log_x_min = -20
    log_x_max = 2
    n_bins_per_decade = 100
    x = np.logspace(log_x_min, log_x_max, n_bins_per_decade * (log_x_max-log_x_min) + 1)
    f = np.zeros_like(x)
    result, error = quad(integrand, x[-1], np.inf, limit=200, epsabs=0)
    if result > 0:
        assert error / result < max_rel_error, f'insufficient accuracy {error/result} for x = {x[-1]}'
    f[-1] = result

    for i in range(len(f)-1, 0, -1):
        result, error = quad(integrand, x[i-1], x[i], limit=1000, epsabs=0)
        if result > 0:
            assert error/result < max_rel_error, f'insufficient accuracy {error/result} for x = {x[i-1]}'
        f[i-1] = result + f[i]

    f *= x

    # check that \Int_0^{x_max} f(x)dx = 1
    checksum = np.sum(f*x) * np.log(10) / n_bins_per_decade
    print('checksum =', checksum, file=stderr)
    checksum = np.sum(f[1:]*(x[1:]-x[:-1]))
    print('checksum2 =', checksum, file=stderr)

    X = np.vstack((x, f)).transpose()
    np.savetxt('sync_spec.dat', X)

if __name__ == '__main__':
    main()

