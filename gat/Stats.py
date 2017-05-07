'''
Stats.py
======

statistical functions
'''

import math
import numpy
import types
from functools import reduce

try:
    import scipy.interpolate
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


class FDRResult:

    def __init__(self):
        pass


def computeQValues(pvalues,
                   vlambda=None,
                   pi0_method="smoother",
                   fdr_level=None,
                   robust=False,
                   smooth_df=3,
                   smooth_log_pi0=False,
                   pi0=None):
    """compute qvalues after the method by Storey et al. (2002)

    The python code derives from the R implementation at 
    http://genomics.princeton.edu/storeylab/qvalue/linux.html.
    """

    if min(pvalues) < 0 or max(pvalues) > 1:
        raise ValueError("p-values out of range")

    m = len(pvalues)
    pvalues = numpy.array(pvalues, dtype=numpy.float)

    if vlambda == None:
        vlambda = numpy.arange(0, 0.95, 0.05)

    if pi0 == None:
        if type(vlambda) == float:
            vlambda = (vlambda,)

        if len(vlambda) > 1 and len(vlambda) < 4:
            raise ValueError(
                " if length of vlambda greater than 1, you need at least 4 values.")

        if len(vlambda) > 1 and (min(vlambda) < 0 or max(vlambda) >= 1):
            raise ValueError("vlambda must be within [0, 1).")

        # estimate pi0
        if len(vlambda) == 1:
            vlambda = vlambda[0]
            if vlambda < 0 or vlambda >= 1:
                raise ValueError("vlambda must be within [0, 1).")

            pi0 = numpy.mean([x >= vlambda for x in pvalues]) / (1.0 - vlambda)
            pi0 = min(pi0, 1.0)
        else:

            pi0 = numpy.zeros(len(vlambda), numpy.float)

            for i in range(len(vlambda)):
                pi0[i] = numpy.mean([x >= vlambda[i]
                                     for x in pvalues]) / (1.0 - vlambda[i])

            if pi0_method == "smoother":
                if smooth_log_pi0:
                    pi0 = numpy.log(pi0)
                if HAS_SCIPY:
                    tck = scipy.interpolate.splrep(
                        vlambda, pi0, k=smooth_df, s=10000)
                    pi0 = scipy.interpolate.splev(max(vlambda), tck)
                else:
                    raise ImportError("pi0_method smoother requires scipy")

                if smooth_log_pi0:
                    pi0 = numpy.exp(pi0)

            elif pi0_method == "bootstrap":
                minpi0 = min(pi0)

                mse = numpy.zeros(len(vlambda), numpy.float)
                pi0_boot = numpy.zeros(len(vlambda), numpy.float)

                for i in range(100):
                    # sample pvalues
                    idx_boot = numpy.random.random_integers(0, m - 1, m)
                    pvalues_boot = pvalues[idx_boot]

                    for x in range(len(vlambda)):
                        # compute number of pvalues larger than lambda[x]
                        pi0_boot[x] = numpy.mean(
                            pvalues_boot > vlambda[x]) / (1.0 - vlambda[x])
                    mse += (pi0_boot - minpi0) ** 2
                pi0 = min(pi0[mse == min(mse)])
            else:
                raise ValueError(
                    "'pi0_method' must be one of 'smoother' or 'bootstrap'.")

            pi0 = min(pi0, 1.0)

    if pi0 <= 0:
        raise ValueError(
            "The estimated pi0 <= 0 (%f). Check that you have valid p-values or use another vlambda method." % pi0)

    if fdr_level != None and (fdr_level <= 0 or fdr_level > 1):
        raise ValueError("'fdr_level' must be within (0, 1].")

    # compute qvalues

    idx = numpy.argsort(pvalues)
    # monotonically decreasing bins, so that bins[i-1] > x >=  bins[i]
    bins = numpy.unique(pvalues)[::-1]

    # v[i] = number of observations less than or equal to pvalue[i]
    # could this be done more elegantly?
    val2bin = len(bins) - numpy.digitize(pvalues, bins)
    v = numpy.zeros(m, dtype=numpy.int)
    lastbin = None
    for x in range(m - 1, -1, -1):
        bin = val2bin[idx[x]]
        if bin != lastbin:
            c = x
        v[idx[x]] = c + 1
        lastbin = bin

    qvalues = pvalues * pi0 * m / v
    if robust:
        qvalues /= (1.0 - (1.0 - pvalues) ** m)

    # bound qvalues by 1 and make them monotonic
    qvalues[idx[m - 1]] = min(qvalues[idx[m - 1]], 1.0)
    for i in range(m - 2, -1, -1):
        qvalues[idx[i]] = min(min(qvalues[idx[i]], qvalues[idx[i + 1]]), 1.0)

    # fill result
    result = FDRResult()
    result.qvalues = qvalues

    if fdr_level != None:
        result.passed = [x <= fdr_level for x in result.qvalues]
    else:
        result.passed = [False for x in result.qvalues]

    result.pvalues = pvalues
    result.pi0 = pi0
    result.vlambda = vlambda
    result.fdr_level = fdr_level

    return result


# {{{ http://code.activestate.com/recipes/511478/ (r1)
def percentile(N, percent, key=lambda x: x):
    """
    Find the percentile of a list of values.

    @parameter N - is a list of values. Note N MUST BE already sorted.
    @parameter percent - a float value from 0.0 to 1.0.
    @parameter key - optional key function to compute value from each element of N.

    @return - the percentile of the values
    """
    if len(N) == 0:
        return None
    k = (len(N) - 1) * percent
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return key(N[int(k)])
    d0 = key(N[int(f)]) * (c - k)
    d1 = key(N[int(c)]) * (k - f)
    return d0 + d1

###################################################################
###################################################################
###################################################################
# adjust P-Value
###################################################################


def adjustPValues(pvalues, method='fdr', n=None):
    '''returns an array of adjusted pvalues

    Reimplementation of p.adjust in the R package.

    p: numeric vector of p-values (possibly with 'NA's).  Any other
    R is coerced by 'as.numeric'.

    method: correction method. Valid values are:

    n: number of comparisons, must be at least 'length(p)'; only set
    this (to non-default) when you know what you are doing

    For more information, see the documentation of the
    p.adjust method in R.
    '''

    if n == None:
        n = len(pvalues)

    if method == "fdr":
        method = "BH"

    # optional, remove NA values
    p = numpy.array(pvalues, dtype=numpy.float)
    lp = len(p)

    assert n <= lp

    if n <= 1:
        return p
    if n == 2 and method == "hommel":
        method = "hochberg"

    if method == "bonferroni":
        p0 = n * p
    elif method == "holm":
        i = numpy.arange(lp)
        o = numpy.argsort(p)
        ro = numpy.argsort(o)
        m = numpy.maximum.accumulate((n - i) * p[o])
        p0 = m[ro]
    elif method == "hommel":
        raise NotImplementedError("hommel method not fully implemented")
    elif method == "hochberg":
        i = numpy.arange(0, lp)[::-1]
        o = numpy.argsort(1 - p)
        ro = numpy.argsort(o)
        m = numpy.minimum.accumulate((n - i) * p[o])
        p0 = m[ro]
    elif method == "BH":
        i = numpy.arange(1, lp + 1)[::-1]
        o = numpy.argsort(1 - p)
        ro = numpy.argsort(o)
        m = numpy.minimum.accumulate(float(n) / i * p[o])
        p0 = m[ro]
    elif method == "BY":
        i = numpy.arange(1, lp + 1)[::-1]
        o = numpy.argsort(1 - p)
        ro = numpy.argsort(o)
        q = numpy.sum(1.0 / numpy.arange(1, n + 1))
        m = numpy.minimum.accumulate(q * float(n) / i * p[o])
        p0 = m[ro]
    elif method == "none":
        p0 = p

    return numpy.minimum(p0, numpy.ones(len(p0)))


class Result(object):

    '''allow both member and dictionary access.'''
    slots = ("_data")

    def __init__(self):
        object.__setattr__(self, "_data", dict())

    def fromR(self, take, r_result):
        '''convert from an *r_result* dictionary using map *take*.

        *take* is a list of tuples mapping a field to the corresponding
        field in *r_result*.
        '''
        for x, y in take:
            if y:
                self._data[x] = r_result.rx(y)[0][0]
            else:
                self._data[x] = r_result.rx(x)[0][0]

            # if y:
            #     self._data[x] = r_result[y]
            # else:
            #     self._data[x] = r_result[x]

        return self

    def __getattr__(self, key):
        if not key.startswith("_"):
            try:
                return object.__getattribute__(self, "_data")[key]
            except KeyError:
                pass
        return getattr(self._data, key)

    def keys(self):
        return list(self._data.keys())

    def values(self):
        return list(self._data.values())

    def __len__(self):
        return self._data.__len__()

    def __str__(self):
        return str(self._data)

    def __contains__(self, key):
        return key in self._data

    def __getitem__(self, key):
        return self._data[key]

    def __delitem__(self, key):
        del self._data[key]

    def __setitem__(self, key, value):
        self._data[key] = value

    def __setattr__(self, key, value):
        if not key.startswith("_"):
            self._data[key] = value
        else:
            object.__setattr__(self, key, value)


class Summary(Result):

    """a collection of distributional parameters. Available properties
    are:

    mean, median, min, max, samplestd, sum, counts
    """

    fields = ("nval", "min", "max", "mean",
              "median", "stddev", "sum", "q1", "q3")

    def __init__(self, values=None,
                 format="%6.4f", mode="float",
                 allow_empty=True):

        Result.__init__(self)
        self._format = format
        self._mode = mode

        # note that this determintes the order of the fields at output
        self.counts, self.min, self.max, self.mean, self.median, self.samplestd, self.sum, self.q1, self.q3 = \
            (0, 0, 0, 0, 0, 0, 0, 0, 0)

        if values != None:

            values = [x for x in values if x != None]

            if len(values) == 0:
                if allow_empty:
                    return
                else:
                    raise ValueError("no data for statistics")

            # convert
            self._nerrors = 0
            if type(values[0]) not in (int, float):
                n = []
                for x in values:
                    try:
                        n.append(float(x))
                    except ValueError:
                        self._nerrors += 1
            else:
                n = values

            # use a non-sort algorithm?
            n.sort()
            if len(n):
                self.q1 = n[len(n) / 4]
                self.q3 = n[len(n) * 3 / 4]
            else:
                self.q1 = self.q3 = 0

            self.counts = len(n)
            self.min = min(n)
            self.max = max(n)
            self.mean = numpy.mean(n)
            self.median = numpy.median(n)
            self.samplestd = numpy.std(n)
            self.sum = reduce(lambda x, y: x + y, n)

    def getHeaders(self):
        """returns header of column separated values."""
        return self.fields

    def getHeader(self):
        """returns header of column separated values."""
        return "\t".join(self.getHeaders())

    def __str__(self):
        """return string representation of data."""

        if self._mode == "int":
            format_vals = "%i"
            format_median = "%.1f"
        else:
            format_vals = self._format
            format_median = self._format

        return "\t".join(("%i" % self.counts,
                          format_vals % self.min,
                          format_vals % self.max,
                          self._format % self.mean,
                          format_median % self.median,
                          self._format % self.samplestd,
                          format_vals % self.sum,
                          format_vals % self.q1,
                          format_vals % self.q3,
                          ))
