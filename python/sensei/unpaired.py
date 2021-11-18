import scipy.stats
from .utils import *
from scipy.stats import mannwhitneyu, ttest_ind, betabinom

def calc_wilcoxon_fn(M, N, m, s, alpha = 0.05, n_sim = 10_000):
    """
    
    :param M: number of patients, as a list
    :param N: number of cells, as a list
    :param m: mean for both groups, as a list
    :param s: std for both groups, as a list
    :param alpha: significance level
    :param n_sim: simulation iterations
    :return: false negative rate, i.e., 1 - power
    """
    N0, N1 = N
    M0, M1 = M
    m0, m1 = m
    s0, s1 = s
    
    a0, b0 = normal_to_beta(m0, s0)
    r0 = betabinom.rvs(N0, a0, b0, size=(M0, n_sim)) / n

    a1, b1 = normal_to_beta(m1, s1)
    r1 = betabinom.rvs(N1, a1, b1, size=(M1, n_sim)) / n
    
    return 1 - sum(mannwhitneyu(r0, r1).pvalue < alpha) / n_sim


def calc_fn_rate_beta(M, N, a, b, alpha=0.05, test_type="one-sided", offset=0, sign=0):
    """
    Calculate false negative rate
    :param M: number of patients
    :param N: number of cells
    :param a: Beta(a, b)
    :param b: Beta(a, b)
    :param alpha: significance level
    :param test_type: one-sided or two-sided
    :param offset:
    :param sign:
    :return: false negative rate, i.e., 1 - power
    """

    if not is_iterable(M):
        M = [M, M]

    if not is_iterable(N):
        N = [N, N]

    Ep = [a[0] / (a[0] + b[0]),
          a[1] / (a[1] + b[1])]

    # Vp = [a[0] * b[0] * (a[0] + b[0] + N[0]) / ((N[0] * (a[0] + b[0]) * (a[0] + b[0])) * (a[0] + b[0] + 1)),
    #      a[1] * b[1] * (a[1] + b[1] + N[1]) / ((N[1] * (a[1] + b[1]) * (a[1] + b[1])) * (a[1] + b[1] + 1))]

    Vp = [var_betabinom_over_n(N[0], a[0], b[0]), var_betabinom_over_n(N[1], a[1], b[1])]

    Et = (Ep[1] - (Ep[0] + offset)) / (Vp[0] / M[0] + Vp[1] / M[1]) ** .5

    if sign == 0:
        Et = abs(Et)
    else:
        Et = sign * Et

    nu = (Vp[0] / M[0] + Vp[1] / M[1]) ** 2 / ((Vp[0] / M[0]) ** 2 / (M[0] - 1) + (Vp[1] / M[1]) ** 2 / (M[1] - 1))

    if test_type == "one-sided":
        t_star = scipy.stats.t.ppf(q=1 - alpha, df=nu)
    elif test_type == "two-sided":
        t_star = scipy.stats.t.ppf(q=1 - alpha / 2, df=nu)
    else:
        raise ValueError("test must be one-sided or two-sided")
    return scipy.stats.t.cdf(t_star - Et, df=nu)


def calc_fn_rate(M, N, m, s, alpha, test_type, offset, sign):
    """

    :param M:
    :param N:
    :param m:
    :param s:
    :param alpha:
    :param test_type:
    :param offset:
    :param sign:
    :return:
    """
    if not is_iterable(s):
        s = [s, s]

    a = [None, None]
    b = [None, None]

    try:
        a[0], b[0] = normal_to_beta(m[0], s[0])
        a[1], b[1] = normal_to_beta(m[1], s[1])
    except ZeroDivisionError:
        return float("nan")
    return calc_fn_rate_beta(M, N, a, b, alpha, test_type, offset, sign)


def calc_fn_rate_override(M, N, m, s, alpha, test_type, override_diff):
    """

    :param M:
    :param N:
    :param m:
    :param s:
    :param alpha:
    :param test_type:
    :param override_diff: overriden difference
    :return:
    """
    if not is_iterable(s):
        s = [s, s]

    a = [None, None]
    b = [None, None]

    try:
        a[0], b[0] = normal_to_beta(m[0], s[0])
        a[1], b[1] = normal_to_beta(m[1], s[1])
    except ZeroDivisionError:
        return float("nan")

    if not is_iterable(M):
        M = [M, M]

    if not is_iterable(N):
        N = [N, N]

    Ep = [a[0] / (a[0] + b[0]),
          a[1] / (a[1] + b[1])]

    # Vp = [a[0] * b[0] * (a[0] + b[0] + N[0]) / ((N[0] * (a[0] + b[0]) * (a[0] + b[0])) * (a[0] + b[0] + 1)),
    #       a[1] * b[1] * (a[1] + b[1] + N[1]) / ((N[1] * (a[1] + b[1]) * (a[1] + b[1])) * (a[1] + b[1] + 1))]

    Vp = [var_betabinom_over_n(N[0], a[0], b[0]), var_betabinom_over_n(N[1], a[1], b[1])]

    Et = override_diff / (Vp[0] / M[0] + Vp[1] / M[1]) ** .5

    Et = abs(Et)

    nu = (Vp[0] / M[0] + Vp[1] / M[1]) ** 2 / ((Vp[0] / M[0]) ** 2 / (M[0] - 1) + (Vp[1] / M[1]) ** 2 / (M[1] - 1))

    if test_type == "one-sided":
        t_star = scipy.stats.t.ppf(q=1 - alpha, df=nu)
    elif test_type == "two-sided":
        t_star = scipy.stats.t.ppf(q=1 - alpha / 2, df=nu)
    else:
        raise ValueError("test must be one-sided or two-sided")
    return scipy.stats.t.cdf(t_star - Et, df=nu)


def calc_fn_rate_baseline(M, m, s, alpha, test_type, offset, sign):
    """

    :param M:
    :param N:
    :param m:
    :param s:
    :param alpha:
    :param test_type:
    :param override_diff: overriden difference
    :return:
    """
    if not is_iterable(s):
        s = [s, s]

    a = [None, None]
    b = [None, None]

    try:
        a[0], b[0] = normal_to_beta(m[0], s[0])
        a[1], b[1] = normal_to_beta(m[1], s[1])
    except ZeroDivisionError:
        return float("nan")

    if not is_iterable(M):
        M = [M, M]

    Ep = m

    Vp = [s[0] ** 2, s[1] ** 2]

    Et = (Ep[1] - (Ep[0] + offset)) / (Vp[0] / M[0] + Vp[1] / M[1]) ** .5

    if sign == 0:
        Et = abs(Et)
    else:
        Et = sign * Et

    nu = (Vp[0] / M[0] + Vp[1] / M[1]) ** 2 / ((Vp[0] / M[0]) ** 2 / (M[0] - 1) + (Vp[1] / M[1]) ** 2 / (M[1] - 1))

    if test_type == "one-sided":
        t_star = scipy.stats.t.ppf(q=1 - alpha, df=nu)
    elif test_type == "two-sided":
        t_star = scipy.stats.t.ppf(q=1 - alpha / 2, df=nu)
    else:
        raise ValueError("test must be one-sided or two-sided")
    return scipy.stats.t.cdf(t_star - Et, df=nu)