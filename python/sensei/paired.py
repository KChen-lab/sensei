import scipy.stats
from .utils import *


def calc_fn_rate_beta(M, N, a, b, rho, alpha=0.05, test_type="one-sided", offset=0, sign=0):
    """

    :param M:
    :param N:
    :param a:
    :param b:
    :param rho:
    :param alpha:
    :param test_type:
    :param offset:
    :param sign:
    :return:
    """
    if is_iterable(M):
        if M[0] == M[1]:
            M = M[0]
        else:
            raise ValueError("Two groups must have the same samples size for paired test.")

    if not is_iterable(N):
        N = [N, N]


    cov = rho * (a[0] * b[0] / ((a[0] + b[0]) * (a[0] + b[0]) * (a[0] + b[0] + 1))) ** .5 * \
                (a[1] * b[1] / ((a[1] + b[1]) * (a[1] + b[1]) * (a[1] + b[1] + 1))) ** .5

    Ep = [a[0] / (a[0] + b[0]), a[1] / (a[1] + b[1])]

    Vp = [var_betabinom_over_n(N[0], a[0], b[0]), var_betabinom_over_n(N[1], a[1], b[1])]

    Et = (Ep[1] - (Ep[0] + offset)) * M ** .5 / (Vp[0] + Vp[1] - 2 * cov) ** .5

    if sign == 0:
        Et = abs(Et)
    else:
        Et = sign * Et

    nu = M - 1

    if test_type == "one-sided":
        t_star = scipy.stats.t.ppf(q=1 - alpha, df=nu)
    elif test_type == "two-sided":
        t_star = scipy.stats.t.ppf(q=1 - alpha / 2, df=nu)
    else:
        raise ValueError("test must be one-sided or two-sided")
    return scipy.stats.t.cdf(t_star - Et, df=nu)


def calc_fn_rate(M, N, m, s, rho, alpha, test_type, offset, sign):
    """

    :param M:
    :param N:
    :param m:
    :param s:
    :param rho:
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

    return calc_fn_rate_beta(M, N, a, b, rho, alpha, test_type, offset, sign)


def calc_fn_rate_override(M, N, m, s, rho, alpha, test_type, override_diff):
    """

    :param M:
    :param N:
    :param m:
    :param s:
    :param rho:
    :param alpha:
    :param test_type:
    :param offset:
    :param sign:
    :param override_diff
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

    if is_iterable(M):
        if M[0] == M[1]:
            M = M[0]
        else:
            raise ValueError("Two groups must have the same samples size for paired test.")

    if not is_iterable(N):
        N = [N, N]

    cov = rho * (a[0] * b[0] / ((a[0] + b[0]) * (a[0] + b[0]) * (a[0] + b[0] + 1))) ** .5 * \
          (a[1] * b[1] / ((a[1] + b[1]) * (a[1] + b[1]) * (a[1] + b[1] + 1))) ** .5

    # Ep = [a[0] / (a[0] + b[0]), a[1] / (a[1] + b[1])]

    Vp = [var_betabinom_over_n(N[0], a[0], b[0]), var_betabinom_over_n(N[1], a[1], b[1])]

    Et = override_diff * M ** .5 / (Vp[0] + Vp[1] - 2 * cov) ** .5

    Et = abs(Et)

    nu = M - 1

    if test_type == "one-sided":
        t_star = scipy.stats.t.ppf(q=1 - alpha, df=nu)
    elif test_type == "two-sided":
        t_star = scipy.stats.t.ppf(q=1 - alpha / 2, df=nu)
    else:
        raise ValueError("test must be one-sided or two-sided")
    return scipy.stats.t.cdf(t_star - Et, df=nu)


def calc_fn_rate_baseline(M, m, s, rho, alpha, test_type, offset, sign):
    """

    :param M:
    :param N:
    :param m:
    :param s:
    :param rho:
    :param alpha:
    :param test_type:
    :param offset:
    :param sign:
    :return:
    """
    if not is_iterable(s):
        s = [s, s]


    if is_iterable(M):
        if M[0] == M[1]:
            M = M[0]
        else:
            raise ValueError("Two groups must have the same samples size for paired test.")

    cov = rho * s[0] * s[1]

    Ep = m

    Vp = [s[0] ** 2, s[1] ** 2]

    Et = (Ep[1] - (Ep[0] + offset)) * M ** .5 / (Vp[0] + Vp[1] - 2 * cov) ** .5

    if sign == 0:
        Et = abs(Et)
    else:
        Et = sign * Et

    nu = M - 1

    if test_type == "one-sided":
        t_star = scipy.stats.t.ppf(q=1 - alpha, df=nu)
    elif test_type == "two-sided":
        t_star = scipy.stats.t.ppf(q=1 - alpha / 2, df=nu)
    else:
        raise ValueError("test must be one-sided or two-sided")
    return scipy.stats.t.cdf(t_star - Et, df=nu)