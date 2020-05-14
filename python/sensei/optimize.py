from . import unpaired
from . import paired

def find_minimum_M(max_M, N, m, s, rho=None, alpha=0.05, beta=0.1, design='unpaired', test_type="one-sided", offset=0,
                   sign=0):
    """
    Find the minimum M satisfies a specific false negative rate (beta)
    :param max_M:
    :param N:
    :param a:
    :param b:
    :param rho:
    :param alpha:
    :param beta:
    :param design: study design: "unpaired" or "paired"
    :param test_type:
    :param offset:
    :param sign:
    :return: Minimum M
    """

    if design == "unpaired":
        test_M = lambda M: unpaired.calc_fn_rate(M, N, m, s, alpha, test_type, offset, sign)
    elif design == "paired":
        test_M = lambda M: paired.calc_fn_rate(M, N, m, s, rho, alpha, test_type, offset, sign)
    else:
        raise ValueError("Unknown study design")

    M = 3
    while test_M(M) > beta:
        if M > max_M:
            return -1
        else:
            M += 1

    return M