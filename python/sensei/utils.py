
def is_iterable(x):
    try:
        iterator = iter(x)
    except TypeError:
        return False
    else:
        return True


def normal_to_beta(m, s):
    v = s ** 2
    a = ((1 - m) / v - 1 / m) * (m ** 2)
    b = a * (1 / m - 1)
    return a, b


def beta_to_normal(a, b):
    m = a / (a + b)
    v = a * b / (a + b) ** 2 / (a + b + 1)
    return m, v ** .5