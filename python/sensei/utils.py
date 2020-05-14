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


def f_bound(f, key, first, last, value, which, *args, **kwargs):
    """

    :param f: the MONOTONIC function
    :param key: keyword arguments or position of argument to be searched
    :param first: [first, last)
    :param last: [first, last)
    :param value: value to be searched
    :param which: lower bound or upper bound
    :param args: arguments for f
    :param kwargs: arguments for f
    :return: first key making f(key, *args, **kwargs) not smaller than value
    """
    args = list(args)

    if which == "lower":
        cmp = lambda x, y: x < y
    elif which == "upper":
        cmp = lambda x, y: x > y
    else:
        raise ValueError("which must be lower or upper")

    if isinstance(key, int):
        args.insert(key, None)
        while first < last:
            mid = first + (last - first) // 2
            args[key] = mid
            if cmp(f(*args, **kwargs), value):
                first = mid + 1
            else:
                last = mid
    elif isinstance(key, str):
        while first < last:
            mid = first + (last - first) // 2
            kwargs[key] = mid
            if cmp(f(*args, **kwargs), value):
                first = mid + 1
            else:
                last = mid
    else:
        raise TypeError("key must be int for a position argument, or str for a keyword argument")

    return first


def f_lower_bound(f, key, first, last, value, *args, **kwargs):
    return f_bound(f, key, first, last, value, "lower", *args, **kwargs)


def f_upper_bound(f, key, first, last, value, *args, **kwargs):
    return f_bound(f, key, first, last, value, "upper", *args, **kwargs)
