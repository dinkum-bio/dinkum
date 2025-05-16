_debug = False


def set_debug(state):
    global _debug
    _debug = bool(state)


def debug(s, *args, **kwargs):
    if args or kwargs:
        print(s, args, kwargs)
    else:
        print(s)
