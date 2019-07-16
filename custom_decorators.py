import time


def measure_time(func):

    def wrapper(*args, **kwargs):
        t1 = time.time()
        x = func()
        t2 = time.time()
        print(f"Time for {func.__name__} is {t2-t1}")
        return x
    return wrapper

