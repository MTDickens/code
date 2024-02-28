from decimal import *
from math import e

# Set the precision to 5
getcontext().prec = 7

# Method 1: J_0 = 1 - 1 / e, J_n = 1 - n J_{n-1}
def naive_J(n):
    assert isinstance(n, int), "n should be an integer"
    assert n >= 0, "n should be non-negative"

    if n == 0:
        return Decimal(1) - Decimal(1) / Decimal(e)
    else:
        return Decimal(1) - Decimal(n) * naive_J(n - 1)
    
# Method 2: J_{n-1} = 1 / n * (1 - J_n), J_n 
max_val = 15

def clever_J_upper(n):
    assert isinstance(n, int), "n should be an integer"
    assert n >= 0, "n should be non-negative"    

    if n > max_val:
        return "n should be less than max_val"
    elif n == max_val:
        return Decimal(e) / Decimal(n+1)
    else:
        return (1 - clever_J_upper(n + 1)) / Decimal(n + 1);

def clever_J_lower(n):
    assert isinstance(n, int), "n should be an integer"
    assert n >= 0, "n should be non-negative"    

    if n > max_val:
        return "n should be less than max_val"
    elif n == max_val:
        return Decimal(1) / Decimal(n+1)
    else:
        return (1 - clever_J_upper(n + 1)) / Decimal(n + 1);

# print func
def print_result(f, count=15):
    assert callable(f), "f should be a function"
    assert isinstance(count, int), "i should be an integer"
    assert count >= 0, "count should be non-negative"

    for i in map(lambda x: f.__name__ + '(' + str(x) + '): ' + str(f(x)), list(range(count + 1))): print(i)