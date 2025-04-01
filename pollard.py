import math
import random

def gcd(a, b):
    """Вычисление НОД двух чисел."""
    while b:
        a, b = b, a % b
    return a

def pollard_rho(N):
    """Алгоритм факторизации Полларда (ρ-алгоритм) с оптимизациями."""
    if N % 2 == 0:
        return 2
    if N % 3 == 0:
        return 3
    if N % 5 == 0:
        return 5

    def f(x):
        return (pow(x, 2, N) + 1) % N

    while True:
        x = random.randint(2, N - 1)
        y = x
        d = 1
        while d == 1:
            x = f(x)
            y = f(f(y))
            d = gcd(abs(x - y), N)
        if d != N:
            return d

def is_prime(n):
    """Проверка, является ли число n простым (тест Миллера-Рабина)."""
    if n < 2:
        return False
    for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
        if n % p == 0:
            return n == p
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
        if a >= n:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

def factorize(N):
    """Рекурсивная функция для нахождения всех делителей числа N."""
    if N == 1:
        return []
    if is_prime(N):
        return [N]
    divisor = pollard_rho(N)
    return factorize(divisor) + factorize(N // divisor)

if __name__ == "__main__":
    N = int(880612586747802414069241)
    factors = factorize(N)
    print(factors)