from prime_generator import generate_N
import math

def is_perfect_square(number):
    return round(number**0.5) * round(number**0.5) == number

def factorization_ferma(N):
    k = 1
    sqrt_N = math.floor(N**0.5)

    y = None
    x = None
    while True:
        if is_perfect_square((sqrt_N + k) ** 2 - N):
             y = ((sqrt_N + k) ** 2 - N) ** 0.5
             x = sqrt_N + k
             break
        k = k + 1
    p, q = int((x + y)), int((x - y))
    return p,q

if __name__ == "__main__":
    bits = 60
    N = generate_N(bits)
    p, q = factorization_ferma(658661)
    print(f"N = {N}")
    print()
    print(f"p, q = {p, q}")
    print(f"p * q = {p * q}")
    print((p* q)==N)