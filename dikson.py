import math
import random
from collections import defaultdict

def get_B_list(B):
    """Генерация списка простых чисел до B с использованием решета Эратосфена."""
    sieve = [True] * (B + 1)
    sieve[0:2] = [False, False]
    for p in range(2, int(math.sqrt(B)) + 1):
        if sieve[p]:
            sieve[p * p:B + 1:p] = [False] * len(range(p * p, B + 1, p))
    return [p for p, is_prime in enumerate(sieve) if is_prime]

def factorize_gladki(a, b_list):
    """Факторизация числа a по списку простых чисел b_list."""
    factors = defaultdict(int)
    for p in b_list:
        while a % p == 0:
            factors[p] += 1
            a //= p
    return factors if a == 1 else None

def gcd(a, b):
    """Вычисление НОД двух чисел."""
    while b:
        a, b = b, a % b
    return a

def factor_dikson(N):
    """Алгоритм факторизации Диксона с адаптивным увеличением параметров."""
    M, B = 1.0, int(math.exp(0.5 * math.sqrt(math.log(N) * math.log(math.log(N)))))
    
    while True:
        b_list = get_B_list(B)
        relations, matrix = [], []
        
        while len(relations) < len(b_list) + 1:
            b = random.randint(2, N - 1)
            a = pow(b, 2, N)
            factors = factorize_gladki(a, b_list)
            if factors:
                relations.append((b, factors))
                matrix.append([factors.get(p, 0) % 2 for p in b_list])
        
        for i in range(len(matrix)):
            for j in range(i + 1, len(matrix)):
                if matrix[i] == matrix[j]:
                    b1, factors1 = relations[i]
                    b2, factors2 = relations[j]
                    X = (b1 * b2) % N
                    Y = 1
                    for p in set(factors1) | set(factors2):
                        Y = (Y * pow(p, (factors1.get(p, 0) + factors2.get(p, 0)) // 2, N)) % N
                    d = gcd(abs(X - Y), N)
                    if 1 < d < N:
                        return d
        
        # Если не найден делитель, увеличиваем параметры
        M *= 1.1
        B = int(B * M)

if __name__ == "__main__":
    N = 428668018913
    divisor = factor_dikson(N)
    print(f"Найден делитель: {divisor}" if divisor else "Делитель не найден.")