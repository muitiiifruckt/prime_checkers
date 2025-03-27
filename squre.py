import math
import numpy as np
from math import log, isqrt, gcd
from prime_generator import generate_N
import time


def legendre_symbol(a, p):
    """Вычисляет символ Лежандра (a/p)."""
    return pow(a, (p - 1) // 2, p)

def generate_factor_base(n, B):
    """Формирует факторную базу: простые p ≤ B, для которых (n/p)=1."""
    sieve = [True] * (B+1)
    sieve[0] = sieve[1] = False
    primes = []
    for p in range(2, B+1):
        if sieve[p]:
            if legendre_symbol(n, p) == 1:
                primes.append(p)
            for i in range(p*p, B+1, p):
                sieve[i] = False
    return primes

def quadratic_polynomial(A, x, n):
    """Возвращает значение Q(x) = (A+x)^2 - n."""
    return (A + x)**2 - n

def sieve_interval(n, A, M, factor_base):
    """Сегментированное решето для интервала x в [-M, M]."""
    interval = range(-M, M+1)
    Q_values = [quadratic_polynomial(A, x, n) for x in interval]
    log_Q = [math.log(abs(q)) if q != 0 else 0 for q in Q_values]
    # Копия логарифмических значений для корректировки
    sieve_array = log_Q.copy()

    # Для каждого простого из базы находим все решения уравнения: (A+x)^2 = n (mod p)
    smooth_candidates = {}
    for p in factor_base:
        # Решаем t^2 ≡ n (mod p). Находим решения методом перебора (p небольшие)
        sols = []
        for t in range(p):
            if (t*t - n) % p == 0:
                sols.append(t)
        # Для каждого решение переводим его к x: A+x ≡ t (mod p)
        for t in sols:
            # Находим минимальное x0 в интервале, удовлетворяющее A+x ≡ t (mod p)
            x0 = (t - A) % p
            # Перебор x = x0 + k*p, попадающих в интервал [-M, M]
            for k in range(-((M - x0) // p) - 1, ((M - x0) // p) + 2):
                x = x0 + k * p
                if -M <= x <= M:
                    idx = x + M  # сдвиг, чтобы индекс начинался с 0
                    sieve_array[idx] -= log(p)
    
    # Выбираем кандидатов, где остаточное значение близко к 0 (с некоторой поправкой)
    threshold = 0.5
    smooth_candidates = []
    for i, val in enumerate(sieve_array):
        if val < threshold:
            x = i - M
            smooth_candidates.append((x, Q_values[i]))
    return smooth_candidates

def trial_division(Q, factor_base):
    """Пробуем разложить Q на простые из факторной базы.
       Возвращает словарь {p: exp} и остаток (в случае, если остался большой множитель)."""
    exponents = {}
    temp = abs(Q)
    for p in factor_base:
        if temp == 1:
            break
        exp = 0
        while temp % p == 0:
            temp //= p
            exp += 1
        if exp > 0:
            exponents[p] = exp
    return exponents, temp

def find_dependencies(matrix):
    """Реализуем гауссовское исключение над GF(2). 
       Возвращает список решений (каждое решение – набор индексов строк, дающих зависимость)."""
    # Преобразуем матрицу в тип int
    mat = np.array(matrix, dtype=int)
    rows, cols = mat.shape
    # Индексы, участвующие в комбинациях
    dependency_sets = []
    # Будем хранить преобразования в виде двоичной матрицы (каждая строка – битовый вектор изначальной позиции)
    transform = np.eye(rows, dtype=int)
    
    col = 0
    for row in range(rows):
        if col >= cols:
            break
        pivot = None
        for i in range(row, rows):
            if mat[i, col] == 1:
                pivot = i
                break
        if pivot is None:
            col += 1
            row -= 1  # повторить итерацию для следующего столбца
            continue
        if pivot != row:
            # Меняем строки местами
            mat[[row, pivot]] = mat[[pivot, row]]
            transform[[row, pivot]] = transform[[pivot, row]]
        # Обнуляем столбец для всех остальных строк
        for i in range(rows):
            if i != row and mat[i, col] == 1:
                mat[i] ^= mat[row]
                transform[i] ^= transform[row]
        col += 1

    # Любая строка, не имеющая ведущего 1, дает зависимость
    # Здесь можно перебрать комбинации, чтобы найти ненулевые зависимости.
    # В данном примере мы просто возвращаем строки преобразования, где матрица нулевая.
    for i in range(rows):
        if not any(mat[i]):
            dependency_sets.append(np.nonzero(transform[i])[0].tolist())
    return dependency_sets

def quadratic_sieve(n):
    # Выбор параметров (значения можно оптимизировать под размер n)
    B = int(math.exp(math.sqrt(math.log(n)*math.log(math.log(n)))) * 0.18)
    factor_base = generate_factor_base(n, B)
    
    A = math.isqrt(n) + 1
    M = B  # интервал можно подобрать эмпирически; часто M ~ B
    smooth_relations = []
    xs = []  # хранение значений x
    exp_vectors = []  # векторы экспонентов по модулю 2
    
    # Сбор достаточного числа гладких соотношений
    # Требуется чуть больше соотношений, чем размер факторной базы  
    required = len(factor_base) + 5
    x0 = -M
    while len(smooth_relations) < required:
        print(f"A, B, M, factor_base: {A},{B},{M},{len(factor_base)}")
        candidates = sieve_interval(n, A, M, factor_base)
        for x, Qx in candidates:
            exponents, rem = trial_division(Qx, factor_base)
            if rem == 1:  # полностью B-гладкий
                xs.append(x)
                # Формируем вектор показателей по модулю 2
                vec = [exponents.get(p, 0) % 2 for p in factor_base]
                exp_vectors.append(vec)
                smooth_relations.append((x, Qx, exponents))
                if len(smooth_relations) >= required:
                    break
        # Если недостаточно соотношений, можно расширить интервал
        M =round(M *1.3)
    
    # Решение системы уравнений по модулю 2
    dependencies = find_dependencies(exp_vectors)
    if not dependencies:
        print("Не найдены зависимости, попробуйте расширить интервал или подобрать параметры.")
        return None
    
    # Пробуем каждую зависимость для нахождения нетривиального делителя
    for dep in dependencies:
        X = 1
        Y = 1
        combined_exponents = {}
        for i in dep:
            x, Qx, exponents = smooth_relations[i]
            X = (X * (A + x)) % n
            Y = (Y * Qx) % n
            for p, exp in exponents.items():
                combined_exponents[p] = combined_exponents.get(p, 0) + exp
        # Каждый показатель делится на 2
        Y_sqrt = 1
        for p, exp in combined_exponents.items():
            Y_sqrt = (Y_sqrt * pow(p, exp // 2, n)) % n
        factor_candidate = gcd(X - Y_sqrt, n)
        if factor_candidate not in (1, n):
            return factor_candidate, n // factor_candidate
    return None

# Пример использования:
if __name__ == "__main__":
    print()
    do = time.time()
    # Пример: число 100-бит (примерное число, около 10^30)
    # Здесь можно подставить своё число, например:
    kolichetvo_bit = 75
    n = generate_N(kolichetvo_bit)  # замените на конкретное 100-битное число
    factors = quadratic_sieve(n)
    if factors:
        print("Найденные делители:", factors)
    else:
        print("Факторизация не удалась, попробуйте изменить параметры.")
    posle = time.time()
    print(f"Время работы при {kolichetvo_bit} bit : {round(posle-do)} sec")
    print()
