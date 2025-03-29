import math
import numpy as np
from math import log, isqrt, gcd
from prime_generator import generate_N
import time


def legendre_symbol(a, p):
    """Вычисляет символ Лежандра (a/p)."""
    return pow(a, (p - 1) // 2, p)

def generate_factor_base(n, B):
    """
    Формирует факторную базу: простые p ≤ B, для которых (n/p)=1.
    Используется векторизация для решета с помощью NumPy.
    """
    sieve = np.ones(B+1, dtype=bool)
    sieve[:2] = False
    for p in range(2, int(B**0.5)+1):
        if sieve[p]:
            sieve[p*p:B+1:p] = False
    primes = np.nonzero(sieve)[0]
    # Фильтрация простых, удовлетворяющих (n/p)==1
    factor_base = [int(p) for p in primes if legendre_symbol(n, int(p)) == 1]
    return factor_base

def quadratic_polynomial(A, x, n):
    """Возвращает значение Q(x) = (A+x)^2 - n."""
    return (A + x)**2 - n

def sieve_interval(n, A, M, factor_base):
    """
    Сегментированное решето для интервала x в [-M, M].
    Вычисляет Q(x) = (A+x)^2 - n и корректирует логарифмические значения для поиска гладких чисел.
    
    Если n достаточно мало (bit_length < 63), используется векторизация с NumPy.
    Для больших чисел (например, до 120 бит) применяется режим со списковыми вычислениями с оптимизированными циклами.
    """
    threshold = 0.5  # порог для выбора кандидатов

    if n.bit_length() < 63:
        # Векторизированный режим с использованием np.int64
        xs = np.arange(-M, M+1, dtype=np.int64)
        Q_values = (A + xs)**2 - n
        abs_Q = np.abs(Q_values)
        log_Q = np.where(abs_Q == 0, 0.0, np.log(abs_Q.astype(np.float64)))
        sieve_array = log_Q.copy()
        
        for p in factor_base:
            logp = math.log(p)
            ts = np.arange(p)
            sols = ts[((ts * ts - n) % p) == 0]
            for t in sols:
                x0 = (t - A) % p
                # Вычисляем k_min и k_max точно, чтобы x = x0 + k*p попадало в [-M, M]
                k_min = int(np.ceil((-M - x0) / p))
                k_max = int(np.floor((M - x0) / p))
                # Вычисляем индексы сразу (без проверки, поскольку границы гарантируют попадание в интервал)
                indices = x0 + np.arange(k_min, k_max+1) * p + M
                sieve_array[indices.astype(int)] -= logp
        
        candidate_indices = np.where(sieve_array < threshold)[0]
        smooth_candidates = [(int(i - M), Q_values[i]) for i in candidate_indices]
    
    else:
        # Режим для больших чисел с использованием списковых вычислений
        xs = list(range(-M, M+1))
        Q_values = [(A + x)**2 - n for x in xs]
        log_Q = [0.0 if q == 0 else math.log(abs(q)) for q in Q_values]
        sieve_array = log_Q.copy()
        
        for p in factor_base:
            logp = math.log(p)
            # Решаем уравнение t^2 ≡ n (mod p) перебором по t в [0, p-1]
            sols = [t for t in range(p) if (t*t - n) % p == 0]
            for t in sols:
                x0 = (t - A) % p
                # Вычисляем границы k таким образом, чтобы x = x0 + k*p попадало в [-M, M]
                k_min = int(math.ceil((-M - x0) / p))
                k_max = int(math.floor((M - x0) / p))
                # Перебор по k без лишней проверки
                for k in range(k_min, k_max+1):
                    idx = x0 + k * p + M
                    sieve_array[idx] -= logp
        
        smooth_candidates = [(i - M, Q_values[i]) for i, v in enumerate(sieve_array) if v < threshold]
    
    return smooth_candidates


def trial_division(Q, factor_base):
    """
    Пытается разложить Q на простые из факторной базы.
    Возвращает словарь {p: exp} и остаток (если остался большой множитель).
    """
    exponents = {}
    temp = abs(Q)
    for p in factor_base:
        if temp == 1:
            break
        exp = 0
        while temp % p == 0:
            temp //= p
            exp += 1
        if exp:
            exponents[p] = exp
    return exponents, temp

def find_dependencies(exp_vectors):
    """
    Гауссово исключение по модулю 2 с использованием битовых масок.
    Каждая строка представлена как целое число, где i-й бит соответствует элементу.
    Возвращает список зависимостей (наборы индексов строк, дающих нулевой вектор).
    """
    n_rows = len(exp_vectors)
    n_cols = len(exp_vectors[0])
    # Преобразуем каждую строку в битовую маску
    bit_rows = []
    for vec in exp_vectors:
        bits = 0
        for bit in vec:
            bits = (bits << 1) | (bit & 1)
        bit_rows.append(bits)
    
    # Запоминаем преобразования: для каждой строки сохраняем битовую маску, где установлен бит i означает, что строка i участвует в комбинации
    transforms = [1 << i for i in range(n_rows)]
    
    pivot_cols = {}
    for col in range(n_cols):
        pivot = None
        for i in range(n_rows):
            # Если столбец уже обработан пропускаем строки без установленного нужного бита
            if (bit_rows[i] >> (n_cols - 1 - col)) & 1:
                if i not in pivot_cols.values():
                    pivot = i
                    break
        if pivot is None:
            continue
        pivot_cols[col] = pivot
        for i in range(n_rows):
            if i != pivot and ((bit_rows[i] >> (n_cols - 1 - col)) & 1):
                bit_rows[i] ^= bit_rows[pivot]
                transforms[i] ^= transforms[pivot]
    
    dependencies = []
    full_mask = (1 << n_cols) - 1
    for i in range(n_rows):
        if bit_rows[i] == 0:
            # Извлекаем индексы строк из битовой маски
            dep = [j for j in range(n_rows) if (transforms[i] >> j) & 1]
            dependencies.append(dep)
    return dependencies

def quadratic_sieve(n):
    """
    Основная функция квадратичного решета.
    Подбирает параметры, собирает соотношения и ищет нетривиальные делители.
    """
    # Подбор параметров; параметр B можно дополнительно оптимизировать
    B = int(math.exp(math.sqrt(math.log(n) * math.log(math.log(n)))) * 0.18)
    factor_base = generate_factor_base(n, B)
    
    A = isqrt(n) + 1
    M = B  # начальный интервал, может быть увеличен
    smooth_relations = []
    xs = []       # храним значения x
    exp_vectors = []  # векторы показателей (по модулю 2)
    
    required = len(factor_base) + 5
    while len(smooth_relations) < required:
        # Можно отключить вывод для повышения производительности
        # print(f"A={A}, B={B}, M={M}, факторов в базе: {len(factor_base)}")
        candidates = sieve_interval(n, A, M, factor_base)
        for x, Qx in candidates:
            exponents, rem = trial_division(Qx, factor_base)
            if rem == 1:  # если Qx полностью B-гладкий
                xs.append(x)
                vec = [exponents.get(p, 0) % 2 for p in factor_base]
                exp_vectors.append(vec)
                smooth_relations.append((x, Qx, exponents))
                if len(smooth_relations) >= required:
                    break
        M = round(M * 1.3)  # расширяем интервал при необходимости
    
    dependencies = find_dependencies(exp_vectors)
    if not dependencies:
        print("Не найдены зависимости, попробуйте расширить интервал или подобрать другие параметры.")
        return None
    
    # Используем найденные зависимости для получения делителей
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
    kolichetvo_bit = 70
    n = generate_N(kolichetvo_bit)  # замените на конкретное 100-битное число
    #n =729148987615913831
    print(n)
    factors = quadratic_sieve(n)
    if factors:
        print("Найденные делители:", factors)
    else:
        print("Факторизация не удалась, попробуйте изменить параметры.")
    posle = time.time()
    print(f"Время работы при {kolichetvo_bit} bit : {round(posle-do,3)} sec")
    print()
