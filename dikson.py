import math
import random
from collections import defaultdict

def get_B_list(B):
    """Генерация списка простых чисел до B с использованием решета Эратосфена."""
    sieve = [True] * (B + 1)
    sieve[0] = sieve[1] = False
    for p in range(2, int(math.sqrt(B)) + 1):
        if sieve[p]:
            for i in range(p * p, B + 1, p):
                sieve[i] = False
    return [p for p, is_prime in enumerate(sieve) if is_prime]

def is_gladki_B(a, b_list):
    """Проверка, является ли число a гладким относительно списка простых чисел b_list."""
    a_copy = a
    for p in b_list:
        while a_copy % p == 0:
            a_copy //= p
    return a_copy == 1

def factorize_gladki(a, b_list):
    """Факторизация числа a по списку простых чисел b_list.
       Возвращает словарь вида {p: степень}."""
    factors = defaultdict(int)
    for p in b_list:
        while a % p == 0:
            factors[p] += 1
            a //= p
    return factors

def gcd(a, b):
    """Вычисление НОД двух чисел."""
    while b:
        a, b = b, a % b
    return a

def factor_dikson(N):
    """Алгоритм факторизации Диксона (упрощённая версия)."""
    # Определяем размер факторной базы по эвристике
    B = int(round(math.exp(0.5 * math.sqrt(math.log(N) * math.log(math.log(N))))))
    b_list = get_B_list(B)  # Факторная база
    relations = []  # Список отношений: (b, факторизация a, где a = b^2 mod N)
    
    # Собираем достаточно отношений: как минимум len(b_list)+1
    while len(relations) < len(b_list) + 1:
        b = random.randint(2, N - 1)
        a = (b ** 2) % N
        if is_gladki_B(a, b_list):  # Если a полностью раскладывается по b_list
            factors = factorize_gladki(a, b_list)
            relations.append((b, factors))
    
    # Построение матрицы (взятие показателей по модулю 2)
    matrix = []
    for _, factors in relations:
        row = []
        for p in b_list:
            row.append(factors.get(p, 0) % 2)
        matrix.append(row)
    
    # Поиск двух отношений с одинаковым вектором по модулю 2.
    # Это гарантирует, что при объединении их факторизации будут иметь чётные показатели.
    for i in range(len(matrix)):
        for j in range(i + 1, len(matrix)):
            if matrix[i] == matrix[j]:
                # Объединяем два отношения
                b1, factors1 = relations[i]
                b2, factors2 = relations[j]
                # X = произведение b_i по модулю N
                X = (b1 * b2) % N
                # Объединяем факторизации
                combined_factors = defaultdict(int)
                for p in set(list(factors1.keys()) + list(factors2.keys())):
                    combined_factors[p] = factors1.get(p, 0) + factors2.get(p, 0)
                # Так как вектора совпадают по mod 2, все показатели четны.
                # Вычисляем Y = произведение p^(exponent/2) (по модулю N)
                Y = 1
                for p, exp in combined_factors.items():
                    Y = (Y * pow(p, exp // 2, N)) % N
                # Вычисляем НОД(X - Y, N)
                d = gcd(X - Y, N)
                if d != 1 and d != N:
                    return d
    return None

if __name__ == "__main__":
    N = 735422377  # Число для факторизации
    divisor = factor_dikson(N)
    if divisor:
        print(f"Найден делитель: {divisor}")
    else:
        print("Делитель не найден.")
