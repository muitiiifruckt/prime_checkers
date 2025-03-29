
"""
Реализация алгоритма Self‑Initializing Quadratic Sieve (SIQS) для факторизации чисел.
Данный алгоритм использует факторную базу, генерацию полиномов, просеивание,
а затем линейную алгебру для поиска зависимостей между найденными гладкими соотношениями.
"""

from math import sqrt, log2, ceil, floor, isqrt
import random
from math import gcd
import sympy

# --------------------------
# Параметры алгоритма
# --------------------------
SIQS_TRIAL_DIVISION_EPS = 25  # поправка для определения порога просеивания
SIQS_MIN_PRIME_POLYNOMIAL = 400  # минимальное простое, используемое для выбора полинома
SIQS_MAX_PRIME_POLYNOMIAL = 4000  # максимальное простое для выбора полинома

# --------------------------
# Классы, используемые в алгоритме
# --------------------------
class FactorBasePrime:
    """
    Класс для хранения простых чисел, входящих в факторную базу.
    Для каждого простого числа p сохраняются:
      - p: само простое число;
      - tmem: значение квадратного корня из n по модулю p (один из корней);
      - lp: приближённое значение log2(p);
      - soln1, soln2: решения сравнения a*x + b (вычисляются для просеивания);
      - ainv: обратный элемент a по модулю p, где a используется в полиномах.
    """
    def __init__(self, p, tmem, lp):
        self.p = p
        self.tmem = tmem
        self.lp = round(log2(p))
        self.soln1 = None
        self.soln2 = None
        self.ainv = None
    def __str__(self):
        return f"{self.p}, {self.tmem}, {self.lp}"
  
class Polynomial:
    """
    Класс для представления квадратичного полинома вида
        f(x) = (a*x + b)^2 - n
    coeff: коэффициенты полинома (например, [b^2 - n, 2*a*b, a^2])
    a, b: параметры, используемые при вычислении полинома.
    Метод eval(x) вычисляет значение полинома для заданного x.
    """
    def __init__(self, coeff, a, b):
        self.coeff = coeff  # список коэффициентов
        self.a = a          # коэффициент a
        self.b = b          # коэффициент b (оригинальное значение b, до возможного преобразования)
    
    def eval(self, x):
        res = 0
        # Вычисляем значение полинома, используя схему Горнера
        for c in reversed(self.coeff):
            res = res * x + c
        return res

# --------------------------
# Вспомогательные функции
# --------------------------
# Функция для вычисления квадратного корня по модулю p с использованием sympy.
# all_roots=False возвращает один из возможных корней.
sqrt_mod_prime = lambda a, p: sympy.sqrt_mod(a, p, all_roots=False)

# Вычисление обратного элемента по модулю: используем встроенный pow для Python 3.8+
modinv = lambda a, m: pow(a, -1, m)

# Проверка, является ли a квадратичным вычетом по модулю простого p,
# используя символ Лежандра из sympy.
is_quad = lambda a, p: sympy.legendre_symbol(a, p) == 1

# Функция lowest_set_bit возвращает индекс младшего установленного (единичного) бита числа a.
lowest_set_bit = lambda a: (a & -a).bit_length() - 1

# --------------------------
# Функции для вычисления зависимостей и делителей
# --------------------------
def siqs_calc_sqrts(indices, rel):
    """
    Вычисляет пару чисел [X, Y] из набора гладких соотношений.
    Для набора индексов indices произведение чисел u (первый элемент в каждой тройке)
    даёт X, а произведение чисел v (второй элемент) даёт Y.
    Y затем извлекается с помощью целочисленного квадратного корня.
    """
    X = 1
    Y = 1
    for i in indices:
        X *= rel[i][0]
        Y *= rel[i][1]
    return [X, isqrt(Y)]

def siqs_factor_from_square(n, indices, rel):
    """
    Из набора гладких соотношений (определённых индексами indices) вычисляет делитель n.
    Фактор f = gcd(|X - Y|, n), где [X, Y] получены функцией siqs_calc_sqrts.
    """
    s1, s2 = siqs_calc_sqrts(indices, rel)
    return gcd(abs(s1 - s2), n)

def siqs_find_factors(n, squares, rel):
    """
    Из набора зависимостей (squares) и гладких соотношений (rel) пытается извлечь
    нетривиальные делители числа n. Возвращает список найденных делителей.
    """
    fac = []
    rem = n
    for inds in squares:
        f = siqs_factor_from_square(n, inds, rel)
        if 1 < f < rem:
            # Пока делитель f делит оставшееся число, добавляем его в список
            while rem % f == 0:
                fac.append(f)
                rem //= f
    if rem > 1:
        fac.append(rem)
    return fac

# --------------------------
# Формирование факторной базы
# --------------------------
def siqs_factor_base_primes(n, nf, small_primes):
    """
    Формирует факторную базу из nf простых чисел.
    Для каждого простого p (кроме 2, которая пропускается) проверяется,
    является ли n квадратичным вычетом по модулю p.
    Если да, то вычисляется sqrt_mod_prime(n mod p, p) и p добавляется в факторную базу.
    """
    fb = []
    for p in small_primes:
        if p == 2:  # пропускаем 2, так как legendre_symbol требует нечетный модуль
            continue
        if is_quad(n, p):
            fb.append(FactorBasePrime(p, sqrt_mod_prime(n % p, p), log2(p)))
            if len(fb) >= nf:
                break
    return fb

# --------------------------
# Генерация полиномов для просеивания
# --------------------------
def siqs_find_first_poly(n, m, fb):
    """
    Генерирует первый полином для просеивания.
    Выбирает подмножество индексов факторной базы, удовлетворяющих условиям
    (простые числа в диапазоне [SIQS_MIN_PRIME_POLYNOMIAL, SIQS_MAX_PRIME_POLYNOMIAL]).
    На основании выбранных простых вычисляется значение a, затем находится b такое, что
         b^2 ≡ n (mod a)
    Возвращает сгенерированные полиномы g и h, а также список B, необходимый для генерации следующих полиномов.
    """
    p_min = p_max = None
    # Определяем индексы простых для выбора множителя a
    for i, f in enumerate(fb):
        if p_min is None and f.p >= SIQS_MIN_PRIME_POLYNOMIAL:
            p_min = i
        if p_max is None and f.p > SIQS_MAX_PRIME_POLYNOMIAL:
            p_max = i - 1
            break
    if p_max is None:
        p_max = len(fb) - 1
    if p_min is None or (p_max - p_min < 20):
        p_min = min(p_min if p_min is not None else 0, 5)
    # Определяем целевое значение для произведения выбранных простых
    target = sqrt(2 * n) / m
    target1 = target / (((fb[p_min].p + fb[p_max].p) / 2) ** 0.5)
    best_a = best_ratio = None
    best_q = None
    # Проводим несколько попыток для выбора оптимального множителя a
    for _ in range(30):
        a = 1
        q = []
        while a < target1:
            p_i = random.randint(p_min, p_max)
            if p_i in q:
                continue
            a *= fb[p_i].p
            q.append(p_i)
        ratio = a / target
        if best_ratio is None or (ratio >= 0.9 and ratio < best_ratio) or (best_ratio < 0.9 and ratio > best_ratio):
            best_ratio = ratio
            best_a = a
            best_q = q
    a, q = best_a, best_q
    s = len(q)
    B = []
    # Для каждого выбранного простого вычисляем вспомогательное значение B[l]
    for l in range(s):
        f = fb[q[l]]
        gamma = (f.tmem * modinv(a // f.p, f.p)) % f.p
        if gamma > f.p // 2:
            gamma = f.p - gamma
        B.append(a // f.p * gamma)
    b = sum(B) % a
    b_orig = b
    if 2 * b > a:
        b = a - b
    # Проверка корректности выбора b: должно выполняться b^2 ≡ n (mod a)
    assert 0 < b and 2 * b <= a and ((b * b - n) % a == 0)
    # Формируем полиномы: g(x) = (a*x + b)^2 - n и h(x) = a*x + b
    g = Polynomial([b * b - n, 2 * a * b, a * a], a, b_orig)
    h = Polynomial([b, a], a, b_orig)
    # Для каждого простого из факторной базы вычисляем обратный элемент и решения для просеивания
    for f in fb:
        if a % f.p:
            f.ainv = modinv(a, f.p)
            f.soln1 = (f.ainv * (f.tmem - b)) % f.p
            f.soln2 = (f.ainv * (-f.tmem - b)) % f.p
    return g, h, B

def siqs_find_next_poly(n, fb, i, g, B):
    """
    Генерирует следующий полином для просеивания, основываясь на предыдущем полиноме g и вспомогательном массиве B.
    Используется значение i для определения корректировки b.
    Возвращает новый полином g и соответствующий ему полином h.
    """
    v = lowest_set_bit(i) + 1
    # z определяется в зависимости от четности (используем ceil для вычислений)
    z = -1 if ceil(i / (2 ** v)) % 2 == 1 else 1
    # Вычисляем новое значение b с поправкой
    b = (g.b + 2 * z * B[v - 1]) % g.a
    a = g.a
    b_orig = b
    if 2 * b > a:
        b = a - b
    assert ((b * b - n) % a == 0)
    # Формируем новый полином g и соответствующий h
    g = Polynomial([b * b - n, 2 * a * b, a * a], a, b_orig)
    h = Polynomial([b, a], a, b_orig)
    # Обновляем решения для просеивания для каждого простого из факторной базы
    for f in fb:
        if a % f.p:
            f.soln1 = (f.ainv * (f.tmem - b)) % f.p
            f.soln2 = (f.ainv * (-f.tmem - b)) % f.p
    return g, h

# --------------------------
# Просеивание (sieving)
# --------------------------
def siqs_sieve(fb, m):
    """
    Функция просеивания для SIQS.
    Создается массив длины 2*m+1, в котором для каждого значения x из интервала [-m, m]
    суммируются логарифмические веса простых, для которых x является решением соответствующего сравнения.
    """
    sieve = [0] * (2 * m + 1)
    # Для каждого простого из факторной базы
    for f in fb:
        # Если решения для просеивания не заданы, пропускаем
        if f.soln1 is None:
            continue
        p = f.p
        # Вычисляем стартовую позицию для первого решения
        i1 = -((m + f.soln1) // p)
        start1 = f.soln1 + i1 * p + m
        # Обновляем массив просеивания для первого решения
        for a in range(start1, 2 * m + 1, p):
            sieve[a] += f.lp
        # Аналогично для второго решения
        i2 = -((m + f.soln2) // p)
        start2 = f.soln2 + i2 * p + m
        for a in range(start2, 2 * m + 1, p):
            sieve[a] += f.lp
    return sieve

def siqs_trial_divide(a, fb):
    """
    Пытается разложить число a на множители, используя простые из факторной базы fb.
    Если a полностью разложимо по элементам fb, возвращает список пар (индекс, степень).
    Если нет – возвращает None.
    """
    div_idx = []
    for i, f in enumerate(fb):
        if a % f.p == 0:
            exp = 0
            while a % f.p == 0:
                a //= f.p
                exp += 1
            div_idx.append((i, exp))
        if a == 1:
            return div_idx
    return None

def siqs_trial_division(n, sieve, fb, smooth, g, h, m, req):
    """
    Пробное деление: для каждого x из просеянного интервала [-m, m],
    если сумма логарифмов (sieve[x+m]) превышает порог, вычисляем значение полинома g(x).
    Затем, если g(x) полностью разлагается по элементам факторной базы fb,
    добавляем найденное гладкое соотношение (u, v и его факторизацию) в список smooth.
    Если количество гладких соотношений достигает req – возвращает True.
    """
    # Вычисляем пороговое значение для гладкости
    limit = log2(m * sqrt(n)) - SIQS_TRIAL_DIVISION_EPS
    for i, sa in enumerate(sieve):
        if sa >= limit:
            x = i - m
            gx = g.eval(x)
            d = siqs_trial_divide(gx, fb)
            if d is not None:
                u = h.eval(x)
                v = gx
                smooth.append((u, v, d))
                if len(smooth) >= req:
                    return True
                    
    return False

# --------------------------
# Линейная алгебра для нахождения зависимостей
# --------------------------
def siqs_build_matrix(fb, smooth):
    """
    Строит бинарную матрицу (по модулю 2) для линейной алгебры.
    Каждая строка соответствует гладкому соотношению, а каждый столбец – простому из факторной базы.
    Элемент равен 1, если степень данного простого нечётна, иначе 0.
    """
    M = []
    for u, v, div in smooth:
        row = [0] * len(fb)
        for j, exp in div:
            row[j] = exp % 2
        M.append(row)
    return M

def siqs_build_matrix_opt(M):
    """
    Преобразует матрицу M в список чисел, представляющих столбцы матрицы.
    Каждый столбец кодируется в двоичном виде как число.
    """
    m = len(M[0])
    cols = ["" for _ in range(m)]
    for row in M:
        for j, bit in enumerate(row):
            cols[j] += "1" if bit else "0"
    # Переворачиваем строку, чтобы младший бит соответствовал первому элементу строки
    return [int(c[::-1], 2) for c in cols], len(M), m

def add_column_opt(M_opt, tgt, src):
    """
    Добавляет (по модулю 2) столбец src к столбцу tgt в оптимизированном представлении.
    """
    M_opt[tgt] ^= M_opt[src]

def find_pivot_column_opt(M_opt, j):
    """
    Для столбца j в оптимизированной матрице возвращает индекс первого ненулевого бита,
    который и будет являться опорным (pivot) для этого столбца.
    Если столбец равен 0, возвращает None.
    """
    return None if M_opt[j] == 0 else lowest_set_bit(M_opt[j])

def siqs_solve_matrix_opt(M_opt, n, m):
    """
    Решает систему линейных уравнений над GF(2) для нахождения зависимостей между строками.
    В результате возвращает список индексов гладких соотношений, из которых можно
    получить квадратное сравнение, способное дать делитель n.
    """
    row_mark = [False] * n  # метки для строк, участвующих в опорных позициях
    pivots = [-1] * m       # для каждого столбца хранится индекс опорной строки
    for j in range(m):
        i = find_pivot_column_opt(M_opt, j)
        if i is not None:
            pivots[j] = i
            row_mark[i] = True
            for k in range(m):
                if k != j and ((M_opt[k] >> i) & 1):
                    add_column_opt(M_opt, k, j)
    squares = []
    # Собираем группы строк, не участвующих в опорных позициях – они дают зависимости
    for i in range(n):
        if not row_mark[i]:
            inds = [i]
            for j in range(m):
                if (M_opt[j] >> i) & 1:
                    inds.append(pivots[j])
            squares.append(inds)
    return squares

# --------------------------
# Выбор параметров просеивания
# --------------------------
def siqs_choose_nf_m(d):
    """
    Выбирает размер факторной базы (nf) и интервал просеивания (m) в зависимости от количества цифр числа n.
    Эти параметры подобраны экспериментально и могут быть изменены для оптимизации.
    """
    if d <= 34: return 200, 65536
    if d <= 36: return 300, 65536
    if d <= 38: return 400, 65536
    if d <= 40: return 500, 65536
    if d <= 42: return 600, 65536
    if d <= 44: return 700, 65536
    if d <= 48: return 1000, 65536
    if d <= 52: return 1200, 65536
    if d <= 56: return 2000, 65536 * 3
    if d <= 60: return 4000, 65536 * 3
    if d <= 66: return 6000, 65536 * 3
    if d <= 74: return 10000, 65536 * 3
    if d <= 80: return 30000, 65536 * 3
    if d <= 88: return 50000, 65536 * 3
    if d <= 94: return 60000, 65536 * 9
    return 100000, 65536 * 9

# --------------------------
# Инициализация маленьких простых чисел с помощью решета Эратосфена
# --------------------------
def trial_div_init_primes(n, upper):
    """
    Выполняет пробное деление числа n, используя все простые до upper.
    Одновременно инициализирует глобальную переменную small_primes – список маленьких простых.
    Возвращает список найденных малых делителей и остаток, если он не равен 1.
    """
    global small_primes
    is_prime = [True] * (upper + 1)
    is_prime[0:2] = [False, False]
    fac = []
    small_primes = []
    rem = n
    maxi = isqrt(upper)
    # Первый проход: проверка простых до sqrt(upper)
    for i in range(2, maxi + 1):
        if is_prime[i]:
            small_primes.append(i)
            while rem % i == 0:
                rem //= i
                fac.append(i)
                # Если остаток является простым, возвращаем его сразу
                if sympy.isprime(rem):
                    fac.append(rem)
                    return fac, 1
            for j in range(i * i, upper + 1, i):
                is_prime[j] = False
    return fac, rem

# --------------------------
# Основная функция SIQS
# --------------------------
def siqs_factorise(n):
    """
    Основная функция факторизации числа n с использованием алгоритма SIQS.
    1. Определяются параметры nf (размер факторной базы) и m (интервал просеивания).
    2. Формируется факторная база.
    3. Ищутся гладкие соотношения с помощью просеивания и пробного деления.
    4. Построение матрицы для линейной алгебры и поиск зависимостей.
    5. Вычисление квадратных корней и получение делителей через НОД.
    Если факторизация успешна (найдено более одного делителя), возвращаются делители.
    В противном случае увеличивается требуемое число гладких соотношений и поиск повторяется.
    """
    d = len(str(n))
    nf, m = siqs_choose_nf_m(d)
    # Формируем факторную базу с использованием глобального списка small_primes
    fb = siqs_factor_base_primes(n, nf, small_primes)
    req_ratio = 1.05  # начальное отношение требуемых гладких соотношений
    smooth = []       # список гладких соотношений (каждый элемент – кортеж (u, v, разложение))
    i_poly = 0        # счетчик полиномов
    prev = 0          # для отслеживания прироста количества гладких соотношений
    # Основной цикл поиска зависимостей
    while True:
        req = round(len(fb) * req_ratio)
        # Цикл для поиска достаточного числа гладких соотношений
        while True:
            # Если это первая итерация, генерируем первый полином, иначе – следующий
            if i_poly == 0:
                g, h, B = siqs_find_first_poly(n, m, fb)
            else:
                g, h = siqs_find_next_poly(n, fb, i_poly, g, B)
            i_poly += 1
            # Сброс счетчика, если число полиномов достигло определенного предела
            if i_poly >= 2 ** (len(B) - 1):
                i_poly = 0
            # Выполняем просеивание для текущей группы полинома
            sieve = siqs_sieve(fb, m)
            # Пробное деление: ищем гладкие соотношения
            if siqs_trial_division(n, sieve, fb, smooth, g, h, m, req):
                break
        # Построение бинарной матрицы для линейной алгебры
        M = siqs_build_matrix(fb, smooth)
        
        M_opt, M_n, M_m = siqs_build_matrix_opt(M)
        
        # Решаем систему линейных уравнений для поиска зависимостей
        squares = siqs_solve_matrix_opt(M_opt, M_n, M_m)
        
        # Пытаемся получить делители из найденных зависимостей
        fac = siqs_find_factors(n, squares, smooth)
        if len(fac) > 1:
            return fac
        # Если делитель не найден, увеличиваем требуемое число гладких соотношений и повторяем поиск
        req_ratio += 0.05


if __name__=='__main__':
    from prime_generator import generate_N
    import time 
    do = time.time()
    kolichetvo_bit = 120
    n = generate_N(kolichetvo_bit)
    #n = 3485198031263125658050790686801
    print(n)
    fac, _ = trial_div_init_primes(n, 1000000)
    if fac:
        print(fac)
        exit()
    print(siqs_factorise(n))
    posle = time.time()
    print(f"Время работы при {kolichetvo_bit} bit : {round(posle-do,3)} sec")
