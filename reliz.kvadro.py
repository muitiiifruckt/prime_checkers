from math import sqrt, log2, ceil, floor, isqrt
import random
from math import gcd
from sympy import isprime
import sympy
from prime_generator import generate_N
import time 

# --------------------------
SIQS_TRIAL_DIVISION_EPS = 25  # поправка для определения порога просеивания
SIQS_MIN_PRIME_POLYNOMIAL = 400  # минимальное простое, используемое для выбора полинома
SIQS_MAX_PRIME_POLYNOMIAL = 4000  # максимальное простое для выбора полинома


# Проверка, является ли a квадратичным вычетом по модулю простого p,
# используя символ Лежандра из sympy.
is_quad = lambda a, p: sympy.legendre_symbol(a, p) == 1
sqrt_mod_prime = lambda a, p: sympy.sqrt_mod(a, p, all_roots=False)
# Вычисление обратного элемента по модулю: используем встроенный pow для Python 3.8+
modinv = lambda a, m: pow(a, -1, m)
lowest_set_bit = lambda a: (a & -a).bit_length() - 1
def find_poly(n,m,factor_base):
    
    p_min = p_max = None
    # Определяем индексы простых для выбора множителя a
    for i, factor_base_element in enumerate(factor_base):
        p = factor_base_element[0]
        if p_min is None and p >= SIQS_MIN_PRIME_POLYNOMIAL:
            p_min = i
        if p_max is None and p > SIQS_MAX_PRIME_POLYNOMIAL:
            p_max = i - 1
            break
    if p_max is None:
        p_max = len(factor_base) - 1
    if p_min is None or (p_max - p_min < 20):
        p_min = min(p_min if p_min is not None else 0, 5)
     # Определяем целевое значение для произведения выбранных простых
    target = sqrt(2 * n) / m
    target1 = target / (((factor_base[p_min][0] + factor_base[p_max][0]) / 2) ** 0.5)
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
            a *= factor_base[p_i][0]
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
        f = factor_base[q[l]]
        gamma = (f[1] * modinv(a // f[0], f[0])) % f.p
        if gamma > f[0] // 2:
            gamma = f[0] - gamma
        B.append(a // f[0] * gamma)
        
    b = sum(B) % a
    b_orig = b
    if 2 * b > a:
        b = a - b
    # Проверка корректности выбора b: должно выполняться b^2 ≡ n (mod a)
    assert 0 < b and 2 * b <= a and ((b * b - n) % a == 0)
    # Формируем полиномы: g(x) = (a*x + b)^2 - n и h(x) = a*x + b
    g = [[b * b - n, 2 * a * b, a * a], a, b_orig]
    h = [[b, a], a, b_orig]
    # Для каждого простого из факторной базы вычисляем обратный элемент и решения для просеивания
    for f in factor_base:
        if a % f.p:
            f[5] = modinv(a, f.p)
            f[3] = (f[5] * (f[1] - b)) % f.p
            f[4] = (f[5] * (-f[1] - b)) % f.p
    return g, h, B

def siqs_sieve(factor_base, m):
    """
    Функция просеивания для SIQS.
    Создается массив длины 2*m+1, в котором для каждого значения x из интервала [-m, m]
    суммируются логарифмические веса простых, для которых x является решением соответствующего сравнения.
    """
    sieve = [0] * (2 * m + 1)
    # Для каждого простого из факторной базы
    for f in factor_base:
        # Если решения для просеивания не заданы, пропускаем
        if f[3] is None:
            continue
        p = f[0]
        # Вычисляем стартовую позицию для первого решения
        i1 = -((m + f[3]) // p)
        start1 = f[3] + i1 * p + m
        # Обновляем массив просеивания для первого решения
        for a in range(start1, 2 * m + 1, p):
            sieve[a] += f[2]
        # Аналогично для второго решения
        i2 = -((m + f[4]) // p)
        start2 = f[4] + i2 * p + m
        for a in range(start2, 2 * m + 1, p):
            sieve[a] += f[1]
    return sieve
def get_small_primes(upper = 100000):
    small_primes = []
    maxi = upper
    # Первый проход: проверка простых до sqrt(upper)
    for i in range(2, maxi + 1):
        if isprime(i):
            small_primes.append(i)
            
    return small_primes

def get_factor_base(n,nf):
    small_primes = get_small_primes()
    
    factor_base = []
    for p in small_primes:
        if p == 2:  # пропускаем 2, так как legendre_symbol требует нечетный модуль
            continue
        if is_quad(n, p):
            factor_base.append([p,sqrt_mod_prime(n % p, p), log2(p),0,0,0])
            if len(factor_base) >= nf:
                break
    return factor_base
def find_poly_2(n, factor_base, i, g, B):
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
    g = [[b * b - n, 2 * a * b, a * a], a, b_orig ]
    h = [[b, a], a, b_orig]
    # Обновляем решения для просеивания для каждого простого из факторной базы
    for f in factor_base:
        if a % f.p:
            f[3] = (f[5] * (f[1] - b)) % f.p
            f[4] = (f[5] * (-f[1] - b)) % f.p
    return g, h
def siqs_trial_divide(a, fb):
    """
    Пытается разложить число a на множители, используя простые из факторной базы fb.
    Если a полностью разложимо по элементам fb, возвращает список пар (индекс, степень).
    Если нет – возвращает None.
    """
    div_idx = []
    for i, f in enumerate(fb):
        if a % f[0] == 0:
            exp = 0
            while a % f[0] == 0:
                a //= f[0]
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

def factor(n, nf = 200,  m = 65536):
    factor_base = get_factor_base(n,nf)
    
    
    req_ratio = 1.05  # начальное отношение требуемых гладких соотношений
    smooth = []       # список гладких соотношений (каждый элемент – кортеж (u, v, разложение))
    i_poly = 0        # счетчик полиномов
    prev = 0          # для отслеживания прироста количества гладких соотношений
    
    while True:
        req = round(len(factor_base) * req_ratio)
        # Цикл для поиска достаточного числа гладких соотношений
        while True:
            # Если это первая итерация, генерируем первый полином, иначе – следующий
            if i_poly == 0:
                g, h, B = find_poly(n, m, factor_base)
            else:
                g, h = find_poly_2(n, factor_base, i_poly, g, B)
            i_poly += 1
            # Сброс счетчика, если число полиномов достигло определенного предела
            if i_poly >= 2 ** (len(B) - 1):
                i_poly = 0
            # Выполняем просеивание для текущей группы полинома
            sieve = siqs_sieve(factor_base, m)
            # Пробное деление: ищем гладкие соотношения
            if siqs_trial_division(n, sieve, factor_base, smooth, g, h, m, req):
                break
    
    pass
if __name__=='__main__':
    do = time.time()
    a= factor(10000)
    print(a)
    
    
    
    posle = time.time()
    print(f"Время работы : {round(posle-do,3)} sec")
