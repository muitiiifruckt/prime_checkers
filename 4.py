# Объединённая программа для лабораторных работ по криптоанализу
# (Лабораторная работа 3 – построение таблицы линейной аппроксимации S-бокса;
#  Лабораторная работа 4 – взлом 2-х раундовой системы шифрования)
#
# Оформление вывода и входные данные приведены согласно методичке.

def int_to_bits(x, n=4):
    """Преобразует число x в список битов длины n (от старшего к младшему)."""
    return [(x >> i) & 1 for i in reversed(range(n))]


def bits_to_int(bits):
    """Преобразует список битов (от старшего к младшему) в число."""
    x = 0
    for b in bits:
        x = (x << 1) | b
    return x


def permutation(bits):
    """Перестановка P: циклический сдвиг вправо для 4-битового слова.
       [b0, b1, b2, b3] → [b3, b0, b1, b2]
       что соответствует перестановке P = [2,3,4,1] (при нумерации с 1)."""
    return [bits[-1]] + bits[:-1]


def permutation_inv(bits):
    """Обратная перестановка P⁻¹: циклический сдвиг влево.
       [b0, b1, b2, b3] → [b1, b2, b3, b0]"""
    return bits[1:] + [bits[0]]


def parity(n: int) -> int:
    """Возвращает чётность числа n (0, если число единичных битов чётно, иначе 1)."""
    return bin(n).count("1") % 2


def int_to_bin_str(x: int, width: int = 4) -> str:
    """Преобразует число в строку двоичного представления фиксированной длины."""
    return format(x, '0{}b'.format(width))


# ================================
# Часть 1. Построение таблицы линейной аппроксимации (LAB‑3)
# ================================

def build_sbox(F1, F2, F3, F4):
    """Строит S‑бокс как список целых чисел от 0 до 15 на основе базисных функций F1...F4.
       S(x) = F1(x)*2^3 + F2(x)*2^2 + F3(x)*2^1 + F4(x)*2^0.
    """
    S = []
    for i in range(len(F1)):
        val = (F1[i] << 3) | (F2[i] << 2) | (F3[i] << 1) | (F4[i] << 0)
        S.append(val)
    return S


def compute_lat(sbox):
    """Вычисляет таблицу линейной аппроксимации (LAT) для S‑бокса.
       Для масок α и β (от 0 до 15):
       LAT(α,β) = (число x, для которых parity(α·x)==parity(β·S(x))) - 2^(n-1)
       при n=4 => вычитаем 8.
    """
    n = 4
    N = 2 ** n
    lat = [[0] * N for _ in range(N)]
    for a in range(N):
        for b in range(N):
            count = 0
            for x in range(N):
                px = parity(a & x)
                py = parity(b & sbox[x])
                if px == py:
                    count += 1
            lat[a][b] = count - (N // 2)
    return lat


def print_lat_table(lat):
    """Выводит таблицу LAT с выравниванием столбцов (маски выводятся в виде 4‑битных строк)."""
    N = len(lat)
    col_width = 4
    print("\nТаблица линейной аппроксимации (LAT):")
    print("        ", end="")
    for b in range(N):
        print(f"{int_to_bin_str(b):>{col_width}}", end=" ")
    print()
    print("      " + "-" * ((col_width + 1) * N))
    for a in range(N):
        print(f"{int_to_bin_str(a):>4} |", end=" ")
        for b in range(N):
            print(f"{lat[a][b]:>{col_width}}", end=" ")
        print()


def find_linear_approximations(lat, threshold=4):
    """Находит все пары (α,β) с |LAT(α,β)| >= threshold и формирует строковое представление линейных аппроксимаций.
       Если LAT(α,β) >= 0, аппроксимация имеет вид:
           (α·x) = (β·S(x))
       Если отрицательно – эквивалентно:
           (α·x) = (β·S(x)) ⊕ 1
       Вероятность совпадения = (8 + |LAT|)/16.
    """
    N = len(lat)
    approximations = []
    for a in range(N):
        for b in range(N):
            if abs(lat[a][b]) >= threshold and not (a == 0 and b == 0):
                if lat[a][b] >= 0:
                    eq_str = f"( {int_to_bin_str(a)} · x ) = ( {int_to_bin_str(b)} · S(x) )"
                else:
                    eq_str = f"( {int_to_bin_str(a)} · x ) = ( {int_to_bin_str(b)} · S(x) ) ⊕ 1"
                prob = (8 + abs(lat[a][b])) / 16
                approximations.append((a, b, eq_str, lat[a][b], prob))
    return approximations


def print_linear_approximations(approximations):
    """Выводит найденную систему линейных уравнений, описывающих аппроксимацию S‑бокса."""
    if not approximations:
        print("Нет найденных линейных аппроксимаций с заданным порогом отклонения.")
    else:
        print("\nНайденные линейные аппроксимации:")
        for a, b, eq_str, bias, prob in approximations:
            print(f"{eq_str}   (LAT = {bias:>3}, вероятность = {prob:.3f})")


# ================================
# Часть 2. Взлом 2‑раундовой системы шифрования (LAB‑4)
# ================================

def sbox_actual(x, S_box):
    """Настоящий s‑бокс: на входе 4‑бит (0..15), возвращает S_box[x]."""
    return S_box[x]


def sbox_actual_inv(x, S_box, S_inv):
    """Обратный s‑бокс: возвращает x такой, что S_box[S_inv[x]] == x."""
    return S_inv[x]


def linear_approximation1(x, S_box):
    """Линейная аппроксимация s‑бокса (функция L(x) для тестовой системы).
       Вычисляет битовые соотношения с использованием четырёх масок:
         (0001 · x) = (0001 · S(x))
         (0011 · x) = (1000 · S(x))
         (0100 · x) = (0100 · S(x))
         (1011 · x) = (0010 · S(x)) ⊕ 1
       Результат – XOR всех четырёх вычисленных битов.
    """
    input_mask1 = 0b0001;
    output_mask1 = 0b0001
    result1 = (bin(x & input_mask1).count('1') % 2) ^ (bin(S_box[x] & output_mask1).count('1') % 2)

    input_mask2 = 0b0011;
    output_mask2 = 0b1000
    result2 = (bin(x & input_mask2).count('1') % 2) ^ (bin(S_box[x] & output_mask2).count('1') % 2)

    input_mask3 = 0b0100;
    output_mask3 = 0b0100
    result3 = (bin(x & input_mask3).count('1') % 2) ^ (bin(S_box[x] & output_mask3).count('1') % 2)

    input_mask4 = 0b1011;
    output_mask4 = 0b0010
    result4 = ((bin(x & input_mask4).count('1') % 2) ^ (bin(S_box[x] & output_mask4).count('1') % 2)) ^ 1

    return result1 ^ result2 ^ result3 ^ result4


def round_encryption_actual(x, K, S_box):
    """Один раунд настоящей системы:
       1. Подставляем x в настоящий s‑бокс.
       2. Складываем (XOR) с ключом K.
       3. Применяем перестановку P.
    """
    y = sbox_actual(x, S_box)
    z = y ^ K
    z_bits = int_to_bits(z, 4)
    out_bits = permutation(z_bits)
    return bits_to_int(out_bits)


def encrypt_actual(x, K, rounds, S_box):
    """Шифрование настоящей системы (количество раундов задаётся параметром rounds)."""
    y = x
    for _ in range(rounds):
        y = round_encryption_actual(y, K, S_box)
    return y


def round_encryption_test(x, K, S_box):
    """Один раунд тестовой системы:
       1. Вычисляем линейную аппроксимацию L(x).
       2. Складываем (XOR) с ключом K.
       3. Применяем перестановку P.
    """
    y = linear_approximation1(x, S_box)
    z = y ^ K
    z_bits = int_to_bits(z, 4)
    out_bits = permutation(z_bits)
    return bits_to_int(out_bits)


def encrypt_test(x, K, rounds, S_box):
    """Шифрование тестовой системы (на основе линейной аппроксимации)."""
    y = x
    for _ in range(rounds):
        y = round_encryption_test(y, K, S_box)
    return y


def round_decryption_actual(y, K, S_box, S_inv):
    """Один обратный раунд настоящей системы:
       1. Применяем обратную перестановку.
       2. Вычитаем ключ (XOR).
       3. Применяем обратный s‑бокс.
    """
    y_bits = int_to_bits(y, 4)
    pre_bits = permutation_inv(y_bits)
    pre_val = bits_to_int(pre_bits)
    x = sbox_actual_inv(pre_val ^ K, S_box, S_inv)
    return x


def decrypt_actual(y, K, rounds, S_box, S_inv):
    """Расшифрование настоящей системы (количество раундов задаётся параметром rounds)."""
    x = y
    for _ in range(rounds):
        x = round_decryption_actual(x, K, S_box, S_inv)
    return x


# ================================
# Главная процедура
# ================================
def main():
    # ЗАДАННЫЕ ВХОДНЫЕ ДАННЫЕ (вариант, используемый в данной программе)
    # Базисные функции S-бокса (длина 16)
    F1 = [0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1]
    F2 = [0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1]
    F3 = [1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0]
    F4 = [0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]

    # Проверяем, что все списки имеют длину 16
    assert len(F1) == 16 and len(F2) == 16 and len(F3) == 16 and len(F4) == 16, "Длина базисных функций должна быть 16."

    # Частичная шифровальная таблица (для x от 0000 до 1011)
    partial_table = {
        0: 5, 1: 11, 2: 7, 3: 9,
        4: 12, 5: 14, 6: 2, 7: 12,
        8: 8, 9: 2, 10: 10, 11: 3  # Дополнительное значение (не используется при переборе ключей, x=12)
    }
    # Для расчёта ошибок используем значения для x=0..11
    partial_values = [partial_table[x] for x in range(12)]

    # Заданная перестановка P – циклический сдвиг вправо, что соответствует P=[2,3,4,1]
    P_str = "[4,3,1,2]"

    # Зашифрованное сообщение (для m от 12 до 15)
    c_message = 3  # согласно тестовой задаче

    # Вывод заголовка и исходных данных
    print("Пример решения лабораторной работы\n")
    print("Исходные данные:")
    print("F1 = (" + ",".join(str(bit) for bit in F1) + ")")
    print("F2 = (" + ",".join(str(bit) for bit in F2) + ")")
    print("F3 = (" + ",".join(str(bit) for bit in F3) + ")")
    print("F4 = (" + ",".join(str(bit) for bit in F4) + ")\n")

    print("Частичная шифровальная таблица 2-раундового шифра от X=0000 до X=1011:")
    pt_str = "(" + ",".join(str(partial_table[x]) for x in range(12)) + ")"
    print("C = " + pt_str + "\n")

    print("Перестановка P = " + P_str)
    print(f"Зашифрованное сообщение c(m) = {c_message} (при 12 <= m < 16)\n")

    # Построение S-бокса
    S_box = build_sbox(F1, F2, F3, F4)
    print("S-box (для x от 0000 до 1111):")
    for x in range(16):
        print(f"x = {int_to_bin_str(x)} -> S(x) = {int_to_bin_str(S_box[x])} (дец. {S_box[x]})")

    # Вычисление и вывод таблицы линейной аппроксимации (LAT)
    lat = compute_lat(S_box)
    print_lat_table(lat)
    approximations = find_linear_approximations(lat, threshold=4)
    print_linear_approximations(approximations)

    # ================================
    # Поиск ключа в тестовой системе (LAB‑4)
    # ================================
    rounds = 2
    error_counts = {}
    for candidate in range(16):
        errors = 0
        for x in range(12):
            test_val = encrypt_test(x, candidate, rounds, S_box)
            if test_val != partial_table[x]:
                errors += 1
        error_counts[candidate] = errors
    print("\nОшибка тестовой системы для каждого ключа candidate:")
    for candidate in range(16):
        print(f"Ключ {int_to_bin_str(candidate)} : ошибок = {error_counts[candidate]}")

    # Выбор ключа с минимальным числом ошибок
    found_key = min(error_counts, key=error_counts.get)
    print(f"\nНайденный ключ: {int_to_bin_str(found_key)} (дец. {found_key})")

    # Вычисление обратного S-бокса S_inv
    S_inv = [0] * 16
    for x in range(16):
        y = S_box[x]
        S_inv[y] = x

    # Расшифровка зашифрованного сообщения:
    print(f"\nРасшифровка зашифрованного сообщения c(m) = {c_message} (при 12 <= m < 16):")
    recovered_m = None
    for m in range(12, 16):
        trial = encrypt_actual(m, found_key, rounds, S_box)
        print(f"m = {int_to_bin_str(m)} → C(m) = {int_to_bin_str(trial)}")
        if trial == c_message:
            recovered_m = m
    if recovered_m is not None:
        print(
            f"\nЗашифрованное сообщение c(m)={c_message} было получено из m = {int_to_bin_str(recovered_m)} ({recovered_m})")
    else:
        print("\nНе найден вход m, для которого шифрование дало заданное зашифрованное сообщение.")


if __name__ == "__main__":
    main()
