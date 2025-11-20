from functools import reduce

class RingElement:
    """
    Элемент кольца: целое число или многочлен.

    Attributes:
        data (int | list[float]):
            - int -> элемент кольца Z
            - list -> многочлен [a0, a1, ..., an] = a0 + a1*x + ... + an*x^n
        is_polynomial (bool): True, если элемент является многочленом.
    """

    def __init__(self, data):
        """
        Инициализация элемента кольца.

        Args:
            data (int | list[float]): данные элемента.
        """
        self.data = data
        self.is_polynomial = isinstance(data, list)

    def __repr__(self):
        """
        Возвращает строковое представление элемента.
        """
        if self.is_polynomial:
            terms = [f"{c}*x^{i}" if i > 0 else str(c)
                     for i, c in enumerate(self.data) if abs(c) > 1e-12]
            return " + ".join(terms) if terms else "0"
        return str(self.data)


def gcd_integers(numbers: list[int]) -> int:
    """
    Вычисляет НОД списка целых чисел.

    Args:
        numbers (list[int]): список целых чисел.

    Returns:
        int: наибольший общий делитель всех чисел списка.
    """

    def gcd_two(a: int, b: int) -> int:
        """
        Вычисляет НОД двух целых чисел с помощью алгоритма Евклида.
        """
        a, b = abs(a), abs(b)
        while b != 0:
            a, b = b, a % b
        return a

    return reduce(gcd_two, numbers)


def gcd_polynomials(polys: list[RingElement]) -> RingElement:
    """
    Вычисляет НОД списка многочленов в K[x].

    Args:
        polys (list[RingElement]): список элементов- многочленов.

    Returns:
        RingElement: многочлен, порождающий главный идеал.
    """
    d = polys[0].data[:]

    for p in polys[1:]:
        g = p.data[:]
        f = d[:]

        while any(g):
            while len(f) >= len(g):
                deg_diff = len(f) - len(g)
                coef = f[-1] / g[-1]
                for i in range(len(g)):
                    f[i + deg_diff] -= coef * g[i]
                while f and abs(f[-1]) < 1e-12:
                    f.pop()
            f, g = g, f

        lc = f[-1]
        d = [c / lc for c in f]

    return RingElement(d)


def gcd_ring_elements(elements: list[RingElement]) -> RingElement:
    """
    Вычисляет порождающий элемент главного идеала для списка элементов кольца.

    Args:
        elements (list[RingElement]): список элементов кольца (ℤ или K[x]).

    Returns:
        RingElement: элемент, порождающий главный идеал.

    Raises:
        ValueError: если список элементов пуст.
    """
    if not elements:
        raise ValueError("Список элементов пуст")
    if not elements[0].is_polynomial:
        nums = [e.data for e in elements]
        return RingElement(gcd_integers(nums))
    else:
        return gcd_polynomials(elements)


def is_in_ideal(f: RingElement, generators: list[RingElement]) -> bool:
    """
    Проверяет, принадлежит ли элемент f идеалу, порожденному generators.

    Args:
        f (RingElement): проверяемый элемент.
        generators (list[RingElement]): список генераторов идеала.

    Returns:
        bool: True, если f принадлежит идеалу, False иначе.
    """
    d = gcd_ring_elements(generators)

    if not f.is_polynomial:  # ℤ
        return f.data % d.data == 0

    f_coef = f.data[:]
    d_coef = d.data[:]

    while len(f_coef) >= len(d_coef):
        coef = f_coef[-1] / d_coef[-1]
        deg_diff = len(f_coef) - len(d_coef)
        for i in range(len(d_coef)):
            f_coef[i + deg_diff] -= coef * d_coef[i]
        while f_coef and abs(f_coef[-1]) < 1e-12:
            f_coef.pop()

    return all(abs(c) < 1e-12 for c in f_coef)

a = RingElement(12)
b = RingElement(18)
gens_int = [a, b]

d_int = gcd_ring_elements(gens_int)
print("Генератор идеала ℤ:", d_int)  # 6

print("15 ∈ I?", is_in_ideal(RingElement(15), gens_int))  # False
print("6 ∈ I?", is_in_ideal(RingElement(6), gens_int))    # True
print("12 ∈ I?", is_in_ideal(RingElement(12), gens_int))  # True

f1 = RingElement([1, -1])       # x - 1
f2 = RingElement([1, -3, 2])    # 2x^2 - 3x + 1
gens_poly = [f1, f2]

d_poly = gcd_ring_elements(gens_poly)
print("Генератор идеала K[x]:", d_poly)  # x - 1

f_test1 = RingElement([1, -1])       # x - 1
f_test2 = RingElement([1, 0, -1])    # x^2 - 1
f_test3 = RingElement([1, 1])        # x + 1

print("x-1 ∈ I?", is_in_ideal(f_test1, gens_poly))   # True
print("x^2-1 ∈ I?", is_in_ideal(f_test2, gens_poly)) # True
print("x+1 ∈ I?", is_in_ideal(f_test3, gens_poly))   # False