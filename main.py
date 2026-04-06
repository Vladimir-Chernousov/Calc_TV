import sys
import math
from typing import List
from PyQt5 import QtWidgets, QtCore, QtGui
from ui_design import Ui_MainWindow

# ===== Вспомогательные функции =====

def _ensure_int(x, name):
    if int(x) != x:
        raise ValueError(f"{name} должен быть целым числом.")
    return int(x)

def _comb(n, k):
    n = _ensure_int(n, "n")
    k = _ensure_int(k, "k")
    if n < 0 or k < 0 or k > n:
        return 0
    k = min(k, n - k)
    numer = 1
    denom = 1
    for i in range(1, k + 1):
        numer *= (n - k + i)
        denom *= i
    return numer // denom

def _phi_standard_normal(x):
    return math.exp(-x*x/2.0) / math.sqrt(2.0 * math.pi)

def _Phi(z):
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))

def _check_prob01(p, name="p"):
    if not (0.0 <= p <= 1.0):
        raise ValueError(f"{name} должна быть в диапазоне [0, 1].")

# ===== Комбинаторика =====

def calc_factorial(n):
    n = _ensure_int(n, "n")
    if n < 0:
        raise ValueError("n должно быть >= 0.")
    return math.factorial(n)

def calc_permutations_no_rep(n):
    return calc_factorial(n)

def calc_permutations_with_rep(n, repeats: List[int]):
    n = _ensure_int(n, "n")
    if n < 0:
        raise ValueError("n должно быть >= 0.")
    if sum(repeats) != n:
        raise ValueError("Сумма повторов должна равняться n.")
    res = math.factorial(n)
    for r in repeats:
        if r < 0:
            raise ValueError("Повторы должны быть неотрицательны.")
        res //= math.factorial(int(r))
    return res

def calc_arrangements_no_rep(n, k):
    n = _ensure_int(n, "n")
    k = _ensure_int(k, "k")
    if not (0 <= k <= n):
        raise ValueError("Должно выполняться 0 <= k <= n.")
    return math.factorial(n) // math.factorial(n - k)

def calc_arrangements_with_rep(n, k):
    n = _ensure_int(n, "n")
    k = _ensure_int(k, "k")
    if n < 0 or k < 0:
        raise ValueError("n, k должны быть >= 0.")
    return n ** k

def calc_combinations_no_rep(n, k):
    return _comb(n, k)

def calc_combinations_with_rep(n, k):
    n = _ensure_int(n, "n")
    k = _ensure_int(k, "k")
    if n <= 0 or k < 0:
        raise ValueError("n должно быть >= 1, k >= 0.")
    return _comb(n + k - 1, k)

# ===== Общая теория вероятностей =====

def calc_laplace(m, n):
    m = _ensure_int(m, "m")
    n = _ensure_int(n, "n")
    if n == 0:
        raise ValueError("n должно быть > 0 (деление на ноль).")
    if not (0 <= m <= n):
        raise ValueError("0 <= m <= n.")
    return m / n

def calc_geometric_prob(s_fav, s_all):
    if s_all <= 0:
        raise ValueError("Мера всей области должна быть > 0.")
    if not (0 <= s_fav <= s_all):
        raise ValueError("Мера благоприятной области должна быть в [0, мера всей области].")
    return s_fav / s_all

def calc_statistical_prob(m, N):
    m = _ensure_int(m, "m")
    N = _ensure_int(N, "N")
    if N <= 0:
        raise ValueError("N должно быть > 0.")
    if not (0 <= m <= N):
        raise ValueError("0 <= m <= N.")
    return m / N

def calc_addition_disjoint(pA, pB):
    _check_prob01(pA, "P(A)")
    _check_prob01(pB, "P(B)")
    s = pA + pB
    if s > 1.0 + 1e-12:
        raise ValueError("Для несовместимых событий сумма вероятностей не должна превышать 1.")
    return s

def calc_addition_general(pA, pB, pAintB):
    _check_prob01(pA, "P(A)")
    _check_prob01(pB, "P(B)")
    _check_prob01(pAintB, "P(A∩B)")
    if pAintB - min(pA, pB) > 1e-12:
        raise ValueError("P(A∩B) не может превышать min(P(A), P(B)).")
    res = pA + pB - pAintB
    if not (-1e-12 <= res <= 1 + 1e-12):
        raise ValueError("Результат вышел из [0,1], проверьте ввод.")
    return max(0.0, min(1.0, res))

def calc_mult_independent(pA, pB):
    _check_prob01(pA, "P(A)")
    _check_prob01(pB, "P(B)")
    return pA * pB

def calc_mult_general(pA_given_B, pB):
    _check_prob01(pA_given_B, "P(A|B)")
    _check_prob01(pB, "P(B)")
    return pA_given_B * pB

def calc_total_probability(pH: List[float], pA_given_H: List[float]):
    if len(pH) != len(pA_given_H):
        raise ValueError("Размерности списков не совпадают.")
    if len(pH) == 0:
        raise ValueError("Нужно задать хотя бы одну гипотезу.")
    for i, (ph, pa) in enumerate(zip(pH, pA_given_H), start=1):
        _check_prob01(ph, f"P(H{i})")
        _check_prob01(pa, f"P(A|H{i})")
    return sum(ph * pa for ph, pa in zip(pH, pA_given_H))

def calc_bayes_idx(i, pH: List[float], pA_given_H: List[float]):
    if not (1 <= i <= len(pH)):
        raise ValueError("Индекс гипотезы i вне диапазона.")
    denom = calc_total_probability(pH, pA_given_H)
    if denom <= 0:
        raise ValueError("Полная вероятность события A равна 0 (деление на ноль).")
    numer = pH[i - 1] * pA_given_H[i - 1]
    return numer / denom

# ===== Независимые повторные испытания =====

def calc_bernoulli(n, k, p):
    n = _ensure_int(n, "n")
    k = _ensure_int(k, "k")
    _check_prob01(p, "p")
    if not (0 <= k <= n):
        raise ValueError("0 <= k <= n.")
    return _comb(n, k) * (p ** k) * ((1 - p) ** (n - k))

def calc_poisson(lmbd, k):
    k = _ensure_int(k, "k")
    if lmbd < 0:
        raise ValueError("λ должно быть >= 0.")
    if k < 0:
        raise ValueError("k должно быть >= 0.")
    return math.exp(-lmbd) * (lmbd ** k) / math.factorial(k)

def calc_laplace_local(n, p, k):
    n = _ensure_int(n, "n")
    k = _ensure_int(k, "k")
    _check_prob01(p, "p")
    q = 1 - p
    if n <= 0:
        raise ValueError("n должно быть > 0.")
    var = n * p * q
    if var == 0:
        raise ValueError("Вариация равна 0 (p=0 или p=1) — приближение неприменимо.")
    mu = n * p
    z = (k - mu) / math.sqrt(var)
    return _phi_standard_normal(z) / math.sqrt(var)

def calc_laplace_integral(n, p, k1, k2):
    n = _ensure_int(n, "n")
    k1 = _ensure_int(k1, "k1")
    k2 = _ensure_int(k2, "k2")
    _check_prob01(p, "p")
    if k2 < k1:
        raise ValueError("k2 должен быть >= k1.")
    q = 1 - p
    var = n * p * q
    if var == 0:
        raise ValueError("Вариация равна 0 (p=0 или p=1) — приближение неприменимо.")
    mu = n * p
    s = math.sqrt(var)
    z2 = (k2 + 0.5 - mu) / s
    z1 = (k1 - 0.5 - mu) / s
    return _Phi(z2) - _Phi(z1)

# ===== Описание формул по вкладкам =====

TABS = {
    "tilesCombinatorics": [
        {"title": "Факториал", "expression": "n! = 1·2·...·n",
         "params": [{"key": "n", "label": "n — целое, n ≥ 0", "type": "int"}],
         "calc": lambda n: calc_factorial(n)},
        {"title": "Перестановки (без повтор.)", "expression": "P(n) = n!",
         "params": [{"key": "n", "label": "n — целое, n ≥ 0", "type": "int"}],
         "calc": lambda n: calc_permutations_no_rep(n)},
        {"title": "Перестановки (с повтор.)", "expression": "P = n! / (n1! · n2! · ...)",
         "params": [
             {"key": "n", "label": "n — всего элементов", "type": "int"},
             {"key": "повторы", "label": "n_i — повторы (через запятую)", "type": "list_int"},
         ],
         "calc": lambda n, reps: calc_permutations_with_rep(n, reps)},
        {"title": "Размещения (без повтор.)", "expression": "A(n,k) = n! / (n-k)!",
         "params": [
             {"key": "n", "label": "n — целое, n ≥ 0", "type": "int"},
             {"key": "k", "label": "k — целое, 0 ≤ k ≤ n", "type": "int"},
         ],
         "calc": lambda n, k: calc_arrangements_no_rep(n, k)},
        {"title": "Размещения (с повтор.)", "expression": "A'(n,k) = n^k",
         "params": [
             {"key": "n", "label": "n — целое, n ≥ 0", "type": "int"},
             {"key": "k", "label": "k — целое, k ≥ 0", "type": "int"},
         ],
         "calc": lambda n, k: calc_arrangements_with_rep(n, k)},
        {"title": "Сочетания (без повтор.)", "expression": "C(n,k) = n! / (k!(n-k)!)",
         "params": [
             {"key": "n", "label": "n — целое, n ≥ 0", "type": "int"},
             {"key": "k", "label": "k — целое, 0 ≤ k ≤ n", "type": "int"},
         ],
         "calc": lambda n, k: calc_combinations_no_rep(n, k)},
        {"title": "Сочетания (с повтор.)", "expression": "C'(n,k) = C(n+k-1, k)",
         "params": [
             {"key": "n", "label": "n — целое, n ≥ 1", "type": "int"},
             {"key": "k", "label": "k — целое, k ≥ 0", "type": "int"},
         ],
         "calc": lambda n, k: calc_combinations_with_rep(n, k)},
    ],
    "tilesProbability": [
        {"title": "Классическая вероятность", "expression": "P(A) = m / n",
         "params": [
             {"key": "m", "label": "m — благоприятные исходы", "type": "int"},
             {"key": "n", "label": "n — все равновозможные исходы", "type": "int"},
         ],
         "calc": lambda m, n: calc_laplace(m, n)},
        {"title": "Геометрическая вероятность", "expression": "P = мера(благоприятной) / мера(всей)",
         "params": [
             {"key": "S_бл", "label": "Мера благоприятной области", "type": "float"},
             {"key": "S_вс", "label": "Мера всей области (>0)", "type": "float"},
         ],
         "calc": lambda s1, s2: calc_geometric_prob(s1, s2)},
        {"title": "Статистическая вероятность", "expression": "P ≈ m / N (при большом N)",
         "params": [
             {"key": "m", "label": "m — число успехов", "type": "int"},
             {"key": "N", "label": "N — число испытаний", "type": "int"},
         ],
         "calc": lambda m, N: calc_statistical_prob(m, N)},
        {"title": "Сложение (несовместимые)", "expression": "P(A∪B) = P(A) + P(B)",
         "params": [
             {"key": "P(A)", "label": "P(A) в [0,1]", "type": "float"},
             {"key": "P(B)", "label": "P(B) в [0,1]", "type": "float"},
         ],
         "calc": lambda pA, pB: calc_addition_disjoint(pA, pB)},
        {"title": "Сложение (общий случай)", "expression": "P(A∪B) = P(A)+P(B)-P(A∩B)",
         "params": [
             {"key": "P(A)", "label": "P(A) в [0,1]", "type": "float"},
             {"key": "P(B)", "label": "P(B) в [0,1]", "type": "float"},
             {"key": "P(A∩B)", "label": "P(A∩B) в [0,1]", "type": "float"},
         ],
         "calc": lambda pA, pB, pAB: calc_addition_general(pA, pB, pAB)},
        {"title": "Умножение (независимые)", "expression": "P(A∩B) = P(A)·P(B)",
         "params": [
             {"key": "P(A)", "label": "P(A) в [0,1]", "type": "float"},
             {"key": "P(B)", "label": "P(B) в [0,1]", "type": "float"},
         ],
         "calc": lambda pA, pB: calc_mult_independent(pA, pB)},
        {"title": "Умножение (общий случай)", "expression": "P(A∩B) = P(A|B)·P(B)",
         "params": [
             {"key": "P(A|B)", "label": "P(A|B) в [0,1]", "type": "float"},
             {"key": "P(B)", "label": "P(B) в [0,1]", "type": "float"},
         ],
         "calc": lambda pAgB, pB: calc_mult_general(pAgB, pB)},
        {"title": "Полная вероятность (до 3 гип.)", "expression": "P(A)=Σ P(H_i)·P(A|H_i)",
         "params": [
             {"key": "P(H1)", "label": "P(H1) в [0,1]", "type": "float_optional"},
             {"key": "P(A|H1)", "label": "P(A|H1) в [0,1]", "type": "float_optional"},
             {"key": "P(H2)", "label": "P(H2) в [0,1]", "type": "float_optional"},
             {"key": "P(A|H2)", "label": "P(A|H2) в [0,1]", "type": "float_optional"},
             {"key": "P(H3)", "label": "P(H3) в [0,1]", "type": "float_optional"},
             {"key": "P(A|H3)", "label": "P(A|H3) в [0,1]", "type": "float_optional"},
         ],
         "calc": lambda ph1, pa1, ph2, pa2, ph3, pa3: calc_total_probability(
             [ph for ph in [ph1, ph2, ph3] if ph is not None],
             [pa for pa in [pa1, pa2, pa3] if pa is not None],
         )},
        {"title": "Формула Байеса (до 3 гип.)", "expression": "P(H_i|A)=P(H_i)P(A|H_i)/Σ P(H_j)P(A|H_j)",
         "params": [
             {"key": "i", "label": "Индекс гипотезы i (1..3)", "type": "int"},
             {"key": "P(H1)", "label": "P(H1) в [0,1]", "type": "float_optional"},
             {"key": "P(A|H1)", "label": "P(A|H1) в [0,1]", "type": "float_optional"},
             {"key": "P(H2)", "label": "P(H2) в [0,1]", "type": "float_optional"},
             {"key": "P(A|H2)", "label": "P(A|H2) в [0,1]", "type": "float_optional"},
             {"key": "P(H3)", "label": "P(H3) в [0,1]", "type": "float_optional"},
             {"key": "P(A|H3)", "label": "P(A|H3) в [0,1]", "type": "float_optional"},
         ],
         "calc": lambda i, ph1, pa1, ph2, pa2, ph3, pa3: calc_bayes_idx(
             i,
             [ph for ph in [ph1, ph2, ph3] if ph is not None],
             [pa for pa in [pa1, pa2, pa3] if pa is not None],
         )},
    ],
    "tilesRepeated": [
        {"title": "Формула Бернулли", "expression": "P(X=k)=C(n,k)·p^k·(1-p)^{n-k}",
         "params": [
             {"key": "n", "label": "n — число испытаний", "type": "int"},
             {"key": "k", "label": "k — число успехов", "type": "int"},
             {"key": "p", "label": "p — вероятность успеха (0..1)", "type": "float"},
         ],
         "calc": lambda n, k, p: calc_bernoulli(n, k, p)},
        {"title": "Лаплас (локальная)", "expression": "≈ φ(z)/√(npq), z=(k-np)/√(npq)",
         "params": [
             {"key": "n", "label": "n — число испытаний", "type": "int"},
             {"key": "p", "label": "p — вероятность успеха (0..1)", "type": "float"},
             {"key": "k", "label": "k — число успехов", "type": "int"},
         ],
         "calc": lambda n, p, k: calc_laplace_local(n, p, k)},
        {"title": "Лаплас (интегральная)", "expression": "≈ Φ(z2)-Φ(z1) с непрерывн. поправкой",
         "params": [
             {"key": "n", "label": "n — число испытаний", "type": "int"},
             {"key": "p", "label": "p — вероятность успеха (0..1)", "type": "float"},
             {"key": "k1", "label": "k1 — нижняя граница", "type": "int"},
             {"key": "k2", "label": "k2 — верхняя граница", "type": "int"},
         ],
         "calc": lambda n, p, k1, k2: calc_laplace_integral(n, p, k1, k2)},
        {"title": "Пуассон", "expression": "P(X=k)=e^{-λ}·λ^k/k!",
         "params": [
             {"key": "λ", "label": "λ — среднее (≥0)", "type": "float"},
             {"key": "k", "label": "k — целое (≥0)", "type": "int"},
         ],
         "calc": lambda lmbd, k: calc_poisson(lmbd, k)},
    ],
}

# ===== Диалог формулы =====

class FormulaDialog(QtWidgets.QDialog):
    def __init__(self, parent, title, expression, params, calc_fn):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.calc_fn = calc_fn
        self.params_meta = params

        self.setModal(True)
        self.setMinimumWidth(480)

        main_layout = QtWidgets.QVBoxLayout(self)
        title_lbl = QtWidgets.QLabel(f"<b>{title}</b>")
        expr_lbl = QtWidgets.QLabel(f"Формула: {expression}")
        expr_lbl.setWordWrap(True)

        form = QtWidgets.QFormLayout()
        self.edits = []

        for p in params:
            edit = QtWidgets.QLineEdit()
            edit.setPlaceholderText(p["key"])
            edit.setClearButtonEnabled(True)
            if p.get("type") == "list_int":
                edit.setToolTip("Список целых через запятую, например: 2,1,3")
            self.edits.append(edit)
            form.addRow(p["label"] + ":", edit)

        self.result_lbl = QtWidgets.QLabel("")
        self.result_lbl.setWordWrap(True)

        btns = QtWidgets.QHBoxLayout()
        self.calc_btn = QtWidgets.QPushButton("Рассчитать")
        self.close_btn = QtWidgets.QPushButton("Закрыть")
        btns.addStretch(1)
        btns.addWidget(self.calc_btn)
        btns.addWidget(self.close_btn)

        main_layout.addWidget(title_lbl)
        main_layout.addWidget(expr_lbl)
        main_layout.addLayout(form)
        main_layout.addWidget(self.result_lbl)
        main_layout.addLayout(btns)

        self.calc_btn.clicked.connect(self._on_calc)
        self.close_btn.clicked.connect(self.accept)

    @staticmethod
    def _to_float(text, allow_empty=False):
        t = text.strip().replace(",", ".")
        if t == "":
            if allow_empty:
                return None
            raise ValueError("Пустое поле.")
        return float(t)

    @staticmethod
    def _to_int(text):
        t = text.strip().replace(",", ".")
        if t == "":
            raise ValueError("Пустое поле.")
        val = float(t)
        if int(val) != val:
            raise ValueError("Ожидалось целое число.")
        return int(val)

    @staticmethod
    def _to_list_int(text):
        t = text.strip()
        if t == "":
            raise ValueError("Пустой список повторов.")
        raw = [x for x in t.replace(";", ",").replace(" ", ",").split(",") if x.strip() != ""]
        try:
            arr = [int(float(x)) for x in raw]
        except Exception:
            raise ValueError("Список повторов должен содержать целые значения.")
        return arr

    def _parse_inputs(self):
        values = []
        for meta, edit in zip(self.params_meta, self.edits):
            typ = meta.get("type", "float")
            txt = edit.text()
            if typ == "int":
                values.append(self._to_int(txt))
            elif typ == "float":
                values.append(self._to_float(txt, allow_empty=False))
            elif typ == "float_optional":
                values.append(self._to_float(txt, allow_empty=True))
            elif typ == "list_int":
                values.append(self._to_list_int(txt))
            else:
                values.append(self._to_float(txt))
        return values

    def _on_calc(self):
        try:
            values = self._parse_inputs()
            result = self.calc_fn(*values)
            if isinstance(result, int) or (isinstance(result, float) and result.is_integer()):
                out = str(int(round(result)))
            else:
                out = f"{result:.10g}"
            self.result_lbl.setText(f"<b>Результат: {out}</b>")
        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Ошибка ввода", str(e))

# ===== Главное окно =====

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        # Контейнеры плиток на вкладках
        self.tiles_combinatorics = self.ui.tilesCombinatorics
        self.tiles_probability   = self.ui.tilesProbability
        self.tiles_repeated      = self.ui.tilesRepeated

        self._setup_style()
        self._populate_tiles(self.tiles_combinatorics, TABS["tilesCombinatorics"])
        self._populate_tiles(self.tiles_probability,   TABS["tilesProbability"])
        self._populate_tiles(self.tiles_repeated,      TABS["tilesRepeated"])

    def _setup_style(self):
        self.setStyleSheet("""
            QWidget#tilesCombinatorics, QWidget#tilesProbability, QWidget#tilesRepeated {
                background: #F5F6F8;
            }
            QPushButton[cssClass="formula-tile"] {
                background-color: #FFFFFF;
                border: 1px solid #D9D9D9;
                border-radius: 10px;
                padding: 18px 14px;
                text-align: center;
                font-size: 15px;
                min-width: 200px;
                min-height: 84px;
            }
            QPushButton[cssClass="formula-tile"]:hover {
                background-color: #E8F4FF;
                border-color: #80BFFF;
            }
        """)
        self.ui.headerLabel.setStyleSheet("color: #222;")

    def _populate_tiles(self, container_widget, formulas):
        layout = container_widget.layout()
        if layout is None:
            layout = QtWidgets.QGridLayout(container_widget)
            container_widget.setLayout(layout)

        while layout.count():
            item = layout.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()

        cols = 3
        r = c = 0
        for f in formulas:
            btn = QtWidgets.QPushButton(f["title"], container_widget)
            btn.setToolTip(f'Формула: {f["expression"]}')
            btn.setProperty("cssClass", "formula-tile")
            btn.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
            btn.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
            btn.clicked.connect(lambda _, f=f: self._open_formula(f))
            layout.addWidget(btn, r, c)
            c += 1
            if c >= cols:
                c = 0
                r += 1

        for i in range(cols):
            layout.setColumnStretch(i, 1)
        for i in range(r + 1):
            layout.setRowStretch(i, 1)

    def _open_formula(self, fdef):
        dlg = FormulaDialog(
            self,
            title=fdef["title"],
            expression=fdef["expression"],
            params=fdef["params"],
            calc_fn=fdef["calc"],
        )
        dlg.exec_()

# ===== Точка входа =====

def main():
    # Атрибуты HiDPI — обязательно до создания QApplication
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)

    app = QtWidgets.QApplication(sys.argv)
    win = MainWindow()
    win.resize(980, 640)
    win.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()