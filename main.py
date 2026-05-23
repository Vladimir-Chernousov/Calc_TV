import sys
import os
import math
import json
import csv
import sqlite3
from datetime import datetime
from typing import List, Any, Tuple

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
    return math.exp(-x * x / 2.0) / math.sqrt(2.0 * math.pi)

def _Phi(z):
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))

def _check_prob01(p, name="p"):
    if not (0.0 <= p <= 1.0):
        raise ValueError(f"{name} должна быть в диапазоне [0, 1].")

def _format_value(value: Any) -> str:
    if isinstance(value, list):
        return ", ".join(_format_value(v) for v in value)
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if value.is_integer():
            return str(int(value))
        return f"{value:.10g}"
    if value is None:
        return ""
    return str(value)

def _collect_hypotheses(*values) -> Tuple[List[float], List[float]]:
    if len(values) % 2 != 0:
        raise ValueError("Некорректное число параметров гипотез.")

    pH = []
    pA_given_H = []

    for idx in range(0, len(values), 2):
        ph = values[idx]
        pa = values[idx + 1]
        h_num = idx // 2 + 1

        if ph is None and pa is None:
            continue
        if ph is None or pa is None:
            raise ValueError(
                f"Для гипотезы H{h_num} нужно заполнить оба поля: P(H{h_num}) и P(A|H{h_num})."
            )

        pH.append(ph)
        pA_given_H.append(pa)

    if not pH:
        raise ValueError("Нужно задать хотя бы одну гипотезу.")

    return pH, pA_given_H


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

# ===== Журнал вычислений =====

def _get_journal_db_path() -> str:
    app_data_dir = QtCore.QStandardPaths.writableLocation(QtCore.QStandardPaths.AppDataLocation)
    if not app_data_dir:
        app_data_dir = os.path.dirname(os.path.abspath(__file__))
    os.makedirs(app_data_dir, exist_ok=True)
    return os.path.join(app_data_dir, "calculations_journal.db")

class ComputationJournal:
    def __init__(self, db_path: str = None):
        self.db_path = db_path or _get_journal_db_path()
        self._init_db()

    def _connect(self):
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        return conn

    def _init_db(self):
        with self._connect() as conn:
            conn.execute("""
                CREATE TABLE IF NOT EXISTS journal (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    created_at TEXT NOT NULL,
                    section_name TEXT NOT NULL,
                    formula_title TEXT NOT NULL,
                    formula_expression TEXT NOT NULL,
                    inputs_json TEXT NOT NULL,
                    result_text TEXT NOT NULL
                )
            """)
            conn.commit()

    def add_entry(self, section_name: str, formula_title: str, formula_expression: str,
                  params_meta: List[dict], values: List[Any], result: Any):
        payload = []
        for meta, value in zip(params_meta, values):
            payload.append({
                "key": meta.get("key", ""),
                "label": meta.get("label", ""),
                "value": _format_value(value),
            })

        with self._connect() as conn:
            conn.execute("""
                INSERT INTO journal (
                    created_at,
                    section_name,
                    formula_title,
                    formula_expression,
                    inputs_json,
                    result_text
                ) VALUES (?, ?, ?, ?, ?, ?)
            """, (
                datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                section_name,
                formula_title,
                formula_expression,
                json.dumps(payload, ensure_ascii=False),
                _format_value(result),
            ))
            conn.commit()

    def fetch_all(self):
        with self._connect() as conn:
            return conn.execute("""
                SELECT *
                FROM journal
                ORDER BY id DESC
            """).fetchall()

    def clear_all(self):
        with self._connect() as conn:
            conn.execute("DELETE FROM journal")
            conn.commit()

    @staticmethod
    def format_inputs_json(inputs_json: str) -> str:
        try:
            items = json.loads(inputs_json)
        except Exception:
            return inputs_json

        parts = []
        for item in items:
            key = item.get("key") or item.get("label") or "param"
            value = item.get("value", "")
            parts.append(f"{key}={value}")
        return "; ".join(parts)

    def export_to_csv(self, file_path: str):
        rows = self.fetch_all()
        with open(file_path, "w", newline="", encoding="utf-8-sig") as f:
            writer = csv.writer(f, delimiter=";")
            writer.writerow([
                "ID",
                "Дата и время",
                "Раздел",
                "Формула",
                "Выражение",
                "Параметры",
                "Результат",
            ])

            for row in rows:
                writer.writerow([
                    row["id"],
                    row["created_at"],
                    row["section_name"],
                    row["formula_title"],
                    row["formula_expression"],
                    self.format_inputs_json(row["inputs_json"]),
                    row["result_text"],
                ])


class JournalDialog(QtWidgets.QDialog):
    def __init__(self, parent, journal: ComputationJournal):
        super().__init__(parent)
        self.journal = journal

        self.setWindowIcon(QtGui.QIcon('src/logo.ico'))
        self.setWindowTitle("Журнал вычислений")
        self.resize(1100, 600)
        self.setModal(True)

        layout = QtWidgets.QVBoxLayout(self)

        title_lbl = QtWidgets.QLabel("<b>Журнал вычислений</b>")
        info_lbl = QtWidgets.QLabel(
            f"Файл базы данных: <span style='color:#555'>{self.journal.db_path}</span>"
        )
        info_lbl.setWordWrap(True)

        self.table = QtWidgets.QTableWidget()
        self.table.setColumnCount(6)
        self.table.setHorizontalHeaderLabels([
            "Дата и время",
            "Раздел",
            "Формула",
            "Выражение",
            "Параметры",
            "Результат",
        ])
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
        self.table.setAlternatingRowColors(True)
        self.table.setWordWrap(True)
        self.table.verticalHeader().setVisible(False)

        header = self.table.horizontalHeader()
        header.setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeToContents)
        header.setSectionResizeMode(1, QtWidgets.QHeaderView.ResizeToContents)
        header.setSectionResizeMode(2, QtWidgets.QHeaderView.ResizeToContents)
        header.setSectionResizeMode(3, QtWidgets.QHeaderView.Stretch)
        header.setSectionResizeMode(4, QtWidgets.QHeaderView.Stretch)
        header.setSectionResizeMode(5, QtWidgets.QHeaderView.ResizeToContents)

        btns = QtWidgets.QHBoxLayout()
        self.refresh_btn = QtWidgets.QPushButton("Обновить")
        self.export_btn = QtWidgets.QPushButton("Экспорт CSV")
        self.clear_btn = QtWidgets.QPushButton("Очистить журнал")
        self.close_btn = QtWidgets.QPushButton("Закрыть")

        btns.addWidget(self.refresh_btn)
        btns.addWidget(self.export_btn)
        btns.addWidget(self.clear_btn)
        btns.addStretch(1)
        btns.addWidget(self.close_btn)

        layout.addWidget(title_lbl)
        layout.addWidget(info_lbl)
        layout.addWidget(self.table)
        layout.addLayout(btns)

        self.refresh_btn.clicked.connect(self.load_data)
        self.export_btn.clicked.connect(self.export_csv)
        self.clear_btn.clicked.connect(self.clear_journal)
        self.close_btn.clicked.connect(self.accept)

        self.load_data()

    def load_data(self):
        rows = self.journal.fetch_all()
        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(rows))

        for row_idx, row in enumerate(rows):
            values = [
                row["created_at"],
                row["section_name"],
                row["formula_title"],
                row["formula_expression"],
                self.journal.format_inputs_json(row["inputs_json"]),
                row["result_text"],
            ]

            for col_idx, value in enumerate(values):
                item = QtWidgets.QTableWidgetItem(str(value))
                item.setTextAlignment(int(QtCore.Qt.AlignLeft | QtCore.Qt.AlignTop))
                self.table.setItem(row_idx, col_idx, item)

        self.table.resizeRowsToContents()
        self.table.setSortingEnabled(True)

    def export_csv(self):
        default_name = f"journal_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Сохранить журнал",
            default_name,
            "CSV Files (*.csv)"
        )
        if not path:
            return

        try:
            self.journal.export_to_csv(path)
            QtWidgets.QMessageBox.information(self, "Готово", "Журнал успешно экспортирован.")
        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Ошибка", f"Не удалось экспортировать журнал:\n{e}")

    def clear_journal(self):
        reply = QtWidgets.QMessageBox.question(
            self,
            "Подтверждение",
            "Очистить весь журнал вычислений?",
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
            QtWidgets.QMessageBox.No
        )
        if reply != QtWidgets.QMessageBox.Yes:
            return

        try:
            self.journal.clear_all()
            self.load_data()
            QtWidgets.QMessageBox.information(self, "Готово", "Журнал очищен.")
        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Ошибка", f"Не удалось очистить журнал:\n{e}")

# ===== Описание формул по вкладкам =====

TABS = {
    "tilesCombinatorics": [
        {
            "title": "Факториал",
            "expression": "n! = 1·2·...·n",
            "params": [{"key": "n", "label": "n — целое, n ≥ 0", "type": "int"}],
            "calc": lambda n: calc_factorial(n)
        },
        {
            "title": "Перестановки (без повтор.)",
            "expression": "P(n) = n!",
            "params": [{"key": "n", "label": "n — целое, n ≥ 0", "type": "int"}],
            "calc": lambda n: calc_permutations_no_rep(n)
        },
        {
            "title": "Перестановки (с повтор.)",
            "expression": "P = n! / (n1! · n2! · ...)",
            "params": [
                {"key": "n", "label": "n — всего элементов", "type": "int"},
                {"key": "повторы", "label": "n_i — повторы (через запятую)", "type": "list_int"},
            ],
            "calc": lambda n, reps: calc_permutations_with_rep(n, reps)
        },
        {
            "title": "Размещения (без повтор.)",
            "expression": "A(n,k) = n! / (n-k)!",
            "params": [
                {"key": "n", "label": "n — целое, n ≥ 0", "type": "int"},
                {"key": "k", "label": "k — целое, 0 ≤ k ≤ n", "type": "int"},
            ],
            "calc": lambda n, k: calc_arrangements_no_rep(n, k)
        },
        {
            "title": "Размещения (с повтор.)",
            "expression": "A'(n,k) = n^k",
            "params": [
                {"key": "n", "label": "n — целое, n ≥ 0", "type": "int"},
                {"key": "k", "label": "k — целое, k ≥ 0", "type": "int"},
            ],
            "calc": lambda n, k: calc_arrangements_with_rep(n, k)
        },
        {
            "title": "Сочетания (без повтор.)",
            "expression": "C(n,k) = n! / (k!(n-k)!)",
            "params": [
                {"key": "n", "label": "n — целое, n ≥ 0", "type": "int"},
                {"key": "k", "label": "k — целое, 0 ≤ k ≤ n", "type": "int"},
            ],
            "calc": lambda n, k: calc_combinations_no_rep(n, k)
        },
        {
            "title": "Сочетания (с повтор.)",
            "expression": "C'(n,k) = C(n+k-1, k)",
            "params": [
                {"key": "n", "label": "n — целое, n ≥ 1", "type": "int"},
                {"key": "k", "label": "k — целое, k ≥ 0", "type": "int"},
            ],
            "calc": lambda n, k: calc_combinations_with_rep(n, k)
        },
    ],
    "tilesProbability": [
        {
            "title": "Классическая вероятность",
            "expression": "P(A) = m / n",
            "params": [
                {"key": "m", "label": "m — благоприятные исходы", "type": "int"},
                {"key": "n", "label": "n — все равновозможные исходы", "type": "int"},
            ],
            "calc": lambda m, n: calc_laplace(m, n)
        },
        {
            "title": "Геометрическая вероятность",
            "expression": "P = мера(благоприятной) / мера(всей)",
            "params": [
                {"key": "S_бл", "label": "Мера благоприятной области", "type": "float"},
                {"key": "S_вс", "label": "Мера всей области (>0)", "type": "float"},
            ],
            "calc": lambda s1, s2: calc_geometric_prob(s1, s2)
        },
        {
            "title": "Статистическая вероятность",
            "expression": "P ≈ m / N (при большом N)",
            "params": [
                {"key": "m", "label": "m — число успехов", "type": "int"},
                {"key": "N", "label": "N — число испытаний", "type": "int"},
            ],
            "calc": lambda m, N: calc_statistical_prob(m, N)
        },
        {
            "title": "Теорема сложения \nдля несовместных событий",
            "expression": "P(A∪B) = P(A) + P(B)",
            "params": [
                {"key": "P(A)", "label": "P(A) в [0,1]", "type": "float"},
                {"key": "P(B)", "label": "P(B) в [0,1]", "type": "float"},
            ],
            "calc": lambda pA, pB: calc_addition_disjoint(pA, pB)
        },
        {
            "title": "Теорема сложения \nдля совместных событий",
            "expression": "P(A∪B) = P(A)+P(B)-P(A∩B)",
            "params": [
                {"key": "P(A)", "label": "P(A) в [0,1]", "type": "float"},
                {"key": "P(B)", "label": "P(B) в [0,1]", "type": "float"},
                {"key": "P(A∩B)", "label": "P(A∩B) в [0,1]", "type": "float"},
            ],
            "calc": lambda pA, pB, pAB: calc_addition_general(pA, pB, pAB)
        },
        {
            "title": "Теорема умножения \nдля независимых событий",
            "expression": "P(A∩B) = P(A)·P(B)",
            "params": [
                {"key": "P(A)", "label": "P(A) в [0,1]", "type": "float"},
                {"key": "P(B)", "label": "P(B) в [0,1]", "type": "float"},
            ],
            "calc": lambda pA, pB: calc_mult_independent(pA, pB)
        },
        {
            "title": "Теорема умножения \nдля зависимых событий",
            "expression": "P(A∩B) = P(A|B)·P(B)",
            "params": [
                {"key": "P(A|B)", "label": "P(A|B) в [0,1]", "type": "float"},
                {"key": "P(B)", "label": "P(B) в [0,1]", "type": "float"},
            ],
            "calc": lambda pAgB, pB: calc_mult_general(pAgB, pB)
        },
        {
            "title": "Полная вероятность (до 3 гип.)",
            "expression": "P(A)=Σ P(H_i)·P(A|H_i)",
            "params": [
                {"key": "P(H1)", "label": "P(H1) в [0,1]", "type": "float_optional"},
                {"key": "P(A|H1)", "label": "P(A|H1) в [0,1]", "type": "float_optional"},
                {"key": "P(H2)", "label": "P(H2) в [0,1]", "type": "float_optional"},
                {"key": "P(A|H2)", "label": "P(A|H2) в [0,1]", "type": "float_optional"},
                {"key": "P(H3)", "label": "P(H3) в [0,1]", "type": "float_optional"},
                {"key": "P(A|H3)", "label": "P(A|H3) в [0,1]", "type": "float_optional"},
            ],
            "calc": lambda ph1, pa1, ph2, pa2, ph3, pa3: calc_total_probability(
                *_collect_hypotheses(ph1, pa1, ph2, pa2, ph3, pa3)
            )
        },
        {
            "title": "Формула Байеса (до 3 гип.)",
            "expression": "P(H_i|A)=P(H_i)P(A|H_i)/Σ P(H_j)P(A|H_j)",
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
                *_collect_hypotheses(ph1, pa1, ph2, pa2, ph3, pa3)
            )
        },
    ],
    "tilesRepeated": [
        {
            "title": "Формула Бернулли",
            "expression": "P(X=k)=C(n,k)·p^k·(1-p)^(n-k)",
            "params": [
                {"key": "n", "label": "n — число испытаний", "type": "int"},
                {"key": "k", "label": "k — число успехов", "type": "int"},
                {"key": "p", "label": "p — вероятность успеха (0..1)", "type": "float"},
            ],
            "calc": lambda n, k, p: calc_bernoulli(n, k, p)
        },
        {
            "title": "Локальная формула Лапласа",
            "expression": "≈ φ(z)/√(npq), z=(k-np)/√(npq)",
            "params": [
                {"key": "n", "label": "n — число испытаний", "type": "int"},
                {"key": "p", "label": "p — вероятность успеха (0..1)", "type": "float"},
                {"key": "k", "label": "k — число успехов", "type": "int"},
            ],
            "calc": lambda n, p, k: calc_laplace_local(n, p, k)
        },
        {
            "title": "Интегральная формула Лапласа",
            "expression": "≈ Φ(z2)-Φ(z1) с непрерывной поправкой",
            "params": [
                {"key": "n", "label": "n — число испытаний", "type": "int"},
                {"key": "p", "label": "p — вероятность успеха (0..1)", "type": "float"},
                {"key": "k1", "label": "k1 — нижняя граница", "type": "int"},
                {"key": "k2", "label": "k2 — верхняя граница", "type": "int"},
            ],
            "calc": lambda n, p, k1, k2: calc_laplace_integral(n, p, k1, k2)
        },
        {
            "title": "Формула Пуассона",
            "expression": "P(X=k)=e^(-λ)·λ^k/k!",
            "params": [
                {"key": "λ", "label": "λ — среднее (≥0)", "type": "float"},
                {"key": "k", "label": "k — целое (≥0)", "type": "int"},
            ],
            "calc": lambda lmbd, k: calc_poisson(lmbd, k)
        },
    ],
}


# ===== Диалог формулы =====

class FormulaDialog(QtWidgets.QDialog):
    def __init__(self, parent, section_name, title, expression, params, calc_fn, journal=None):
        super().__init__(parent)
        self.setWindowIcon(QtGui.QIcon('src/logo.ico'))
        self.setWindowTitle(title)

        self.section_name = section_name
        self.formula_title = title
        self.formula_expression = expression
        self.calc_fn = calc_fn
        self.params_meta = params
        self.journal = journal

        self.setModal(True)
        self.setMinimumWidth(500)

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
                edit.setToolTip("Список целых через запятую, например: 1,2,3")
            self.edits.append(edit)
            form.addRow(p["label"] + ":", edit)

        self.result_lbl = QtWidgets.QLabel("")
        self.result_lbl.setWordWrap(True)

        self.saved_lbl = QtWidgets.QLabel("")
        self.saved_lbl.setWordWrap(True)

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
        main_layout.addWidget(self.saved_lbl)
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
            out = _format_value(result)

            self.result_lbl.setText(f"<b>Результат: {out}</b>")
            self.saved_lbl.setText("")

            if self.journal is not None:
                try:
                    self.journal.add_entry(
                        section_name=self.section_name,
                        formula_title=self.formula_title,
                        formula_expression=self.formula_expression,
                        params_meta=self.params_meta,
                        values=values,
                        result=out,
                    )
                    self.saved_lbl.setText(
                        "<span style='color:#2E7D32;'>Результат сохранён в журнал вычислений.</span>"
                    )
                except Exception as journal_error:
                    self.saved_lbl.setText(
                        "<span style='color:#C62828;'>Результат посчитан, но не удалось сохранить запись.</span>"
                    )
                    QtWidgets.QMessageBox.warning(
                        self,
                        "Предупреждение",
                        f"Результат вычислен, но журнал не был обновлён:\n{journal_error}"
                    )

        except Exception as e:
            self.saved_lbl.setText("")
            QtWidgets.QMessageBox.critical(self, "Ошибка ввода", str(e))


# ===== Главное окно =====

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        self.journal = ComputationJournal()

        self.tiles_combinatorics = self.ui.tilesCombinatorics
        self.tiles_probability = self.ui.tilesProbability
        self.tiles_repeated = self.ui.tilesRepeated

        self._setup_style()
        self._setup_menu()
        self._populate_tiles(self.tiles_combinatorics, TABS["tilesCombinatorics"], "Комбинаторика")
        self._populate_tiles(self.tiles_probability, TABS["tilesProbability"], "Теория вероятностей")
        self._populate_tiles(self.tiles_repeated, TABS["tilesRepeated"], "Независимые повторные испытания")

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

    def _setup_menu(self):
        journal_menu = self.menuBar().addMenu("Журнал")

        self.open_journal_action = QtWidgets.QAction("Открыть журнал", self)
        self.open_journal_action.setShortcut("Ctrl+J")
        self.open_journal_action.triggered.connect(self._show_journal)

        self.export_journal_action = QtWidgets.QAction("Экспорт в CSV", self)
        self.export_journal_action.triggered.connect(self._export_journal)

        self.clear_journal_action = QtWidgets.QAction("Очистить журнал", self)
        self.clear_journal_action.triggered.connect(self._clear_journal)

        journal_menu.addAction(self.open_journal_action)
        journal_menu.addAction(self.export_journal_action)
        journal_menu.addSeparator()
        journal_menu.addAction(self.clear_journal_action)

        self.statusBar().showMessage("Готово. Ctrl+J — открыть журнал вычислений.", 5000)

    def _populate_tiles(self, container_widget, formulas, section_name):
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
            btn.clicked.connect(lambda _, f=f, sec=section_name: self._open_formula(sec, f))
            layout.addWidget(btn, r, c)
            c += 1
            if c >= cols:
                c = 0
                r += 1

        for i in range(cols):
            layout.setColumnStretch(i, 1)
        for i in range(r + 1):
            layout.setRowStretch(i, 1)

    def _open_formula(self, section_name, fdef):
        dlg = FormulaDialog(
            self,
            section_name=section_name,
            title=fdef["title"],
            expression=fdef["expression"],
            params=fdef["params"],
            calc_fn=fdef["calc"],
            journal=self.journal,
        )
        dlg.exec_()

    def _show_journal(self):
        dlg = JournalDialog(self, self.journal)
        dlg.exec_()

    def _export_journal(self):
        default_name = f"journal_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Сохранить журнал",
            default_name,
            "CSV Files (*.csv)"
        )
        if not path:
            return

        try:
            self.journal.export_to_csv(path)
            QtWidgets.QMessageBox.information(self, "Готово", "Журнал успешно экспортирован.")
        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Ошибка", f"Не удалось экспортировать журнал:\n{e}")

    def _clear_journal(self):
        reply = QtWidgets.QMessageBox.question(
            self,
            "Подтверждение",
            "Очистить весь журнал вычислений?",
            QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
            QtWidgets.QMessageBox.No
        )
        if reply != QtWidgets.QMessageBox.Yes:
            return

        try:
            self.journal.clear_all()
            QtWidgets.QMessageBox.information(self, "Готово", "Журнал очищен.")
        except Exception as e:
            QtWidgets.QMessageBox.critical(self, "Ошибка", f"Не удалось очистить журнал:\n{e}")


# ===== Точка входа =====

def main():
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)

    app = QtWidgets.QApplication(sys.argv)
    app.setApplicationName("Калькулятор теории вероятностей")
    app.setOrganizationName("ProbabilityApp")

    win = MainWindow()
    win.resize(980, 640)
    win.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()