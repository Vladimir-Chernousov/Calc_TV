from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1280, 640)
        MainWindow.setWindowTitle("Калькулятор теории вероятностей")
        MainWindow.setWindowIcon(QtGui.QIcon('src/logo.ico'))

        self.centralwidget = QtWidgets.QWidget(MainWindow)
        MainWindow.setCentralWidget(self.centralwidget)

        self.verticalLayout = QtWidgets.QVBoxLayout(self.centralwidget)
        self.verticalLayout.setContentsMargins(12, 12, 12, 12)
        self.verticalLayout.setSpacing(8)

        # Заголовок
        self.headerLabel = QtWidgets.QLabel(self.centralwidget)
        self.headerLabel.setObjectName("headerLabel")
        font = QtGui.QFont()
        font.setPointSize(16)
        font.setBold(True)
        self.headerLabel.setFont(font)
        self.headerLabel.setAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
        self.headerLabel.setText("Калькулятор теории вероятностей")
        self.verticalLayout.addWidget(self.headerLabel)

        # Вкладки
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setObjectName("tabWidget")
        self.verticalLayout.addWidget(self.tabWidget)

        # ---------- Вкладка 1: Комбинаторика ----------
        self.tabCombinatorics = QtWidgets.QWidget()
        self.tabCombinatorics.setObjectName("tabCombinatorics")
        self.tabWidget.addTab(self.tabCombinatorics, "Комбинаторика")

        vboxComb = QtWidgets.QVBoxLayout(self.tabCombinatorics)
        vboxComb.setContentsMargins(0, 0, 0, 0)
        vboxComb.setSpacing(0)

        self.scrollCombinatorics = QtWidgets.QScrollArea(self.tabCombinatorics)
        self.scrollCombinatorics.setObjectName("scrollCombinatorics")
        self.scrollCombinatorics.setWidgetResizable(True)
        vboxComb.addWidget(self.scrollCombinatorics)

        self.tilesCombinatorics = QtWidgets.QWidget()
        self.tilesCombinatorics.setObjectName("tilesCombinatorics")
        self.scrollCombinatorics.setWidget(self.tilesCombinatorics)

        self.gridCombinatorics = QtWidgets.QGridLayout(self.tilesCombinatorics)
        self.gridCombinatorics.setContentsMargins(8, 8, 8, 8)
        self.gridCombinatorics.setHorizontalSpacing(12)
        self.gridCombinatorics.setVerticalSpacing(12)

        # ---------- Вкладка 2: Теория Вероятности ----------
        self.tabProbability = QtWidgets.QWidget()
        self.tabProbability.setObjectName("tabProbability")
        self.tabWidget.addTab(self.tabProbability, "Теория Вероятностей")

        vboxProb = QtWidgets.QVBoxLayout(self.tabProbability)
        vboxProb.setContentsMargins(0, 0, 0, 0)
        vboxProb.setSpacing(0)

        self.scrollProbability = QtWidgets.QScrollArea(self.tabProbability)
        self.scrollProbability.setObjectName("scrollProbability")
        self.scrollProbability.setWidgetResizable(True)
        vboxProb.addWidget(self.scrollProbability)

        self.tilesProbability = QtWidgets.QWidget()
        self.tilesProbability.setObjectName("tilesProbability")
        self.scrollProbability.setWidget(self.tilesProbability)

        self.gridProbability = QtWidgets.QGridLayout(self.tilesProbability)
        self.gridProbability.setContentsMargins(8, 8, 8, 8)
        self.gridProbability.setHorizontalSpacing(12)
        self.gridProbability.setVerticalSpacing(12)

        # ---------- Вкладка 3: Независимые повторные испытания ----------
        self.tabRepeated = QtWidgets.QWidget()
        self.tabRepeated.setObjectName("tabRepeated")
        self.tabWidget.addTab(self.tabRepeated, "Независимые повторные испытания")

        vboxRep = QtWidgets.QVBoxLayout(self.tabRepeated)
        vboxRep.setContentsMargins(0, 0, 0, 0)
        vboxRep.setSpacing(0)

        self.scrollRepeated = QtWidgets.QScrollArea(self.tabRepeated)
        self.scrollRepeated.setObjectName("scrollRepeated")
        self.scrollRepeated.setWidgetResizable(True)
        vboxRep.addWidget(self.scrollRepeated)

        self.tilesRepeated = QtWidgets.QWidget()
        self.tilesRepeated.setObjectName("tilesRepeated")
        self.scrollRepeated.setWidget(self.tilesRepeated)

        self.gridRepeated = QtWidgets.QGridLayout(self.tilesRepeated)
        self.gridRepeated.setContentsMargins(8, 8, 8, 8)
        self.gridRepeated.setHorizontalSpacing(12)
        self.gridRepeated.setVerticalSpacing(12)

        # Меню/статусбар (по желанию)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        MainWindow.setStatusBar(self.statusbar)