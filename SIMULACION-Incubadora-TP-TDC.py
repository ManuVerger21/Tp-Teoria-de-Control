"""
SIMULACION INTERACTIVA - Control de Temperatura en Incubadora Neonatal
Teoria de Control - UTN FRBA - Prof. Mgtr. Omar Civale
=====================================================================
App de escritorio (PyQt5 + pyqtgraph). Mismo formato que el TP del
router, pero con modelo termico de 3 nodos (aire, estructura, neonato).

Requisitos:  pip install PyQt5 pyqtgraph numpy
Ejecutar:    python SIMULACION-Incubadora-TP-TDC.py
"""
import sys
import numpy as np
from PyQt5 import QtWidgets, QtCore
import pyqtgraph as pg

# ==========================================
# PARAMETROS FISICOS (caso nominal)
# ==========================================
DT = 0.5            # paso de integracion [s]
SPEED = 20          # pasos de simulacion por refresco (acelera el tiempo)

# Capacidades termicas [J/°C]
CAIR, CSTR, CBABY = 2289.0, 8000.0, 5250.0
# Resistencias termicas [°C/W]
RAS, RBA, RENV, RSTR = 0.020, 0.090, 0.180, 0.400
QMET = 4.0          # calor metabolico neonato [W]
PMAX = 250.0        # potencia maxima calefactor [W]
L_D  = 8.0          # retardo medicion+distribucion [s]
TAMB_NOM = 24.0     # sala nominal [°C]

# PID por Ziegler-Nichols (curva de reaccion)
KP, KI, KD = 147.6, 9.85, 553.1

HISTORY_LEN = 800


class SimuladorIncubadora(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Simulación TPI - Control de Temperatura en Incubadora Neonatal")
        self.resize(1100, 1000)

        # Estado fisico inicial: todo en equilibrio con la sala
        self.Ta = TAMB_NOM
        self.Ts = TAMB_NOM
        self.Tb = TAMB_NOM
        self.integral = 0.0
        self.prev_err = 0.0
        self.u_buffer = [0.0] * int(L_D / DT)

        # Entradas del usuario
        self.referencia = 36.5      # setpoint [°C]
        self.tamb = TAMB_NOM        # temperatura de sala [°C]
        self.qp = 0.0               # perdida por portillo [W]
        self.qevap = 0.8            # evaporacion [W]

        # Buffers de datos
        self.d_ref   = np.full(HISTORY_LEN, TAMB_NOM)
        self.d_ta    = np.full(HISTORY_LEN, TAMB_NOM)
        self.d_tb    = np.full(HISTORY_LEN, TAMB_NOM)
        self.d_u     = np.zeros(HISTORY_LEN)
        self.d_err   = np.zeros(HISTORY_LEN)

        self.setup_ui()

        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.update_simulation)
        self.timer.start(30)

    def setup_ui(self):
        pg.setConfigOption('background', 'k')
        pg.setConfigOption('foreground', 'w')

        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        layout = QtWidgets.QVBoxLayout(central)

        ctrl = QtWidgets.QHBoxLayout()

        # --- Referencia (setpoint) ---
        gb_ref = QtWidgets.QGroupBox("Referencia: Temp. deseada (Tref)")
        v = QtWidgets.QVBoxLayout()
        self.sld_ref = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.sld_ref.setRange(300, 380)   # 30.0 .. 38.0 °C (x10)
        self.sld_ref.setValue(365)
        self.lbl_ref = QtWidgets.QLabel("36.5 °C")
        v.addWidget(self.sld_ref); v.addWidget(self.lbl_ref)
        gb_ref.setLayout(v)

        # --- Temperatura de sala ---
        gb_amb = QtWidgets.QGroupBox("Perturbación: Temp. de sala")
        v = QtWidgets.QVBoxLayout()
        self.sld_amb = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.sld_amb.setRange(180, 280)   # 18 .. 28 °C
        self.sld_amb.setValue(240)
        self.lbl_amb = QtWidgets.QLabel("24.0 °C")
        v.addWidget(self.sld_amb); v.addWidget(self.lbl_amb)
        gb_amb.setLayout(v)

        # --- Portillo ---
        gb_door = QtWidgets.QGroupBox("Perturbación: Apertura portillo")
        v = QtWidgets.QVBoxLayout()
        self.sld_door = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.sld_door.setRange(0, 30)     # 0..30 W de perdida
        self.lbl_door = QtWidgets.QLabel("0 W")
        v.addWidget(self.sld_door); v.addWidget(self.lbl_door)
        gb_door.setLayout(v)

        # --- PID en vivo ---
        gb_pid = QtWidgets.QGroupBox("Sintonización PID (en vivo)")
        form = QtWidgets.QFormLayout()
        self.spin_kp = QtWidgets.QDoubleSpinBox(); self.spin_kp.setRange(0, 500); self.spin_kp.setValue(KP)
        self.spin_ki = QtWidgets.QDoubleSpinBox(); self.spin_ki.setRange(0, 200); self.spin_ki.setDecimals(2); self.spin_ki.setValue(KI)
        self.spin_kd = QtWidgets.QDoubleSpinBox(); self.spin_kd.setRange(0, 2000); self.spin_kd.setValue(KD)
        form.addRow("Kp:", self.spin_kp)
        form.addRow("Ki:", self.spin_ki)
        form.addRow("Kd:", self.spin_kd)
        gb_pid.setLayout(form)

        ctrl.addWidget(gb_ref); ctrl.addWidget(gb_amb)
        ctrl.addWidget(gb_door); ctrl.addWidget(gb_pid)
        layout.addLayout(ctrl)

        self.sld_ref.valueChanged.connect(self.update_labels)
        self.sld_amb.valueChanged.connect(self.update_labels)
        self.sld_door.valueChanged.connect(self.update_labels)

        # --- Graficos ---
        gw = pg.GraphicsLayoutWidget(); layout.addWidget(gw)

        p1 = gw.addPlot(title="Temperatura: Referencia (verde), Aire Ta (amarillo), Neonato Tb (naranja)")
        p1.setYRange(22, 39); p1.showGrid(x=True, y=True)
        p1.addLine(y=37.0, pen=pg.mkPen('g', style=QtCore.Qt.DotLine))
        p1.addLine(y=36.0, pen=pg.mkPen('g', style=QtCore.Qt.DotLine))
        self.c_ref = p1.plot(pen=pg.mkPen('g', width=2, style=QtCore.Qt.DashLine))
        self.c_ta  = p1.plot(pen=pg.mkPen('y', width=3))
        self.c_tb  = p1.plot(pen=pg.mkPen((255,140,0), width=2))

        gw.nextRow()
        p2 = gw.addPlot(title="Señal de control: Potencia del calefactor u [%]")
        p2.setYRange(-5, 105); p2.showGrid(x=True, y=True)
        p2.addLine(y=100, pen=pg.mkPen('r', style=QtCore.Qt.DotLine))
        self.c_u = p2.plot(pen=pg.mkPen('c', width=2))

        gw.nextRow()
        p3 = gw.addPlot(title="Error  e = Tref − Ta  [°C]")
        p3.setYRange(-3, 3); p3.showGrid(x=True, y=True)
        p3.addLine(y=0.5, pen=pg.mkPen('g', style=QtCore.Qt.DotLine))
        p3.addLine(y=-0.5, pen=pg.mkPen('g', style=QtCore.Qt.DotLine))
        self.c_err = p3.plot(pen=pg.mkPen('m', width=2))

    def update_labels(self):
        self.referencia = self.sld_ref.value() / 10.0
        self.tamb = self.sld_amb.value() / 10.0
        self.qp = float(self.sld_door.value())
        self.lbl_ref.setText(f"{self.referencia:.1f} °C")
        self.lbl_amb.setText(f"{self.tamb:.1f} °C")
        self.lbl_door.setText(f"{int(self.qp)} W")

    def paso_fisico(self, kp, ki, kd):
        # --- PID discreto con anti-windup condicional ---
        e = self.referencia - self.Ta
        d = (e - self.prev_err) / DT
        self.prev_err = e
        integ_tent = self.integral + e * DT
        u_raw = kp * e + ki * integ_tent + kd * d
        u_sat = float(np.clip(u_raw / 100.0, 0.0, 1.0))
        saturado = (u_raw / 100.0 > 1.0) or (u_raw / 100.0 < 0.0)
        if (not saturado) or (e * u_raw < 0):
            self.integral = integ_tent
        u = u_sat

        # --- retardo de transporte ---
        self.u_buffer.append(u)
        u_appl = self.u_buffer.pop(0)

        # --- modelo termico de 3 nodos ---
        Qh = PMAX * u_appl
        dTa = (Qh + (self.Ts - self.Ta)/RAS + (self.Tb - self.Ta)/RBA
               - (self.Ta - self.tamb)/RENV - self.qp) / CAIR
        dTs = ((self.Ta - self.Ts)/RAS - (self.Ts - self.tamb)/RSTR) / CSTR
        dTb = (QMET - (self.Tb - self.Ta)/RBA - self.qevap) / CBABY
        self.Ta += dTa * DT
        self.Ts += dTs * DT
        self.Tb += dTb * DT
        return u, e

    def update_simulation(self):
        kp = self.spin_kp.value(); ki = self.spin_ki.value(); kd = self.spin_kd.value()
        for _ in range(SPEED):          # acelera el tiempo simulado
            u, e = self.paso_fisico(kp, ki, kd)

        def roll(buf, val):
            buf[:-1] = buf[1:]; buf[-1] = val
        roll(self.d_ref, self.referencia)
        roll(self.d_ta, self.Ta)
        roll(self.d_tb, self.Tb)
        roll(self.d_u, u * 100)
        roll(self.d_err, e)

        self.c_ref.setData(self.d_ref)
        self.c_ta.setData(self.d_ta)
        self.c_tb.setData(self.d_tb)
        self.c_u.setData(self.d_u)
        self.c_err.setData(self.d_err)


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    win = SimuladorIncubadora()
    win.show()
    sys.exit(app.exec_())
