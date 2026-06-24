# Simulador de Control de Temperatura (PID) – Incubadora Neonatal

Trabajo Práctico Integrador Final – **Teoría de Control** – UTN FRBA
Docente: Prof. Mgtr. Omar Civale
Tema: Control de temperatura en incubadoras neonatales (modo aire)

## Descripción del sistema

La simulación implementa un sistema de control en **lazo cerrado** que regula la
temperatura del aire interno de una incubadora neonatal mediante un controlador
**PID digital**. A diferencia de un modelo de un solo bloque, la planta se modela
con **tres nodos térmicos** acoplados:

- **Aire interno (Ta)** – variable controlada.
- **Estructura + colchón + cúpula (Ts)** – masa térmica de acople.
- **Neonato (Tb)** – carga principal (1,5 kg).

Esto permite mostrar que la respuesta depende de la carga térmica y no solo del
controlador: el aire cambia rápido, pero el bebé y la estructura actúan como
masas que almacenan y liberan calor.

## Identificación de la planta y sintonía

A partir de la **curva de reacción** (escalón de potencia en lazo abierto) se
extraen los parámetros de un modelo de primer orden con retardo:

| Parámetro | Valor | Significado |
|-----------|-------|-------------|
| K   | ≈ 0,28 °C/% | ganancia estática |
| τ   | ≈ 258 s     | constante de tiempo |
| L   | ≈ 7,5 s     | retardo (medición + distribución) |
| L/τ | ≈ 0,03      | < 0,3 → fácilmente controlable |

Con el método de **Ziegler-Nichols (curva de reacción)** se obtienen:

- **Kp = 147,6**
- **Ki = 9,85**  (Ti = 15 s)
- **Kd = 553,1** (Td = 3,75 s)

El controlador incorpora **anti-windup condicional**: la acción integral se
congela cuando el calefactor satura al 100 %, evitando el sobreimpulso excesivo
típico de Ziegler-Nichols (sin él, el pico supera los 39 °C; con él, queda en
36,53 °C, dentro de la banda ±0,5 °C).

## Desempeño (cambio de setpoint 35 → 36,5 °C)

- Sobreimpulso: **0,03 °C**
- Error en estado estable: **≈ 0,001 °C**
- Tiempo de establecimiento: **≈ 537 s (8,9 min)**
- Rechazo de perturbaciones (portillo, sala a 22 °C, evaporación): desviación
  máxima < 0,05 °C.

## Archivos

- `SIMULACION-Incubadora-TP-TDC.py` – **App interactiva** (PyQt5 + pyqtgraph),
  con sliders en vivo (referencia, sala, portillo) y sintonización PID editable.
- `simulacion_incubadora.py` – **Script de análisis** (matplotlib): genera la
  curva de reacción y la respuesta de lazo cerrado con todos los escenarios.

## Requisitos

```
pip install PyQt5 pyqtgraph numpy matplotlib scipy
```

## Ejecución

App interactiva:
```
python SIMULACION-Incubadora-TP-TDC.py
```

Análisis con gráficos:
```
python simulacion_incubadora.py
```

## Uso de la app interactiva

1. Mover el slider **Referencia** para fijar la temperatura deseada y ver cómo
   el aire (amarillo) la alcanza, mientras el neonato (naranja) sube más lento.
2. Bajar el slider **Temp. de sala** para simular una sala fría: el PID compensa
   subiendo la potencia.
3. Subir el slider **Apertura portillo** para inyectar una pérdida brusca y
   observar la recuperación.
4. Editar **Kp, Ki, Kd** en vivo para ver el efecto de cada acción.

## Autores

- Pagés, Rubén
- Verger, Manuel Augusto
