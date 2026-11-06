# Farmacocinética Ambiental: Degradación de Medicamentos Expuestos al Sol en Ferias Populares

> **Desafío Universitario de Métodos Numéricos**  
> Simulación interactiva de degradación de medicamentos en suspensión expuestos a radiación solar ultravioleta y calor extremo en la **Feria 16 de Julio (El Alto, Bolivia, 4.000 m.s.n.m.)**.  
> Comparación de métodos de un solo paso: **Método de Heun (RK2)** vs. **Runge-Kutta de 4.° Orden (RK4)** frente a la **Solución Analítica Exacta**.

---

## 1. Descripción del Modelo y Contexto Sanitario

En ferias populares al aire libre en ciudades de gran altitud —como la Feria 16 de Julio de El Alto—, medicamentos antibióticos en suspensión (por ejemplo, amoxicilina pediátrica) son comercializados directamente sobre tarimas o bajo carpas plásticas, expuestos a **radiación solar ultravioleta extrema (Índice UV > 12)** y temperaturas superiores a 30 °C.

Bajo estas condiciones ambientales, el principio activo sufre **fotólisis y termólisis acelerada**, reduciendo su concentración eficaz antes de la fecha de caducidad. Cuando un niño o adulto vulnerable ingiere un antibiótico degradado por debajo del **90% de su concentración original ($C_0$)**, recibe una **subdosis terapéutica** que conduce a dos consecuencias graves:
1. **Fracaso del tratamiento médico:** la infección respiratoria o gastrointestinal no cede.
2. **Selección de bacterias multirresistentes (RAM):** una de las principales amenazas de salud pública global reconocidas por la OMS.

---

## 2. Ecuación Diferencial Ordinaria (EDO) y Solución Analítica

El fenómeno de degradación química de primer orden se modela mediante el siguiente problema de valor inicial (PVI):

$$\frac{dC}{dt} = -k \cdot C, \quad C(0) = C_0$$

Donde:
* $C(t)$: Concentración de principio activo en el tiempo $t$ $[\text{mg/mL}]$.
* $t$: Tiempo de exposición solar continua $[\text{horas}]$.
* $k$: Constante cinética de degradación aparente $[\text{h}^{-1}]$.
* $C_0$: Concentración inicial al momento de la reconstitución o apertura $[\text{mg/mL}]$.

### Parámetros por Defecto (Hipotéticos)
* $C_0 = 50.0 \text{ mg/mL}$
* $k = 0.025 \text{ h}^{-1}$ (Caso base representativo)
* $T = 8.0 \text{ h}$ (Jornada típica de feria: 08:00 a 16:00)
* $h = 1.0 \text{ h}$ (Paso de integración primario)
* $h_2 = 0.5 \text{ h}$ (Segundo paso para análisis de convergencia)
* Umbral de seguridad farmacopeico: $90\%$ de $C_0$ ($45.0 \text{ mg/mL}$)

> **Nota:** Todos los parámetros son **hipotéticos** con propósitos didácticos para la materia de Métodos Numéricos.

---

### Deducción de la Solución Analítica (Separación de Variables)

1. **Separación de variables:**
   $$\frac{1}{C} \, dC = -k \, dt$$

2. **Integración en ambos miembros:**
   $$\int \frac{1}{C} \, dC = \int -k \, dt \implies \ln|C| = -k \cdot t + A$$

3. **Aplicación de exponencial:**
   $$C(t) = e^{-kt + A} = e^A \cdot e^{-kt} = K \cdot e^{-kt}$$

4. **Aplicación de la condición inicial $C(0) = C_0$:**
   $$C(0) = K \cdot e^0 = K = C_0 \implies C(t) = C_0 \cdot e^{-kt}$$

5. **Comprobación por sustitución directa:**
   $$\frac{d}{dt}\left[C_0 \cdot e^{-kt}\right] = C_0 \cdot (-k) \cdot e^{-kt} = -k \cdot \left(C_0 \cdot e^{-kt}\right) = -k \cdot C(t) \quad \checkmark$$
   $$C(0) = C_0 \cdot e^0 = C_0 \quad \checkmark$$

6. **Deducción del Tiempo Crítico $t_{90}$:**
   $$C(t_{90}) = 0.90 \cdot C_0 \implies C_0 \cdot e^{-k \cdot t_{90}} = 0.90 \cdot C_0 \implies e^{-k \cdot t_{90}} = 0.90$$
   $$-k \cdot t_{90} = \ln(0.90) \implies t_{90} = \frac{\ln(1 / 0.90)}{k} = -\frac{\ln(0.90)}{k} \approx \frac{0.1053605}{k}$$

Para $k = 0.025 \text{ h}^{-1}$:
$$t_{90} = \frac{0.1053605}{0.025} \approx 4.2144 \text{ horas} \quad (\approx 4 \text{ h } 13 \text{ min})$$

---

## 3. Métodos Numéricos Implementados

Ambos métodos se implementaron como **funciones puras y reutilizables**, con ajuste automático del último paso si $T$ no es múltiplo exacto de $h$:

### 1. Método de Heun (RK2 Predictor-Corrector)
* **Predictor (Euler hacia adelante):**
  $$k_1 = f(t_n, C_n) = -k \cdot C_n$$
  $$\tilde{C} = C_n + h \cdot k_1$$
* **Corrector (Trapecio):**
  $$k_2 = f(t_n + h, \tilde{C}) = -k \cdot \tilde{C}$$
  $$C_{n+1} = C_n + \frac{h}{2} \cdot (k_1 + k_2)$$
* **Orden de convergencia:** Global $\mathcal{O}(h^2)$, local $\mathcal{O}(h^3)$.

### 2. Método de Runge-Kutta de 4.° Orden (RK4 Clásico)
$$k_1 = f(t_n, C_n)$$
$$k_2 = f\left(t_n + \frac{h}{2}, C_n + \frac{h}{2} k_1\right)$$
$$k_3 = f\left(t_n + \frac{h}{2}, C_n + \frac{h}{2} k_2\right)$$
$$k_4 = f(t_n + h, C_n + h \cdot k_3)$$
$$C_{n+1} = C_n + \frac{h}{6} \cdot (k_1 + 2k_2 + 2k_3 + k_4)$$
* **Orden de convergencia:** Global $\mathcal{O}(h^4)$, local $\mathcal{O}(h^5)$.

---

## 4. Estructura del Repositorio

El proyecto está diseñado bajo la filosofía de **cero dependencias locales**, sin procesos de *build*, empaquetado ni dependencias de `npm`:

```text
├── index.html                           # Aplicación web completa (HTML5 + CSS3 + JavaScript puro)
├── README.md                            # Documentación técnica, modelo matemático y guía de uso
├── simulador_alerta_desapariciones.html # Módulo complementario de seguridad comunitaria
└── preguntas_alucinaciones_v2.csv       # Archivo complementario del proyecto
```

---

## 5. Instrucciones de Ejecución y Publicación

### Opción A: Ejecución Local Inmediata (Sin Internet)
1. Clona o descarga este repositorio:
   ```bash
   git clone https://github.com/tu-usuario/farmacocinetica-feria.git
   ```
2. Haz doble clic en el archivo `index.html` o ábrelo directamente en tu navegador preferido (Chrome, Edge, Firefox, Safari).
3. No requiere Node.js, Python, extensiones ni configuración de servidores.

### Opción B: Publicación en GitHub Pages (Paso a Paso)
1. Crea un nuevo repositorio en [GitHub](https://github.com) (público).
2. Sube los archivos `index.html` y `README.md` a la rama principal (`main`):
   ```bash
   git init
   git add .
   git commit -m "Publicación inicial del simulador de farmacocinética"
   git branch -M main
   git remote add origin https://github.com/tu-usuario/nombre-del-repo.git
   git push -u origin main
   ```
3. En GitHub, dirígete a:  
   **Settings** (pestaña superior) $\rightarrow$ **Pages** (menú lateral izquierdo).
4. En **Build and deployment $\rightarrow$ Branch**:
   * Selecciona `main`.
   * Carpeta: `/ (root)`.
   * Presiona **Save**.
5. Espera unos 60 segundos. GitHub te proporcionará la URL pública:  
   `https://tu-usuario.github.io/nombre-del-repo/`

---

## 6. Guía Rápida para Probar la Aplicación

Para apreciar claramente la diferencia de orden de convergencia entre Heun $\mathcal{O}(h^2)$ y RK4 $\mathcal{O}(h^4)$, prueba los siguientes tamaños de paso en el simulador:

| Paso $h$ [h] | Comportamiento Heun (RK2) | Comportamiento RK4 | Observación Pedagógica |
| :---: | :---: | :---: | :--- |
| **$h = 4.0$** | Error visible ($\approx 0.013$ mg/mL) | Error ínfimo ($\approx 7 \times 10^{-7}$ mg/mL) | **Excelente contraste:** Heun se desvía visiblemente; RK4 se mantiene sobre la curva analítica. |
| **$h = 2.0$** | Error moderado ($\approx 3.4 \times 10^{-3}$ mg/mL) | Error microscópico ($\approx 4.4 \times 10^{-8}$ mg/mL) | El error de Heun se reduce a la cuarta parte ($2^2=4$). |
| **$h = 1.0$** *(Defecto)* | Error pequeño ($\approx 8.7 \times 10^{-4}$ mg/mL) | Error despreciable ($\approx 2.7 \times 10^{-8}$ mg/mL) | RK4 es $\approx 31.900$ veces más preciso que Heun. |
| **$h = 0.5$** | Error residual ($\approx 2.1 \times 10^{-4}$ mg/mL) | Error límite ($\approx 1.7 \times 10^{-9}$ mg/mL) | Permite validar la tabla de orden de convergencia ($p \approx 2$ y $p \approx 4$). |

### Pruebas de Escenarios y Elementos Interactivos
1. **Frasco animado:** Pulsa **Reproducir (▶)** para ver el líquido decolorarse en tiempo real de frambuesa/rosa a naranja y luego a gris turbio. Observa el semáforo cambiar a rojo cuando se cruza el umbral $t_{90}$.
2. **Escenarios rápidos:** Prueba el botón **"Día Muy Soleado ($k = 0.05$)"** para ver el tiempo crítico reducirse a **2 h 06 min**, o **"A la Sombra ($k = 0.01$)"** para extender la vida útil a más de **10 horas**.
3. **Escala Logarítmica:** En la gráfica de errores, pulsa **"Escala: Lineal / Logarítmica"** para apreciar con claridad los errores de RK4 en el orden de $10^{-8}$.
4. **Consola del Navegador:** Presiona `F12` y abre la **Consola**: verás la autoverificación con `console.assert` validando la condición inicial y que el error de RK4 es menor que Heun.

---

## 7. Limitaciones y Descargo de Responsabilidad

* Este software ha sido desarrollado con fines **estrictamente académicos y pedagógicos** para la materia de Métodos Numéricos.
* Las cinéticas de degradación no representan formulaciones comerciales específicas ni sustituyen estudios de estabilidad según directrices ICH Q1A(R2).
* **Consejo Sanitario:** Adquiera siempre medicamentos en farmacias autorizadas y manténgalos almacenados en lugares frescos y protegidos de la luz solar directa.

