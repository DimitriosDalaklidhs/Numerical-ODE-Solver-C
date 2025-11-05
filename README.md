# Numerical ODE Solver in C
This project compares three numerical integration methods (Euler, Improved Euler, and Runge–Kutta) implemented in C.
Includes automated testing, performance analysis, and reproducible build scripts.
Demonstrates algorithmic reasoning, modular C design, and CI automation.

Implements numerical methods for solving first-order **Ordinary Differential Equations (ODEs)** of the form:

> y'(t) = f(t, y),  y(a) = y₀

It provides three common numerical solvers:
- **FE (Forward Euler)**
- **CD (Central Difference / Modified Euler)**
- **IE (Implicit Euler)**

---

 Example:

 Enter initial value problem data y'(t)=f(t,y), y(a)=y0
a: 0
b: 1
y0: 1
Do you want to run for h=0.2 and h=0.1? [y/n]: y




**Example output:**

Max error FE: 2.010816e+000
Max error CD: 6.156161e-001
Max error IE: 2.882349e-001


---

##  Files
| File | Description |

| `main.c` | Main driver program (runs simulations, displays results) |
| `f_function.c` | Defines `f(t, y)` and the exact solution `y_exact(t)` |
| `f_function.h` | Function prototypes |
| `Makefile.win` | Dev-C++ build configuration |

---

## Mathematical Model
For **PROBLEM_ID = 1**,  
f(t, y) = 2y → exact solution y(t) = e^(2t)  
For **PROBLEM_ID = 2**,  
f(t, y) = 1 − 2πsin(2πt) → exact solution y(t) = t + cos(2πt)

---

## How to Compile
You can compile the project with **Dev-C++** or GCC manually:
```bash
gcc main.c f_function.c -o Project2.exe -lm


