import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# 1. DOMAIN AND PHYSICAL PARAMETERS
# ============================================================

L = 1.0
U = 1.0

nx = 101
ny = 101

dx = L / (nx - 1)
dy = L / (ny - 1)

x = np.linspace(0, L, nx)
y = np.linspace(0, L, ny)

X, Y = np.meshgrid(x, y)

# Reynolds number = U*L/nu = 100
nu = 0.01

dt = 0.001

max_time_steps = 10000

# Convergence tolerances
poisson_tol = 1e-6
flow_tol = 1e-6

max_poisson_iterations = 5000


# ============================================================
# 2. INITIALIZE STREAMFUNCTION AND VORTICITY
# ============================================================

# IMPORTANT:
# array[j,i] = (y,x)

s = np.zeros((ny, nx))       # streamfunction
o = np.zeros((ny, nx))       # vorticity

u = np.zeros((ny, nx))
v = np.zeros((ny, nx))


# ============================================================
# 3. STREAMFUNCTION BOUNDARY CONDITION
# ============================================================

# All cavity walls are impermeable.
# Therefore streamfunction is constant on every wall.
# We choose that constant to be zero.

s[0, :] = 0.0        # bottom
s[-1, :] = 0.0       # top
s[:, 0] = 0.0       # left
s[:, -1] = 0.0       # right


# ============================================================
# 4. MAIN TIME-STEPPING LOOP
# ============================================================

for n in range(max_time_steps):

    # --------------------------------------------------------
    # Save old vorticity
    # --------------------------------------------------------

    o_old = o.copy()


    # ========================================================
    # 5. SOLVE POISSON EQUATION
    #
    #       ∇²ψ = -ω
    #
    # ========================================================

    for poisson_iter in range(max_poisson_iterations):

        s_old = s.copy()
        s_new = s.copy()

        # Interior points only
        for j in range(1, ny - 1):
            for i in range(1, nx - 1):

                s_new[j, i] = 0.25 * (
                    s_old[j, i + 1]
                    + s_old[j, i - 1]
                    + s_old[j + 1, i]
                    + s_old[j - 1, i]
                    + o[j, i] * dx**2
                )

        # Keep ψ = 0 on all walls
        s_new[0, :] = 0.0
        s_new[-1, :] = 0.0
        s_new[:, 0] = 0.0
        s_new[:, -1] = 0.0

        # Poisson convergence
        poisson_error = np.max(np.abs(s_new - s_old))

        s = s_new.copy()

        if poisson_error < poisson_tol:
            break


    # ========================================================
    # 6. CALCULATE VELOCITY FROM STREAMFUNCTION
    #
    #       u = ψ_y
    #       v = -ψ_x
    #
    # ========================================================

    for j in range(1, ny - 1):
        for i in range(1, nx - 1):

            u[j, i] = (
                s[j + 1, i] - s[j - 1, i]
            ) / (2 * dy)

            v[j, i] = -(
                s[j, i + 1] - s[j, i - 1]
            ) / (2 * dx)


    # --------------------------------------------------------
    # Wall velocities
    # --------------------------------------------------------

    # Bottom wall
    u[0, :] = 0.0
    v[0, :] = 0.0

    # Top moving lid
    u[-1, :] = U
    v[-1, :] = 0.0

    # Left wall
    u[:, 0] = 0.0
    v[:, 0] = 0.0

    # Right wall
    u[:, -1] = 0.0
    v[:, -1] = 0.0


    # ========================================================
    # 7. APPLY GHIA WALL VORTICITY
    # ========================================================

    # --------------------------------------------------------
    # TOP WALL
    # --------------------------------------------------------

    j = ny - 1

    for i in range(1, nx - 1):

        s_j = s[j, i]
        s_j1 = s[j - 1, i]
        s_j2 = s[j - 2, i]

        # Ghost streamfunction above the lid
        s_ghost = (
            3 * U * dy
            - 1.5 * s_j
            + 3.0 * s_j1
            - 0.5 * s_j2
        )

        # CD-2 second derivative
        o[j, i] = -(
            s_ghost
            - 2.0 * s_j
            + s_j1
        ) / dy**2


    # --------------------------------------------------------
    # BOTTOM WALL
    # --------------------------------------------------------

    j = 0

    for i in range(1, nx - 1):

        s_j = s[j, i]
        s_j1 = s[j + 1, i]
        s_j2 = s[j + 2, i]

        # Stationary wall: U_wall = 0
        s_ghost = (
            -1.5 * s_j
            + 3.0 * s_j1
            - 0.5 * s_j2
        )

        o[j, i] = -(
            s_j1
            - 2.0 * s_j
            + s_ghost
        ) / dy**2


    # --------------------------------------------------------
    # LEFT WALL
    # --------------------------------------------------------

    i = 0

    for j in range(1, ny - 1):

        s_i = s[j, i]
        s_i1 = s[j, i + 1]
        s_i2 = s[j, i + 2]

        s_ghost = (
            -1.5 * s_i
            + 3.0 * s_i1
            - 0.5 * s_i2
        )

        o[j, i] = -(
            s_i1
            - 2.0 * s_i
            + s_ghost
        ) / dx**2


    # --------------------------------------------------------
    # RIGHT WALL
    # --------------------------------------------------------

    i = nx - 1

    for j in range(1, ny - 1):

        s_i = s[j, i]
        s_i1 = s[j, i - 1]
        s_i2 = s[j, i - 2]

        s_ghost = (
            -1.5 * s_i
            + 3.0 * s_i1
            - 0.5 * s_i2
        )

        o[j, i] = -(
            s_ghost
            - 2.0 * s_i
            + s_i1
        ) / dx**2


    # ========================================================
    # 8. UPDATE INTERIOR VORTICITY
    #
    #       ω_t + uω_x + vω_y = ν∇²ω
    #
    # ========================================================

    o_new = o.copy()

    for j in range(1, ny - 1):
        for i in range(1, nx - 1):

            # ------------------------------------------------
            # UPWIND dω/dx
            # ------------------------------------------------

            if u[j, i] > 0:

                omega_x = (
                    o[j, i] - o[j, i - 1]
                ) / dx

            else:

                omega_x = (
                    o[j, i + 1] - o[j, i]
                ) / dx


            # ------------------------------------------------
            # UPWIND dω/dy
            # ------------------------------------------------

            if v[j, i] > 0:

                omega_y = (
                    o[j, i] - o[j - 1, i]
                ) / dy

            else:

                omega_y = (
                    o[j + 1, i] - o[j, i]
                ) / dy


            # ------------------------------------------------
            # CD-2 diffusion
            # ------------------------------------------------

            omega_xx = (
                o[j, i + 1]
                - 2.0 * o[j, i]
                + o[j, i - 1]
            ) / dx**2

            omega_yy = (
                o[j + 1, i]
                - 2.0 * o[j, i]
                + o[j - 1, i]
            ) / dy**2


            # ------------------------------------------------
            # Explicit Euler
            # ------------------------------------------------

            o_new[j, i] = (
                o[j, i]
                - dt * (
                    u[j, i] * omega_x
                    + v[j, i] * omega_y
                )
                + nu * dt * (
                    omega_xx
                    + omega_yy
                )
            )


    # Update vorticity
    o = o_new.copy()


    # ========================================================
    # 9. FLOW CONVERGENCE
    # ========================================================

    flow_error = np.max(np.abs(o - o_old))


    # Print progress
    if n % 100 == 0:

        print(
            "Step =", n,
            " | Poisson error =", poisson_error,
            " | Flow error =", flow_error
        )


    # Check steady state
    if flow_error < flow_tol:

        print()
        print("==========================================")
        print("FLOW CONVERGED!")
        print("Time step =", n)
        print("Flow error =", flow_error)
        print("Poisson error =", poisson_error)
        print("==========================================")

        break


# ============================================================
# 10. FINAL VELOCITY
# ============================================================

for j in range(1, ny - 1):
    for i in range(1, nx - 1):

        u[j, i] = (
            s[j + 1, i] - s[j - 1, i]
        ) / (2 * dy)

        v[j, i] = -(
            s[j, i + 1] - s[j, i - 1]
        ) / (2 * dx)


# ============================================================
# 11. STREAMFUNCTION CONTOUR
# ============================================================

plt.figure()

plt.contourf(
    X,
    Y,
    s,
    levels=30
)

plt.colorbar(label="Streamfunction ψ")

plt.xlabel("x")
plt.ylabel("y")

plt.title("Lid-Driven Cavity: Streamfunction")

plt.axis("equal")

plt.show()


# ============================================================
# 12. VORTICITY CONTOUR
# ============================================================

plt.figure()

plt.contourf(
    X,
    Y,
    o,
    levels=30
)

plt.colorbar(label="Vorticity ω")

plt.xlabel("x")
plt.ylabel("y")

plt.title("Lid-Driven Cavity: Vorticity")

plt.axis("equal")

plt.show()


# ============================================================
# 13. STREAMLINES
# ============================================================

plt.figure()

plt.streamplot(
    X,
    Y,
    u,
    v,
    density=1.5
)

plt.xlabel("x")
plt.ylabel("y")

plt.title("Lid-Driven Cavity: Streamlines")

plt.axis("equal")

plt.show()


# ============================================================
# 14. VELOCITY VECTORS
# ============================================================

skip = 5

plt.figure()

plt.quiver(
    X[::skip, ::skip],
    Y[::skip, ::skip],
    u[::skip, ::skip],
    v[::skip, ::skip]
)

plt.xlabel("x")
plt.ylabel("y")

plt.title("Lid-Driven Cavity: Velocity Field")

plt.axis("equal")

plt.show()
