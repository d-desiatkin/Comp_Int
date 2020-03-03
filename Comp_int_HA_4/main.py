from scipy import integrate
import numpy as np
from math import cos
from math import sin
from cvxopt import solvers
from cvxopt import matrix
import matplotlib.pyplot as plt
from time import sleep
import matplotlib


# Robot parameters
l = [1, 0.7]
lcm = [x/2 for x in l]
mass = [10, 7]
I_1 = 1 / 12 * mass[0] * lcm[0] ** 2
I_2 = 1 / 12 * mass[1] * lcm[1] ** 2

# Gravity constant
g = 0

# Boundary conditions of the task
U1_max = 1000
U1_min = -1000
U2_max = 700
U2_min = -700

# Coefficients of PD controler
kp = 100
kd = 10000

# Optimization parameters
# solvers.options['show_progress'] = False
solvers.options['maxiters'] = 200


def fk(q):
    x = l[0] * cos(q[0].item()) + l[1] * cos(q[0].item() + q[1].item())
    y = l[0] * sin(q[0].item()) + l[1] * sin(q[0].item() + q[1].item())
    return x, y


def jac_1(q):
    j = np.array([[-lcm[0] * sin(q[0]), 0.],
                  [lcm[0] * cos(q[0]), 0.],
                  [0., 0.]])
    return j


def jac_2(q):
    j = np.array([[-l[0] * sin(q[0]) - lcm[1] * sin(q[0] + q[1]), -lcm[1] * sin(q[0] + q[1])],
                 [l[0] * cos(q[0]) + lcm[1] * cos(q[0] + q[1]), lcm[1] * cos(q[0] + q[1])],
                 [0., 0.]])
    return j


def optimisation(U_w):

    H = np.array([[1, 0],
                  [0, 1]], dtype=np.float64)

    F = np.zeros((2, 1), dtype=np.float64)

    A = np.array([[1, 0],
                  [0, 1],
                  [-1, 0],
                  [0, -1]], dtype=np.float64)

    B = np.array([[U1_max-U_w[0,0]],
                  [U2_max-U_w[1,0]],
                  [U_w[0,0]-U1_min],
                  [U_w[1,0]-U2_min]], dtype=np.float64)

    H = matrix(H)
    F = matrix(F)
    A = matrix(A)
    B = matrix(B)

    sol = solvers.qp(H, F, A, B, None, None, 'dsdp', None)

    U = np.array(sol['x']) + U_w

    return U


def dynamics(time, x, q_w):
    q = np.reshape(x[0:4], (4,1))
    q = np.concatenate((q, np.zeros((2,1))))
    U_w = kp * (q_w[0:2] - q[0:2]) + kd * (q_w[2:4] - q[2:4])
    U = optimisation(U_w)
    D = mass[0] * (jac_1(q[0:2, :]).T @ jac_1(q[0:2, :])) + mass[1] * (jac_2(q[0:2, :]).T @ jac_2(q[0:2, :])) + \
        np.array([[I_1 + I_2, I_2], [I_2, I_2]])

    h = - mass[1] * l[0] * lcm[1] * sin(q[1])

    C = np.array([[h*q[3].item(), h*q[3].item()+h*q[2].item()], [-h*q[2].item(), 0]])

    gr = np.array([[(mass[0] * lcm[0] + mass[1] * l[0]) * g * cos(q[0]) + mass[1] * lcm[1] * g * cos(q[0] + q[1])],
                   [mass[1] * lcm[1] * g * cos(q[0] + q[1])]])

    q[4:6] = np.linalg.pinv(D)@(U - C @ q[2:4] - gr)
    # q[4] = q[4].item()
    # q[5] = q[5].item()

    x_der = np.zeros((6,1))

    x_der[0:2] = q[2:4]
    x_der[2:4] = q[4:6]
    x_der[4:6] = U
    x_der = x_der.flatten()
    return x_der



time = np.linspace(0, 2000, 20001)

#  Initial position
q0 = np.array([[np.pi/4], [np.pi/4], [0], [0], [0], [0]])
x0 = np.concatenate((q0[0:4], np.zeros((2, 1)))).flatten()

# Wanted position
q_w = np.array([[np.pi/3], [np.pi/3], [0], [0], [0], [0]])

# dynamics(x0, time, q_w)
ans = integrate.odeint(dynamics, x0, time, (q_w,), tfirst=True)
# hh_ans = integrate.ode(dynamics)
# hh_ans.set_initial_value(x0, 0.0)
# hh_ans.set_f_params(q_w)
# hh_ans.set_integrator("dopri5", rtol=0.0001, max_step=5000)
# ans = []
# for i in time:
#     ans.append(hh_ans.integrate(i))

plt.figure()
x_s, y_s = fk(x0.flatten())
x_t, y_t = fk(q_w[0:2, 0])
x = []
y = []
for i in range(len(ans)):
    x_, y_ = fk(ans[i, 0:2])
    x.append(x_)
    y.append(y_)
plt.plot(x, y, c='b')
plt.scatter(x_t, y_t, s=20, c='r')
plt.scatter(x_s, y_s, s=20, c='black')
plt.show()
plt.figure()
plt.plot(ans[:, 4])
plt.plot(ans[:, 5])
plt.show()

# print(np.array([1, 2]) @ np.array([[1],[2]]))
# print(ans)
# print(q_w.flatten())