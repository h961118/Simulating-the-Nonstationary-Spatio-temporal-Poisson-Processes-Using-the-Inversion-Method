# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 16:49:49 2020

@author: Administrator
"""
import math
import random
import CubicEquationSolver

# input the number of realizations
N = 1

# input the parameters of the model
# input the partition of x-dimension as [x_0, x_1, ..., x_K]
x = [1, 2, 3]
# input the partition of y-dimension as [y_0, y_1, ..., y_J]
y = [1, 2, 3]
# input the partition of time interval as [s_0, s_1, ..., s_p]
t = [0, 1, 2]

# input the intensity at each knot as
# q = [[[z_0,0,0, z_0,1,0, ..., z_0,J,0], [z_1,0,0, z_1,1,0, ..., z_1,J,0], ..., [z_K,0,0, z_K,1,0, ..., z_K,J,0]],
#      [[z_0,0,1, z_0,1,1, ..., z_0,J,1], [z_1,0,1, z_1,1,1, ..., z_1,J,1], ..., [z_K,0,1, z_K,1,1, ..., z_K,J,1]],
#      ...
#      [[z_0,0,p, z_0,1,p, ..., z_0,J,p], [z_1,0,p, z_1,1,p, ..., z_1,J,p], ..., [z_K,0,p, z_K,1,p, ..., z_K,J,p]]]
q = [[[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]],
     [[1, 1, 1], [1, 1, 1], [1, 1, 1]]]


# define the functions, preparing for the simulation of our method
# define the integrated intensity m(t), and calculate the coefficients of m(t)
def cal_para(x, y, t, q):
    N_p = len(t)
    N_j = len(y)
    N_k = len(x)

    zeros = []
    for z in range(0, N_j):
        zeros.append(0)
    a = [zeros]
    b = [zeros]
    c = [zeros]
    d = [zeros]
    s = [zeros]
    for k in range(1, N_k):
        a_k = [0]
        b_k = [0]
        c_k = [0]
        d_k = [0]
        s_k = [0]
        xd = x[k] - x[k - 1]
        for j in range(1, N_j):
            yd = y[j] - y[j - 1]
            a_kj = (1 / 6) * xd * yd * (2 * x[k - 1] + x[k])
            b_kj = (1 / 6) * xd * yd * (2 * y[j - 1] + y[j])
            c_kj = (1 / 6) * xd * yd * (x[k - 1] + 2 * x[k])
            d_kj = (1 / 6) * xd * yd * (y[j - 1] + 2 * y[j])
            s_kj = xd * yd
            a_k.append(a_kj)
            b_k.append(b_kj)
            c_k.append(c_kj)
            d_k.append(d_kj)
            s_k.append(s_kj)
        a.append(a_k)
        b.append(b_k)
        c.append(c_k)
        d.append(d_k)
        s.append(s_k)

    f = []
    for p in range(0, N_p):
        f_p = []
        for k in range(1, N_k):
            f_k = []
            xd = x[k] - x[k - 1]
            for j in range(1, N_j):
                yd = y[j] - y[j - 1]
                # sum the integrals up
                f_pkj = (-a[k][j] * yd - b[k][j] * xd + s[k][j] * (x[k] * y[j] - x[k - 1] * y[j - 1]) / 2) * \
                    q[p][k - 1][j - 1] + \
                    (a[k][j] * yd + s[k][j] * (x[k - 1] * y[j - 1] - x[k - 1] * y[j]) / 2) * q[p][k][j - 1] + \
                    (b[k][j] * xd + s[k][j] * (x[k - 1] * y[j - 1] - x[k] * y[j - 1]) / 2) * q[p][k - 1][j] + \
                    (-c[k][j] * yd + s[k][j] * (x[k] * y[j] - x[k] * y[j - 1]) / 2) * q[p][k - 1][j] + \
                    (-d[k][j] * xd + s[k][j] * (x[k] * y[j] - x[k - 1] * y[j]) / 2) * q[p][k][j - 1] + \
                    (c[k][j] * yd + d[k][j] * xd + s[k][j] *
                     (x[k - 1] * y[j - 1] - x[k] * y[j]) / 2) * q[p][k][j]
                f_k.append(f_pkj)
            f_p.append(f_k)
        f.append(f_p)

    m = []
    for p in range(0, N_p):
        m_p = 0
        for k in range(1, N_k):
            for j in range(1, N_j):
                m_p = m_p + f[p][k - 1][j - 1] / s[k][j]
        m.append(m_p)

    para = []
    for p in range(0, N_p - 1):
        k = (m[p + 1] - m[p]) / (t[p + 1] - t[p])
        d = m[p] - k * t[p]
        para.append([k, d])

    return para


# attain the parameters of m(t)
linpara = cal_para(x, y, t, q)


# attain the value of m(t) at any time t
def cal_m(ti):
    for p in range(1, len(t)):
        if ti <= t[p]:
            value = linpara[p - 1][0] * ti + linpara[p - 1][1]
            break
        p = p + 1
    return value


# calculate the intensity function at vertices at any time t, \lambda(x_k,y_j,t)
def cal_g(ti, k, j):
    N_p = len(t)
    p = 1
    for p in range(1, N_p):
        if ti <= t[p]:
            break
        p = p + 1
    p = p - 1
    fraction = (ti - t[p]) / (t[p + 1] - t[p])
    g_value = q[p][k][j] + fraction * (q[p + 1][k][j] - q[p][k][j])
    return g_value


# calculate the p.d.f. of y conditional on time t
def cal_h(ti, yi):
    N_p = len(t)
    N_j = len(y)
    N_k = len(x)
    m_ti = cal_m(ti)

    p = 1
    for p in range(1, N_p):
        if ti <= t[p]:
            break
        p = p + 1
    p = p - 1

    j = 0
    while j < N_j:
        if yi <= y[j]:
            break
        j = j + 1

    h = 0
    for k in range(1, N_k):
        xd = x[k] - x[k - 1]
        yd = y[j] - y[j - 1]
        y2d = (y[j] - y[j - 1]) ** 2
        yid = yi - y[j - 1]
        akj = (xd / (2 * y2d)) * (yi - y[j]) * (yi * xd +
                                                2 * y[j - 1] * x[k - 1] - y[j] * (x[k] + x[k - 1]))
        bkj = (-xd / yd) * (yi - y[j]) * yi
        ckj = (-xd / (2 * y2d)) * yid * (yi * xd - 2 *
                                         y[j] * x[k] + y[j - 1] * (x[k] + x[k - 1]))
        dkj = (xd / yd) * yid * yi
        skj = xd * yd
        hkj = (-akj * yd - bkj * xd + ((-xd / yd) * (yi - y[j])) * (x[k] * y[j] - x[k - 1] * y[j - 1])) * cal_g(
            ti,
            k - 1,
            j - 1) + \
            (akj * yd + ((-xd / yd) * (yi - y[j])) * (x[k - 1] * y[j - 1] - x[k - 1] * y[j])) * cal_g(ti, k,
                                                                                                      j - 1) + \
            (bkj * xd + ((-xd / yd) * (yi - y[j])) * (x[k - 1] * y[j - 1] - x[k] * y[j - 1])) * cal_g(ti, k - 1,
                                                                                                      j) + \
            (-ckj * yd + ((xd / yd) * yid) * (x[k] * y[j] - x[k] * y[j - 1])) * cal_g(ti, k - 1, j) + \
            (-dkj * xd + ((xd / yd) * yid) * (x[k] * y[j] - x[k - 1] * y[j])) * cal_g(ti, k, j - 1) + \
            (ckj * yd + dkj * xd + ((xd / yd) * yid) *
             (x[k - 1] * y[j - 1] - x[k] * y[j])) * cal_g(ti, k, j)
        h = h + hkj / skj
    h = h / m_ti
    return h


# calculate the c.d.f. of y conditional on t
def cal_ycdf(ti):
    N_j = len(y)
    N_k = len(x)
    R_2 = [0]
    R_1 = [0]
    R_0 = [0]
    m_ti = cal_m(ti)
    for j in range(1, N_j):
        R_2j = 0
        R_1j = 0
        R_0j = 0
        for k in range(1, N_k):
            xd = x[k] - x[k - 1]
            yd = y[j] - y[j - 1]
            y2d = (y[j] - y[j - 1]) ** 2
            x2d = (x[k] - x[k - 1]) ** 2
            akj2 = x2d / (2 * y2d)
            akj1 = xd * (x[k - 1] * y[j - 1] - x[k] * y[j]) / y2d
            akj0 = (y[j] ** 2 * (x[k] + x[k - 1]) - 2 * x[k - 1] *
                    y[j - 1] * y[j]) * (x[k] - x[k - 1]) / (2 * y2d)
            bkj2 = -xd / yd
            bkj1 = xd * y[j] / yd
            bkj0 = 0
            ckj2 = -x2d / (2 * y2d)
            ckj1 = -xd * (x[k - 1] * y[j - 1] - x[k] * y[j]) / y2d
            ckj0 = (y[j - 1] ** 2 * (x[k] + x[k - 1]) - 2 * x[k] *
                    y[j - 1] * y[j]) * (x[k] - x[k - 1]) / (2 * y2d)
            dkj2 = xd / yd
            dkj1 = -xd * y[j - 1] / yd
            dkj0 = 0
            skj = xd * yd

            qua_kj = ((-akj2 * yd - bkj2 * xd) * cal_g(ti, k - 1, j - 1) +
                      (akj2 * yd) * cal_g(ti, k, j - 1) +
                      (bkj2 * xd) * cal_g(ti, k - 1, j) +
                      (-ckj2 * yd) * cal_g(ti, k - 1, j) +
                      (-dkj2 * xd) * cal_g(ti, k, j - 1) +
                      (ckj2 * yd + dkj2 * xd) * cal_g(ti, k, j)) / (m_ti)

            lin_kj = ((-akj1 * yd - bkj1 * xd + bkj2 * (x[k] * y[j] - x[k - 1] * y[j - 1])) * cal_g(ti, k - 1,
                                                                                                    j - 1) +
                      (akj1 * yd + bkj2 * (x[k - 1] * y[j - 1] - x[k - 1] * y[j])) * cal_g(ti, k, j - 1) +
                      (bkj1 * xd + bkj2 * (x[k - 1] * y[j - 1] - x[k] * y[j - 1])) * cal_g(ti, k - 1, j) +
                      (-ckj1 * yd + dkj2 * (x[k] * y[j] - x[k] * y[j - 1])) * cal_g(ti, k - 1, j) +
                      (-dkj1 * xd + dkj2 * (x[k] * y[j] - x[k - 1] * y[j])) * cal_g(ti, k, j - 1) +
                      (ckj1 * yd + dkj1 * xd + dkj2 * (x[k - 1] * y[j - 1] - x[k] * y[j])) * cal_g(ti, k, j)) / (
                m_ti)

            con_kj = ((-akj0 * yd - bkj0 * xd + bkj1 * (x[k] * y[j] - x[k - 1] * y[j - 1])) * cal_g(ti, k - 1,
                                                                                                    j - 1) +
                      (akj0 * yd + bkj1 * (x[k - 1] * y[j - 1] - x[k - 1] * y[j])) * cal_g(ti, k, j - 1) +
                      (bkj0 * xd + bkj1 * (x[k - 1] * y[j - 1] - x[k] * y[j - 1])) * cal_g(ti, k - 1, j) +
                      (-ckj0 * yd + dkj1 * (x[k] * y[j] - x[k] * y[j - 1])) * cal_g(ti, k - 1, j) +
                      (-dkj0 * xd + dkj1 * (x[k] * y[j] - x[k - 1] * y[j])) * cal_g(ti, k, j - 1) +
                      (ckj0 * yd + dkj0 * xd + dkj1 * (x[k - 1] * y[j - 1] - x[k] * y[j])) * cal_g(ti, k, j)) / (
                m_ti)

            R_2j = R_2j + qua_kj / skj
            R_1j = R_1j + lin_kj / skj
            R_0j = R_0j + con_kj / skj

        R_2.append(R_2j)
        R_1.append(R_1j)
        R_0.append(R_0j)

    C = [0]
    C_1 = -(R_2[1] * y[0] ** 3 / 3 + R_1[1] * y[0] ** 2 / 2 + R_0[1] * y[0])
    C.append(C_1)
    for j in range(1, N_j - 1):
        C_j1 = C[j] + y[j] ** 3 * (R_2[j] - R_2[j + 1]) / 3 + y[j] ** 2 * (R_1[j] - R_1[j + 1]) / 2 + y[j] * (
            R_0[j] - R_0[j + 1])
        C.append(C_j1)

    return R_0, R_1, R_2, C


# simulate the y-coordinate conditional on t with inversion method
def inverse_y(ti, u):
    N_j = len(y)
    R_0, R_1, R_2, C = cal_ycdf(ti)
    j = 1
    yloc = 0.5
    while j < N_j:
        if u < R_2[j] * y[j] ** 3 / 3 + R_1[j] * y[j] ** 2 / 2 + R_0[j] * y[j] + C[j]:
            break
        j = j + 1
    if abs(R_2[j]) > 0.00001:
        y_list = CubicEquationSolver.solve(
            R_2[j] / 3, R_1[j] / 2, R_0[j], C[j] - u)
    else:
        if abs(R_1[j]) > 0.00001:
            y_list = [(-R_0[j] + math.sqrt(R_0[j] ** 2 - 2 * R_1[j] * (C[j] - u))) / R_1[j],
                      (-R_0[j] - math.sqrt(R_0[j] ** 2 - 2 * R_1[j] * (C[j] - u))) / R_1[j], 0]
        else:
            y_list = [(u - C[j]) / R_0[j], 0, 0]

    k = 0
    while k <= 2:
        if abs((y_list[k] - complex(y_list[k]).real)) * 100000 - 1 < 0 and y[j - 1] <= complex(y_list[k]).real <= y[j]:
            yloc = complex(y_list[k]).real
            break
        k = k + 1
    return yloc


# calculate the pdf of x conditional on y and t
def xpdf(ti, yi):
    N_k = len(x)
    N_j = len(y)
    m_ti = cal_m(ti)
    hy = cal_h(ti, yi)

    j = 1
    while j < N_j:
        if yi <= y[j]:
            break
        j = j + 1

    z = []
    for element in x:
        z.append(element)

    q_cap = [[0, 0, 0]]
    for k in range(1, N_k):
        xkp = ((yi - y[j - 1]) / (y[j - 1] - y[j])) * (x[k] - x[k - 1]) + x[k]
        z.insert(2 * k - 1, xkp)
        xd = x[k] - x[k - 1]
        yd = y[j] - y[j - 1]
        skj = xd * yd
        xi = x[k - 1]
        hx_start = (1 / skj) * ((1 / (m_ti * hy))) * (
            ((-yd * xi - yi * xd + (x[k] * y[j] - x[k - 1] * y[j - 1])) * cal_g(ti, k - 1, j - 1)) +
            ((yd * xi + x[k - 1] * y[j - 1] - x[k - 1] * y[j]) * cal_g(ti, k, j - 1)) +
            ((yi * xd - y[j - 1] * xd) * cal_g(ti, k - 1, j)))
        xi = xkp
        hx_med = (1 / skj) * ((1 / (m_ti * hy))) * (
            ((-yd * xi - yi * xd + (x[k] * y[j] - x[k - 1] * y[j - 1])) * cal_g(ti, k - 1, j - 1)) +
            ((yd * xi + x[k - 1] * y[j - 1] - x[k - 1] * y[j]) * cal_g(ti, k, j - 1)) +
            ((yi * xd - y[j - 1] * xd) * cal_g(ti, k - 1, j)))
        xi = x[k]
        hx_end = (1 / skj) * ((1 / (m_ti * hy))) * (
            ((-yd * xi + (x[k] * y[j] - x[k] * y[j - 1])) * cal_g(ti, k - 1, j)) +
            ((-xd * yi + x[k] * y[j] - x[k - 1] * y[j]) * cal_g(ti, k, j - 1)) +
            ((yd * xi + yi * xd + x[k - 1] * y[j - 1] - x[k] * y[j]) * cal_g(ti, k, j)))
        k1 = (hx_med - hx_start) / (xkp - x[k - 1])
        q1 = [0, k1, hx_start - k1 * x[k - 1]]
        k2 = (hx_end - hx_med) / (x[k] - xkp)
        q2 = [0, k2, hx_end - k2 * x[k]]
        q_cap.append(q1)
        q_cap.append(q2)
    return q_cap, z


# simulate the x-coordinate with inversion method conditional on y and t
def inverse_x(ti, yi, u):
    Q, z = xpdf(ti, yi)

    D = [0]
    D_1 = - Q[1][1] * z[0] ** 2 / 2 - Q[1][2] * z[0]
    D.append(D_1)
    for k in range(1, len(z) - 1):
        D_k1 = D[k] + z[k] ** 2 * \
            (Q[k][1] - Q[k + 1][1]) / 2 + z[k] * (Q[k][2] - Q[k + 1][2])
        D.append(D_k1)

    z_cdf_value = []
    for k in range(0, len(z)):
        z_cdf_value.append(Q[k][1] * z[k] ** 2 / 2 + Q[k][2] * z[k] + D[k])

    k = 0
    while 1:
        if u <= z_cdf_value[k]:
            break
        k = k + 1

    # solve the quadratic equation
    if abs(Q[k][1]) >= 0.00001:
        if z[k - 1] <= (-Q[k][2] + math.sqrt(Q[k][2] ** 2 - 2 * Q[k][1] * (D[k] - u))) / Q[k][1] <= z[k]:
            xi = (-Q[k][2] + math.sqrt(Q[k][2] ** 2 -
                                       2 * Q[k][1] * (D[k] - u))) / Q[k][1]
        else:
            xi = (-Q[k][2] - math.sqrt(Q[k][2] ** 2 -
                                       2 * Q[k][1] * (D[k] - u))) / Q[k][1]
    else:
        xi = (u - D[k]) / Q[k][2]
    return xi


# simulate N realizations
# define the null set to store all the arrivals in N realizations
all_arrivals = []

for realization in range(0, N):
    # define the null set to store the arrivals in i-th realization
    arr_i = []
    t_i = 0
    k = 1
    while k <= len(t) - 1:

        u = random.random()
        uk = 1 - math.exp(-(linpara[k - 1][0] / 2) * (t[k] **
                                                      2 - t_i ** 2) - linpara[k - 1][1] * (t[k] - t_i))

        if u > uk:
            u = (u - uk) / (1 - uk)
            t_i = t[k]
            k = k + 1
        else:
            if abs(linpara[k - 1][0]) > 0.00001:
                t_i = (-linpara[k - 1][1] + math.sqrt(
                    linpara[k - 1][1] ** 2 + (linpara[k - 1][0] ** 2) * (t_i ** 2) + 2 * linpara[k - 1][0] *
                    linpara[k - 1][1] * t_i - 2 *
                    linpara[k - 1][0] * math.log(1 - u))) / linpara[k - 1][0]
            else:
                t_i = t_i - math.log(1 - u) / linpara[k - 1][1]

            # generate y conditional on t
            u = random.random()
            yloc = inverse_y(t_i, u)
            # generate x conditional on t and y
            v = random.random()
            xloc = inverse_x(t_i, yloc, v)
            arr_i.append((t_i, yloc, xloc))

    all_arrivals.append(arr_i)

# the simulated arrivals in N realizations is denoted as all_arrivals, which is a list of N lists
# all_arrivals = [[(t_0,y_0,x_0),(t_1,y_1,x_1),...],...]
# the representation of one arrival above is based on the order of simulation
# in each realization(list), the sequence is ascending with the epoch of arrival t_i
print(all_arrivals)
