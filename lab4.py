import random
import numpy as np
from scipy.stats import t,f


def kohren(dispersion, m):
    fisher = fisher_t(0.95, 1, (m - 1) * 4)
    gt = fisher/(fisher+(m-1)-2)
    gp = max(dispersion) / sum(dispersion)
    return gp < gt


def student(dispersion_reproduction, m, y_mean, xn):
    tt = 0
    f3 = (m - 1) * 4
    prob = 0.95

    x_vec = [i * 0.0001 for i in range(int(5 / 0.0001))]
    par = 0.5 + prob / 0.1 * 0.05
    for i in x_vec:
        if abs(t.cdf(i, f3) - par) < 0.000005:
            tt = i
            break

    dispersion_statistic_mark = (dispersion_reproduction / (4 * m)) ** 0.5

    beta = [1 / 4 * sum(y_mean[j] for j in range(4))]
    for i in range(3):
        b = 0
        for j in range(4):
            b += y_mean[j] * xn[j][i]
        beta.append(1 / 4 * b)

    te = []
    for i in beta:
        te.append(abs(i) / dispersion_statistic_mark)

    return te[0] > tt, te[1] > tt, te[2] > tt, te[3] > tt


def fisher_t(prob, d, f3):
    x_vec = [i*0.001 for i in range(int(10/0.001))]
    for i in x_vec:
        if abs(f.cdf(i, 4-d, f3)-prob) < 0.0001:
            return i


def fisher(m, d, y_mean, yo, dispersion_reproduction):
    dispersion_ad = 0
    for i in range(4):
        dispersion_ad += (yo[i] - y_mean[i]) ** 2

    dispersion_ad = dispersion_ad * m / (4 - d)

    fp = dispersion_ad / dispersion_reproduction
    f3 = (m - 1) * 4

    return fp < fisher_t(0.95, d, f3)


while True:
    x1_min = 10
    x1_max = 50
    x2_min = 25
    x2_max = 65
    x3_min = 50
    x3_max = 65
    xm_min = (x1_min + x2_min + x3_min) / 3
    xm_max = (x1_max + x2_max + x3_max) / 3
    y_min = 200 + xm_min
    y_max = 200 + xm_max
    x_norm = [[-1, -1, -1],
              [-1, 1, 1],
              [1, -1, 1],
              [1, 1, -1]]

    x_arr = [[20, 5, 20],
             [20, 40, 45],
             [70, 5, 45],
             [70, 40, 20]]

    m = 2
    y_arr = [[random.randint(int(y_min), int(y_max)) for i in range(m)] for j in range(4)]

    print("X1   X2  X3  Y1   Y2")
    for i in range(4):
        print(x_arr[i], y_arr[i])

    while True:
        y_aver = []
        for i in range(4):
            y_aver.append(sum(y_arr[i]) / m)

        dispersion = []
        for i in range(len(y_arr)):
            dispersion.append(0)
            for j in range(m):
                dispersion[i] += (y_aver[i] - y_arr[i][j]) ** 2
            dispersion[i] /= m

        dispersion_reproduction = sum(dispersion) / 4

        if kohren(dispersion, m):
            break
        else:
            m += 1
            for i in range(4):
                y_arr[i].append(random.randint(int(y_min), int(y_max)))

    k = student(dispersion_reproduction, m, y_aver, x_norm)
    d = sum(k)

    mx1 = (x_arr[0][0] + x_arr[1][0] + x_arr[2][0] + x_arr[3][0]) / 4
    mx2 = (x_arr[0][1] + x_arr[1][1] + x_arr[2][1] + x_arr[3][1]) / 4
    mx3 = (x_arr[0][2] + x_arr[1][2] + x_arr[2][2] + x_arr[3][2]) / 4
    my = sum(y_aver) / 4

    a11 = (x_arr[0][0] ** 2 + x_arr[1][0] ** 2 + x_arr[2][0] ** 2 + x_arr[3][0] ** 2) / 4
    a22 = (x_arr[0][1] ** 2 + x_arr[1][1] ** 2 + x_arr[2][1] ** 2 + x_arr[3][1] ** 2) / 4
    a33 = (x_arr[0][2] ** 2 + x_arr[1][2] ** 2 + x_arr[2][2] ** 2 + x_arr[3][2] ** 2) / 4
    a12 = (x_arr[0][0] * x_arr[0][1] + x_arr[1][0] * x_arr[1][1] + x_arr[2][0] * x_arr[2][1] + x_arr[3][0] * x_arr[3][1]) / 4
    a13 = (x_arr[0][0] * x_arr[0][2] + x_arr[1][0] * x_arr[1][2] + x_arr[2][0] * x_arr[2][2] + x_arr[3][0] * x_arr[3][2]) / 4
    a23 = (x_arr[0][1] * x_arr[0][2] + x_arr[1][1] * x_arr[1][2] + x_arr[2][1] * x_arr[2][2] + x_arr[3][1] * x_arr[3][2]) / 4

    a1 = (x_arr[0][0] * y_aver[0] + x_arr[1][0] * y_aver[1] + x_arr[2][0] * y_aver[2] + x_arr[3][0] * y_aver[3]) / 4
    a2 = (x_arr[0][1] * y_aver[0] + x_arr[1][1] * y_aver[1] + x_arr[2][1] * y_aver[2] + x_arr[3][1] * y_aver[3]) / 4
    a3 = (x_arr[0][2] * y_aver[0] + x_arr[1][2] * y_aver[1] + x_arr[2][2] * y_aver[2] + x_arr[3][2] * y_aver[3]) / 4

    b_arr = []
    b_arr.append(np.linalg.det(np.array([[my, mx1, mx2, mx3],
                                         [a1, a11, a12, a13],
                                         [a2, a12, a22, a23],
                                         [a3, a13, a23, a33]]))/np.linalg.det(np.array([[1, mx1, mx2, mx3],
                                                                                        [mx1, a11, a12, a13],
                                                                                        [mx2, a12, a22, a23],
                                                                                        [mx3, a13, a23, a33]])))
    b_arr.append(np.linalg.det(np.array([[1, my, mx2, mx3],
                                         [mx1, a1, a12, a13],
                                         [mx2, a2, a22, a23],
                                         [mx3, a3, a23, a33]])) / np.linalg.det(np.array([[1, mx1, mx2, mx3],
                                                                                          [mx1, a11, a12, a13],
                                                                                          [mx2, a12, a22, a23],
                                                                                          [mx3, a13, a23, a33]])))
    b_arr.append(np.linalg.det(np.array([[1, mx1, my, mx3],
                                         [mx1, a11, a1, a13],
                                         [mx2, a12, a2, a23],
                                         [mx3, a13, a3, a33]])) / np.linalg.det(np.array([[1, mx1, mx2, mx3],
                                                                                          [mx1, a11, a12, a13],
                                                                                          [mx2, a12, a22, a23],
                                                                                          [mx3, a13, a23, a33]])))
    b_arr.append(np.linalg.det(np.array([[1, mx1, mx2, my],
                                         [mx1, a11, a12, a1],
                                         [mx2, a12, a22, a2],
                                         [mx3, a13, a23, a3]])) / np.linalg.det(np.array([[1, mx1, mx2, mx3],
                                                                                          [mx1, a11, a12, a13],
                                                                                          [mx2, a12, a22, a23],
                                                                                          [mx3, a13, a23, a33]])))
    """b_arr = ([b_arr[i] * k[i] for i in range(4)])"""
    """Попередній рядок можна роздокументовати, це відкине незначущі коефіцієнти, але збільшить похибку"""

    yo = []
    for i in range(4):
        yo.append(b_arr[0] + b_arr[1] * x_arr[i][0] + b_arr[2] * x_arr[i][1] + b_arr[3] * x_arr[i][2])

    if d == 4:
        m += 1
        for i in range(4):
            y_arr[i].append(random.randint(int(y_min), int(y_max)))

    elif fisher(m, d, y_aver, yo, dispersion_reproduction):    #Дуже часто адекватно, тому для перевірки коду 4 лаби краще поставити: not fisher(...)   
        print("Рівняння лінійної регресії без еффекту взаємодії адекватно оригіналу при рівні значимості 0,05")

        '''checking'''
        checks = []
        errors = 0
        for i in range(4):
            print("Y average", i, " = ", y_aver[i])
            print("Y found  ", i, " = ", yo[i])
            checks.append(round(y_aver[i], 10) == round(yo[i], 10))
        for j in range(4):
            if not checks[j]:
                errors += 1
                print("Our test failed!")

        if errors == 0:
            print("Successful check!")
        break
        """Тут программа переходить до рівняння взаємодії, коли не адекватно"""
    else: 
        print("Рівняння лінійної регресії без еффекту взаємодії неадекватно оригіналу при рівні значимості 0,05")
        print("Переходимо до рівняння регресії з ефектом взаємодії") 
        X1_MIN, X1_MAX = -10, 50
        X2_MIN, X2_MAX = 25, 65
        X3_MIN, X3_MAX = -10, 15
        m = 3
        N = 8

        mx_max = (X1_MAX + X2_MAX + X3_MAX) / 3
        mx_min = (X1_MIN + X2_MIN + X3_MIN) / 3
        Y_MAX = mx_max + 200
        Y_MIN = mx_min + 200

        x_matrix = [
            [X1_MIN, X2_MIN, X3_MIN],
            [X1_MIN, X2_MIN, X3_MAX],
            [X1_MIN, X2_MAX, X3_MIN],
            [X1_MIN, X2_MAX, X3_MAX],
            [X1_MAX, X2_MIN, X3_MIN],
            [X1_MAX, X2_MIN, X3_MAX],
            [X1_MAX, X2_MAX, X3_MIN],
            [X1_MAX, X2_MAX, X3_MAX]
        ]

        while True:  # цикл на проверку однородности дисперсии
            y_list = [[random.randint(int(Y_MIN), int(Y_MAX)) for i in range(m)] for j in range(N)]
            my_list = []
            mx1 = 0
            mx2 = 0
            mx3 = 0

            for obj in y_list:  # создание списка my
                my_list.append(sum(obj) / len(obj))

            for obj in x_matrix:
                mx1 += obj[0]
                mx2 += obj[1]
                mx3 += obj[2]

            mx1 /= 8
            mx2 /= 8
            mx3 /= 8
            my = sum(my_list) / 8

            x123i = 0; x1i = 0; x2i = 0; x3i = 0; x12i = 0; x13i = 0; x23i = 0
            m00 = 0; m01 = 0; m02 = 0; m03 = 0; m04 = 0; m05 = 0; m06 = 0; m07 = 0
            m10 = 0; m11 = 0; m12 = 0; m13 = 0; m14 = 0; m15 = 0; m16 = 0; m17 = 0
            m20 = 0; m21 = 0; m22 = 0; m23 = 0; m24 = 0; m25 = 0; m26 = 0; m27 = 0
            m30 = 0; m31 = 0; m32 = 0; m33 = 0; m34 = 0; m35 = 0; m36 = 0; m37 = 0
            m40 = 0; m41 = 0; m42 = 0; m43 = 0; m44 = 0; m45 = 0; m46 = 0; m47 = 0
            m50 = 0; m51 = 0; m52 = 0; m53 = 0; m54 = 0; m55 = 0; m56 = 0; m57 = 0
            m60 = 0; m61 = 0; m62 = 0; m63 = 0; m64 = 0; m65 = 0; m66 = 0; m67 = 0
            m70 = 0; m71 = 0; m72 = 0; m73 = 0; m74 = 0; m75 = 0; m76 = 0; m77 = 0
            k0 = 0; k1 = 0; k2 = 0;  k3 = 0;   k4 = 0; k5 = 0;    k6 = 0;  k7 = 0

            """Coefficients"""
            for i in range(0, 8):
                x123i = x_matrix[i][0] * x_matrix[i][1] * x_matrix[i][2]
                x1i = x_matrix[i][0]
                x2i = x_matrix[i][1]
                x3i = x_matrix[i][2]
                x12i = x1i * x2i
                x13i = x1i * x3i
                x23i = x2i * x3i
                m00 = N
                m01 += x1i
                m02 += x2i
                m03 += x3i
                m04 += x12i
                m05 += x13i
                m06 += x23i
                m07 += x123i

                m10 += x1i
                m11 += x1i ^ 2
                m12 += x12i
                m13 += x13i
                m14 += x1i ^ 2 * x2i
                m15 += x1i ^ 2 * x3i
                m16 += x123i
                m17 += x1i ^ 2 * x23i

                m20 += x2i
                m21 += x12i
                m22 += x2i ^ 2
                m23 += x23i
                m24 += x2i ^ 2 * x1i
                m25 += x123i
                m26 += x2i ^ 2 * x3i
                m27 += x2i ^ 2 * x13i

                m30 += x3i
                m31 += x13i
                m32 += x23i
                m33 += x3i ^ 2
                m34 += x123i
                m35 += x3i ^ 2 * x1i
                m36 += x3i ^ 2 * x2i
                m37 += x3i ^ 2 * x12i

                m40 += x12i
                m41 += x1i ^ 2 * x2i
                m42 += x1i * x2i ^ 2
                m43 += x123i
                m44 += x1i ^ 2 * x2i ^ 2
                m45 += x1i ^ 2 * x23i
                m46 += x13i * x2i ^ 2
                m47 += x2i ^ 2 * x1i ^ 2 * x3i

                m50 += x13i
                m51 += x1i ^ 2 * x3i
                m52 += x123i
                m53 += x1i * x3i ^ 2
                m54 += x1i ^ 2 * x23i
                m55 += x1i ^ 2 * x3i ^ 2
                m56 += x3i ^ 2 * x12i
                m57 += x3i ^ 2 * x1i ^ 2 * x2i

                m60 += x23i
                m61 += x123i
                m62 += x3i * x2i ^ 2
                m63 += x2i * x3i ^ 2
                m64 += x13i * x2i ^ 2
                m65 += x12i * x3i ^ 2
                m66 += x2i ^ 2 * x3i ^ 2
                m67 += x2i ^ 2 * x3i ^ 2 * x1i

                m70 += x123i
                m71 += x1i ^ 2 * x23i
                m72 += x13i * x2i ^ 2
                m73 += x12i * x3i ^ 2
                m74 += x1i ^ 2 * x2i ^ 2 * x3i
                m75 += x1i ^ 2 * x3i ^ 2 * x2i
                m76 += x1i * x2i ^ 2 * x3i ^ 2
                m77 += x2i ^ 2 * x3i ^ 2 * x1i ^ 2

                k0 += my_list[i]
                k1 += my_list[i] * x1i
                k2 += my_list[i] * x2i
                k3 += my_list[i] * x3i
                k4 += my_list[i] * x12i
                k5 += my_list[i] * x13i
                k6 += my_list[i] * x23i
                k7 += my_list[i] * x123i

            Det = np.linalg.det([
                [m00, m10, m20, m30, m40, m50, m60, m70],
                [m01, m11, m21, m31, m41, m51, m61, m71],
                [m02, m12, m22, m32, m42, m52, m62, m72],
                [m03, m13, m23, m33, m43, m53, m63, m73],
                [m04, m14, m24, m34, m44, m54, m64, m74],
                [m05, m15, m25, m35, m45, m55, m65, m75],
                [m06, m16, m26, m36, m46, m56, m66, m76],
                [m07, m17, m27, m37, m47, m57, m67, m77]
            ])

            b0 = np.linalg.det([
                [k0, m10, m20, m30, m40, m50, m60, m70],
                [k1, m11, m21, m31, m41, m51, m61, m71],
                [k2, m12, m22, m32, m42, m52, m62, m72],
                [k3, m13, m23, m33, m43, m53, m63, m73],
                [k4, m14, m24, m34, m44, m54, m64, m74],
                [k5, m15, m25, m35, m45, m55, m65, m75],
                [k6, m16, m26, m36, m46, m56, m66, m76],
                [k7, m17, m27, m37, m47, m57, m67, m77],
            ]) / Det

            b1 = np.linalg.det([
                [m00, k0, m20, m30, m40, m50, m60, m70],
                [m01, k1, m21, m31, m41, m51, m61, m71],
                [m02, k2, m22, m32, m42, m52, m62, m72],
                [m03, k3, m23, m33, m43, m53, m63, m73],
                [m04, k4, m24, m34, m44, m54, m64, m74],
                [m05, k5, m25, m35, m45, m55, m65, m75],
                [m06, k6, m26, m36, m46, m56, m66, m76],
                [m07, k7, m27, m37, m47, m57, m67, m77]]) / Det

            b2 = np.linalg.det([
                [m00, m10, k0, m30, m40, m50, m60, m70],
                [m01, m11, k1, m31, m41, m51, m61, m71],
                [m02, m12, k2, m32, m42, m52, m62, m72],
                [m03, m13, k3, m33, m43, m53, m63, m73],
                [m04, m14, k4, m34, m44, m54, m64, m74],
                [m05, m15, k5, m35, m45, m55, m65, m75],
                [m06, m16, k6, m36, m46, m56, m66, m76],
                [m07, m17, k7, m37, m47, m57, m67, m77]
            ]) / Det

            b3 = np.linalg.det([
                [m00, m10, m20, k0, m40, m50, m60, m70],
                [m01, m11, m21, k1, m41, m51, m61, m71],
                [m02, m12, m22, k2, m42, m52, m62, m72],
                [m03, m13, m23, k3, m43, m53, m63, m73],
                [m04, m14, m24, k4, m44, m54, m64, m74],
                [m05, m15, m25, k5, m45, m55, m65, m75],
                [m06, m16, m26, k6, m46, m56, m66, m76],
                [m07, m17, m27, k7, m47, m57, m67, m77]
            ]) / Det

            b12 = np.linalg.det([
                [m00, m10, m20, m30, k0, m50, m60, m70],
                [m01, m11, m21, m31, k1, m51, m61, m71],
                [m02, m12, m22, m32, k2, m52, m62, m72],
                [m03, m13, m23, m33, k3, m53, m63, m73],
                [m04, m14, m24, m34, k4, m54, m64, m74],
                [m05, m15, m25, m35, k5, m55, m65, m75],
                [m06, m16, m26, m36, k6, m56, m66, m76],
                [m07, m17, m27, m37, k7, m57, m67, m77]
            ]) / Det

            b13 = np.linalg.det([
                [m00, m10, m20, m30, m40, k0, m60, m70],
                [m01, m11, m21, m31, m41, k1, m61, m71],
                [m02, m12, m22, m32, m42, k2, m62, m72],
                [m03, m13, m23, m33, m43, k3, m63, m73],
                [m04, m14, m24, m34, m44, k4, m64, m74],
                [m05, m15, m25, m35, m45, k5, m65, m75],
                [m06, m16, m26, m36, m46, k6, m66, m76],
                [m07, m17, m27, m37, m47, k7, m67, m77]
            ]) / Det

            b23 = np.linalg.det([
                [m00, m10, m20, m30, m40, m50, k0, m70],
                [m01, m11, m21, m31, m41, m51, k1, m71],
                [m02, m12, m22, m32, m42, m52, k2, m72],
                [m03, m13, m23, m33, m43, m53, k3, m73],
                [m04, m14, m24, m34, m44, m54, k4, m74],
                [m05, m15, m25, m35, m45, m55, k5, m75],
                [m06, m16, m26, m36, m46, m56, k6, m76],
                [m07, m17, m27, m37, m47, m57, k7, m77]
            ]) / Det

            b123 = np.linalg.det([
                [m00, m10, m20, m30, m40, m50, m60, k0],
                [m01, m11, m21, m31, m41, m51, m61, k1],
                [m02, m12, m22, m32, m42, m52, m62, k2],
                [m03, m13, m23, m33, m43, m53, m63, k3],
                [m04, m14, m24, m34, m44, m54, m64, k4],
                [m05, m15, m25, m35, m45, m55, m65, k5],
                [m06, m16, m26, m36, m46, m56, m66, k6],
                [m07, m17, m27, m37, m47, m57, m67, k7]
            ]) / Det

            print("b\u2080:", "%.2f" % b0, " b\u2081:", "%.2f" % b1, " b\u2082:", "%.2f" % b2, " b\u2083:", "%.2f" % b3,
                  " b\u2081\u2082:", "%.2f" % b12,
                  " b\u2081\u2083:", "%.2f" % b13, " b\u2081\u2082\u2083:", "%.2f" % b123)

            print(
                f"Рівняння регресії: y = {b0:.2f}{b1:+.2f}*x\u2081{b2:+.2f}*x\u2082{b3:+.2f}*x\u2083{b12:+.2f}*x\u2081\u2082{b13:+.2f}*x\u2081\u2083{b123:+.2f}*x\u2081\u2082\u2083")

            # find dispersion
            S2 = []
            for i in range(len(y_list)):
                S2.append(((y_list[i][0] - my_list[i]) ** 2 + (y_list[i][1] - my_list[i]) ** 2 + (
                        y_list[i][2] - my_list[i]) ** 2) / 3)

            """KOHREN"""
            Gp = max(S2) / sum(S2)
            print("Gp:", Gp)
            m = len(y_list[0])
            f1 = m - 1
            f2 = N  # N=8
            q = 0.05

            Gt = [None, 0.68, 0.516, 0.438, 0.391, 0.3595, 0.3365, 0.3185, 0.3043, 0.2926, 0.2829, 0.2462, 0.2022,
                  0.1616,
                  0.1250]
            print("Gt:", Gt[f1])

            if Gp < Gt[f1]:
                print("Дисперсія однорідна")
                break
            else:
                print("Дисперсія не однорідна")
                m += 1
                for obj in y_list:
                    obj.append(random.randint(int(Y_MIN), int(Y_MAX)))

        x_matrix_normal = [
            [1, -1, -1, -1],
            [1, -1, -1, 1],
            [1, -1, 1, -1],
            [1, -1, 1, 1],
            [1, 1, -1, -1],
            [1, 1, -1, 1],
            [1, 1, 1, -1],
            [1, 1, 1, 1],
        ]

        """STUDENT"""


        def getBeta(i):
            sum = 0
            for j in range(N):
                sum += my_list[j] * x_matrix_normal[j][i]
            sum /= N
            return sum


        S2B = sum(S2) / N
        S2beta = S2B / (N * m)
        Sbeta = np.sqrt(S2beta)

        beta0 = getBeta(0)
        beta1 = getBeta(1)
        beta2 = getBeta(2)
        beta3 = getBeta(3)

        t_criterion = []
        t_criterion.append(abs(beta0) / Sbeta, )
        t_criterion.append(abs(beta1) / Sbeta)
        t_criterion.append(abs(beta2) / Sbeta)
        t_criterion.append(abs(beta3) / Sbeta)

        t0 = abs(beta0) / Sbeta
        t1 = abs(beta1) / Sbeta
        t2 = abs(beta2) / Sbeta
        t3 = abs(beta3) / Sbeta

        f3 = f1 * f2

        t_tab = t.ppf((1 + (1 - q)) / 2, f3)
        print("t табличне:", t_tab)
        if t0 < t_tab:
            b0 = 0
            print("t\u2080:", t0, " t0<t_tab; b0=0")
        if t1 < t_tab:
            b1 = 0
            print("t\u2081:", t1, " t\u2081<t_tab; b\u2081=0")
        if t2 < t_tab:
            b2 = 0
            print("t\u2082:", t2, " t\u2082<t_tab; b\u2082=0")
        if t3 < t_tab:
            b3 = 0
            print("t\u2083:", t3, " t\u2083<t_tab; b\u2083=0")

        y_hat = []
        for i in range(N):
            y_hat.append(
                b0 + b1 * x_matrix[i][0] + b2 * x_matrix[i][1] + b3 * x_matrix[i][2] + b12 * x_matrix[i][0] *
                x_matrix[i][1] +
                b13 * x_matrix[i][0] * x_matrix[i][2] + b123 * x_matrix[i][0] * x_matrix[i][1] * x_matrix[i][2])

            print(
                f"^y{chr(8321 + i)} = {b0:.2f}{b1:+.2f}*x{chr(8321 + i)}\u2081{b2:+.2f}*x{chr(8321 + i)}\u2082{b3:+.2f}*x{chr(8321 + i)}\u2083{b12:+.2f}*x{chr(8321 + i)}\u2081"
                f"*x{chr(8321 + i)}\u2082{b13:+.2f}*x{chr(8321 + i)}\u2081*x{chr(8321 + i)}\u2083{b123:+.2f}*x{chr(8321 + i)}\u2081*x{chr(8321 + i)}\u2082*x{chr(8321 + i)}\u2083 "
                f"= {y_hat[i]:.2f}")

        """FISHER"""

        d = 2
        f4 = N - d
        S2_ad = 0
        for i in range(N):
            S2_ad += (m / (N - d) * ((y_hat[i] - my_list[i]) ** 2))

        Fp = S2_ad / S2B
        Ft = f.ppf(1 - q, f4, f3)
        print("Fp:", Fp)
        print("Ft:", Ft)
        if Fp > Ft:
            print("Рівняння регресії не адекватно оригіналу при рівні значимості 0,05")
            print("Щоб почати з початку введіть довільний символ в консоль\n")
            if input() == 0:
                break
        else:
            print("Рівняння регресії адекватно оригіналу при рівні значимості 0,05")
            break