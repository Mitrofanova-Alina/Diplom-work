import numpy as np
import matplotlib.pyplot as plt
from mpmath import *


def calculate_area_ellipse_circle_intersect(r, a, b, eps):
    if eps < 0:
        raise ValueError(f"'{eps}' is incorrect value for epsilon. It must be greater than or equal to 0")
    if a <= 0 or b <= 0:
        raise ValueError(f"'{a}' or '{b}' is incorrect value for a or b. It must be greater than 0")
    if a < b:
        a, b = b, a
    interval_analysis = IntervalAnalysis()
    n = 500
    x = np.linspace(0, a + eps * 5, n)
    y = np.linspace(0, r + eps * 5, n)
    X_tmp = []
    for i in range(len(x) - 1):
        for j in range(len(y) - 1):
            X_0 = iv.matrix([[mpi(x[i], x[i + 1])], [mpi(y[j], y[j + 1])]])
            X = interval_analysis.nonlinear_kravchik(r, a, b, eps, X_0)
            if X[0]:
                X_tmp.append(X)
                # print("(", i, ", ", j, "): ", X)
    if len(X_tmp) < 1:
        area = iv.matrix([mpf(0.0), mpf(0.0)])
        return [mpf(area[0]), mpf(area[1])]
    x1, x2, y1, y2 = [], [], [], []
    for i in X_tmp:
        x1.append(mpf(i[0].a))
        x2.append(mpf(i[0].b))
        y1.append(mpf(i[1].a))
        y2.append(mpf(i[1].b))
    answer = iv.matrix([[mpi(min(x1), max(x2))], [mpi(min(y1), max(y2))]])
    print("Answer X = ", answer)
    points = create_interval_points(r, a, b, eps, answer)
    area = shoelace_formula(points)
    print("Area S = ", area)
    # interval_analysis.draw_circle_and_ellipse(r, a, b, eps, X_tmp, answer)
    return area


def create_interval_points(r, a, b, eps, answer):
    n = 500
    phi = np.linspace(- np.pi / 2, np.pi / 2, n)
    circle1 = np.array([(r - eps) * np.cos(phi), (r - eps) * np.sin(phi)])
    circle2 = np.array([(r + eps) * np.cos(phi), (r + eps) * np.sin(phi)])
    ellipse1 = np.array([(a - eps) * np.cos(phi), (b - eps) * np.sin(phi)])
    ellipse2 = np.array([(a + eps) * np.cos(phi), (b + eps) * np.sin(phi)])
    points_circle = []
    points_ellipse = []
    find_first_circle_intersect = False
    find_second_circle_intersect = False
    find_first_ellipse_intersect = False
    find_second_ellipse_intersect = False
    curr_circle_intersect = None
    curr_ellipse_intersect = None
    for i in range(0, len(phi) - 1, 2):
        circle_11, circle_12, circle_21, circle_22 = circle1[:, i], circle1[:, i + 1], circle2[:, i], circle2[:, i + 1]
        circle_x1 = min([circle_11[0], circle_12[0], circle_21[0], circle_22[0]])
        circle_x2 = max([circle_11[0], circle_12[0], circle_21[0], circle_22[0]])
        circle_y1 = min([circle_11[1], circle_12[1], circle_21[1], circle_22[1]])
        circle_y2 = max([circle_11[1], circle_12[1], circle_21[1], circle_22[1]])
        ellipse_11, ellipse_12, ellipse_21, ellipse_22 = ellipse1[:, i], ellipse1[:, i + 1], ellipse2[:, i], ellipse2[:, i + 1]
        ellipse_x1 = min([ellipse_11[0], ellipse_12[0], ellipse_21[0], ellipse_22[0]])
        ellipse_x2 = max([ellipse_11[0], ellipse_12[0], ellipse_21[0], ellipse_22[0]])
        ellipse_y1 = min([ellipse_11[1], ellipse_12[1], ellipse_21[1], ellipse_22[1]])
        ellipse_y2 = max([ellipse_11[1], ellipse_12[1], ellipse_21[1], ellipse_22[1]])
        # find circle in second answer
        x1, x2, y1, y2 = mpf(answer[0].a), mpf(answer[0].b), - mpf(answer[1].b), - mpf(answer[1].a)
        if ellipse_x1 >= x1 and ellipse_x2 <= x2 and ellipse_y1 >= y1 and ellipse_y2 <= y2:
            curr_ellipse_intersect = iv.matrix([[mpi(ellipse_x1, ellipse_x2)], [mpi(ellipse_y1, ellipse_y2)]])
        elif len(points_ellipse) == 0 and curr_ellipse_intersect:
            points_ellipse.append(iv.matrix([curr_ellipse_intersect[0], curr_ellipse_intersect[1]]))
            find_first_ellipse_intersect = True
        if circle_x1 >= x1 and circle_x2 <= x2 and circle_y1 >= y1 and circle_y2 <= y2:
            curr_circle_intersect = iv.matrix([[mpi(circle_x1, circle_x2)], [mpi(circle_y1, circle_y2)]])
        elif len(points_circle) == 0 and curr_circle_intersect:
            points_circle.append(iv.matrix([curr_circle_intersect[0], curr_circle_intersect[1]]))
            find_first_circle_intersect = True

        x1, x2, y1, y2 = mpf(answer[0].a), mpf(answer[0].b), mpf(answer[1].a), mpf(answer[1].b)
        if ellipse_x1 >= x1 and ellipse_x2 <= x2 and ellipse_y1 >= y1 and ellipse_y2 <= y2 and not find_second_ellipse_intersect:
            points_ellipse.append(iv.matrix([[mpi(ellipse_x1, ellipse_x2)], [mpi(ellipse_y1, ellipse_y2)]]))
            find_second_ellipse_intersect = True
        if circle_x1 >= x1 and circle_x2 <= x2 and circle_y1 >= y1 and circle_y2 <= y2 and not find_second_circle_intersect:
            points_circle.append(iv.matrix([[mpi(circle_x1, circle_x2)], [mpi(circle_y1, circle_y2)]]))
            find_second_circle_intersect = True

        if find_first_ellipse_intersect and not find_second_ellipse_intersect:
            points_ellipse.append(iv.matrix([[mpi(ellipse_x1, ellipse_x2)], [mpi(ellipse_y1, ellipse_y2)]]))
        if find_first_circle_intersect and not find_second_circle_intersect:
            points_circle.append(iv.matrix([[mpi(circle_x1, circle_x2)], [mpi(circle_y1, circle_y2)]]))

    # fig, ax = plt.subplots()
    # ax.set_title("Intersection of circle and ellipse, eps = " + str(eps) + "\n Circle r = " + str(
    #     r) + ", Ellipse a = " + str(a) + ", b = " + str(b))
    # ax.set_xlabel('x')
    # ax.set_ylabel('y')
    # ax.axis('equal')
    # ax.grid('on')
    # ax.plot(circle1[0, :], circle1[1, :], color='tab:green')
    # ax.plot(circle2[0, :], circle2[1, :], color='tab:green')
    # ax.plot(ellipse1[0, :], ellipse1[1, :], color='tab:orange')
    # ax.plot(ellipse2[0, :], ellipse2[1, :], color='tab:orange')
    # for i in points_circle:
    #     circle_x1, circle_x2, circle_y1, circle_y2 = mpf(i[0].a), mpf(i[0].b), mpf(i[1].a), mpf(i[1].b)
    #     ax.plot([circle_x1, circle_x2, circle_x1, circle_x1, circle_x1, circle_x2, circle_x2, circle_x2], [circle_y1, circle_y1, circle_y1, circle_y2, circle_y2, circle_y2, circle_y1, circle_y2], color='tab:red')
    # for i in points_ellipse:
    #     circle_x1, circle_x2, circle_y1, circle_y2 = mpf(i[0].a), mpf(i[0].b), mpf(i[1].a), mpf(i[1].b)
    #     ax.plot([circle_x1, circle_x2, circle_x1, circle_x1, circle_x1, circle_x2, circle_x2, circle_x2], [circle_y1, circle_y1, circle_y1, circle_y2, circle_y2, circle_y2, circle_y1, circle_y2], color='tab:blue')
    # circle_x1, circle_x2, circle_y1, circle_y2 = mpf(answer[0].a), mpf(answer[0].b), mpf(answer[1].a), mpf(answer[1].b)
    # ax.plot([circle_x1, circle_x2, circle_x1, circle_x1, circle_x1, circle_x2, circle_x2, circle_x2], [circle_y1, circle_y1, circle_y1, circle_y2, circle_y2, circle_y2, circle_y1, circle_y2], color='tab:cyan')
    # ax.plot([circle_x1, circle_x2, circle_x1, circle_x1, circle_x1, circle_x2, circle_x2, circle_x2], [-circle_y1, -circle_y1, -circle_y1, -circle_y2, -circle_y2, -circle_y2, -circle_y1, -circle_y2], color='tab:cyan')
    # # ax.plot([-circle_x1, -circle_x2, -circle_x1, -circle_x1, -circle_x1, -circle_x2, -circle_x2, -circle_x2], [-circle_y1, -circle_y1, -circle_y1, -circle_y2, -circle_y2, -circle_y2, -circle_y1, -circle_y2], color='tab:cyan')
    # # ax.plot([-circle_x1, -circle_x2, -circle_x1, -circle_x1, -circle_x1, -circle_x2, -circle_x2, -circle_x2], [circle_y1, circle_y1, circle_y1, circle_y2, circle_y2, circle_y2, circle_y1, circle_y2], color='tab:cyan')
    # ax.set_xlim((0, 7))
    # ax.set_ylim((-6, 6))
    # fig.show()
    # plt.close(fig)
    points_circle.reverse()
    return points_ellipse + points_circle



# A = 1 / 2 * |sum1 + x_n * y_1 - sum2 - x_1 * y_n|
# sum1 = sum_{1, n-1} (x_i * y_{i+1})
# sum2 = sum_{1, n-1} (y_i * x_{i+1})
def shoelace_formula(points):
    n = len(points)
    area = points[-1][0] * points[0][1] - points[0][0] * points[-1][1]
    for i in range(0, n - 1):
        area += points[i][0] * points[i + 1][1] - points[i + 1][0] * points[i][1]
    return [mpf(abs(area.a)), mpf(abs(area.b))]

def canonical_form(ellipse, rot):
    A = rot[0][0] ** 2 + rot[0][1] ** 2
    C = rot[1][0] ** 2 + rot[1][1] ** 2
    B = rot[0][0] * rot[0][1] + rot[1][0] * rot[1][1]
    if A - C == 0:
        alpha = np.pi / 4
    else:
        alpha = 1 / 2 * np.arctan((2 * B) / (A - C))
    # print(alpha)
    new_ellipse = np.array([ellipse[0, :] * np.cos(alpha) + ellipse[1, :] * np.sin(alpha),
                            - ellipse[0, :] * np.sin(alpha) + ellipse[1, :] * np.cos(alpha)])
    a = max(new_ellipse[0])
    b = max(new_ellipse[1])
    if a < b:
        a, b = b, a
    return a, b

class IntervalAnalysis:
    # Task is:
    #   x^2 + y^2 = r^2 + [-eps, eps]
    #   x^2 / a^2 + y^2 / b^2 = 1 + [-eps, eps]
    # Then for method Kravchik we have:
    #   ( x^2 + y^2 - r^2 - [-eps, eps]           )( x_1 ) = ( 0 )
    #   ( x^2 / a^2 + y^2 / b^2 - 1 - [-eps, eps] )( x_2 ) = ( 0 )
    def nonlinear_kravchik(self, r, a, b, eps, X):
        I = iv.matrix([[mpi(1., 1.), mpi(0., 0.)], [mpi(0., 0.), mpi(1., 1.)]])
        X1 = iv.matrix([[mpi(0, 0)], [mpi(0, 0)]])
        tol = 1e-16
        num_iter = 0
        var_lambda = self.__lambda(X, r, a, b)
        var_bracket = I - var_lambda * self.__jacobian(X, r, a, b)
        while 1:
            num_iter += 1
            count = 0
            if abs(iv.det(self.__jacobian(self.__mid(X), r, a, b))) < 0.000001:
                X1 = [[], []]
                return X1
            F_value = self.__F(self.__mid(X), r, a, b, eps)
            tmp2 = X - self.__mid(X)
            tmp3 = var_lambda * F_value
            tmp4 = var_bracket * tmp2
            X1WithoutIntersection = self.__mid(X) - tmp3 + tmp4

            for i in range(0, 2):
                if mpf(X1WithoutIntersection[i].a) > mpf(X[i].b) or mpf(X1WithoutIntersection[i].b) < mpf(X[i].a):
                    X1 = [[], []]
                    return X1
                else:
                    inf = mpf(X[i].a) if mpf(X1WithoutIntersection[i].a) < mpf(X[i].a) else mpf(
                        X1WithoutIntersection[i].a)
                    sup = mpf(X[i].b) if mpf(X1WithoutIntersection[i].b) > mpf(X[i].b) else mpf(
                        X1WithoutIntersection[i].b)
                    X1[i] = mpi(inf, sup)

            for i in range(0, 2):
                if abs(mpf(X1[i].a) - mpf(X[i].a)) < tol and abs(mpf(X1[i].b) - mpf(X[i].b)) < tol:
                    count += 1

            if count == 2:
                X = iv.matrix([X1[0], X1[1]])
                break
            X = iv.matrix([X1[0], X1[1]])
        # print("num iter = ", num_iter)
        return X

    def __F(self, X, r, a, b, eps):
        return iv.matrix([[X[0] ** 2 + X[1] ** 2 - r ** 2 - mpi(- eps, eps)],
                          [X[0] ** 2 / a ** 2 + X[1] ** 2 / b ** 2 - 1 - mpi(- eps, eps)]])

    def __jacobian(self, X, r, a, b):
        return iv.matrix([[2 * X[0], 2 * X[1]],
                          [2 * X[0] / a ** 2, 2 * X[1] / b ** 2]])

    def __lambda(self, X, r, a, b):
        J = self.__jacobian(self.__mid(X), r, a, b)
        return J ** -1

    def __mid(self, X):
        return iv.matrix([[X[0].mid], [X[1].mid]])

    def draw_circle_and_ellipse(self, r, a, b, eps, X_tmp, answer):
        fig, ax = plt.subplots()
        ax.set_title("Intersection of circle and ellipse, eps = " + str(eps) + "\n Circle r = " + str(
            r) + ", Ellipse a = " + str(a) + ", b = " + str(b))
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.axis('equal')
        ax.grid('on')
        phi = np.linspace(0, 2 * np.pi, 100)
        # circle = np.array([r * np.cos(phi), r * np.sin(phi)])
        circle1 = np.array([(r - eps) * np.cos(phi), (r - eps) * np.sin(phi)])
        circle2 = np.array([(r + eps) * np.cos(phi), (r + eps) * np.sin(phi)])
        ellipse1 = np.array([(a - eps) * np.cos(phi), (b - eps) * np.sin(phi)])
        ellipse2 = np.array([(a + eps) * np.cos(phi), (b + eps) * np.sin(phi)])
        ax.plot(circle1[0, :], circle1[1, :], color='tab:green')
        ax.plot(circle2[0, :], circle2[1, :], color='tab:green')
        ax.plot(ellipse1[0, :], ellipse1[1, :], color='tab:orange')
        ax.plot(ellipse2[0, :], ellipse2[1, :], color='tab:orange')
        colors = ['tab:brown', 'tab:pink', 'tab:purple', 'tab:gray', 'tab:olive', 'tab:cyan', 'tab:green', 'tab:blue',
                  'tab:red', 'tab:orange']
        j = 0
        # for i in X_tmp:
        #     x1 = mpf(i[0].a)
        #     x2 = mpf(i[0].b)
        #     y1 = mpf(i[1].a)
        #     y2 = mpf(i[1].b)
        #     ax.plot([x1, x2], [y1, y1], color=colors[j])
        #     ax.plot([x1, x1], [y1, y2], color=colors[j])
        #     ax.plot([x1, x2], [y2, y2], color=colors[j])
        #     ax.plot([x2, x2], [y1, y2], color=colors[j])
        #     j += 1
        x1, x2, y1, y2 = mpf(answer[0].a), mpf(answer[0].b), mpf(answer[1].a), mpf(answer[1].b)
        ax.plot([x1, x2, x1, x1, x1, x2, x2, x2], [y1, y1, y1, y2, y2, y2, y1, y2], color='tab:red')
        ax.plot([x1, x2, x1, x1, x1, x2, x2, x2], [-y1, -y1, -y1, -y2, -y2, -y2, -y1, -y2], color='tab:red')
        ax.plot([-x1, -x2, -x1, -x1, -x1, -x2, -x2, -x2], [-y1, -y1, -y1, -y2, -y2, -y2, -y1, -y2], color='tab:red')
        ax.plot([-x1, -x2, -x1, -x1, -x1, -x2, -x2, -x2], [y1, y1, y1, y2, y2, y2, y1, y2], color='tab:red')
        ax.set_xlim((3.5, 4.5))
        ax.set_ylim((2.5, 3.5))
        fig.show()
        plt.close(fig)
