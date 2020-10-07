from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pylab
import matplotlib.pyplot as plt
import matplotlib.patches as ptc
from celluloid import Camera
from interval_analysis import calculate_area_ellipse_circle_intersect, canonical_form
from mpmath import *

#               M
#              /
#             /
#        l_3 /
#           /
#        B /
#         /|
#    l_3 / |
#       /  | l_4
#    A |   |
# l_2  |   |
#    O  l_1  C
# l_1 : l_2 : l_3 : l_4 = 2 : 1 : 2.5 : 2.5


class LambdaMechanism:
    def __init__(self, size, phi, delta_axis_normal, delta_axis_angle, Rb, Hb, Rbh):
        self.l1 = 2 * size
        self.l2 = size
        self.l3 = 2.5 * size
        self.l4 = 2.5 * size
        self.phi = phi
        self.delta_axis_normal = delta_axis_normal
        self.delta_axis_angle = delta_axis_angle
        self.delta_B_normal = None
        self.delta_B_angle = None
        self.Rb = Rb
        self.Hb = Hb
        self.Rbh = Rbh
        self.Rh = Rb + Rbh
        self.A = None
        self.M = None
        self.B = None
        self.cylinder = None
        self.cylinder_rot = None
        self.cylinder_rot_out = None

    def calculate_lambda_motion(self):
        rot = self.__rotation(self.delta_axis_normal, self.delta_axis_angle)
        A = (rot * self.l2).dot(np.array([np.cos(self.phi), np.sin(self.phi), np.zeros(len(self.phi))]))
        AB_xy = np.sqrt(self.l3 ** 2 - A[2, :] ** 2)
        AC = np.sqrt((self.l1 - A[0, :]) ** 2 + A[1, :] ** 2)
        angle_ACB = np.arccos((AC ** 2 + self.l4 ** 2 - AB_xy ** 2) / (2 * AC * self.l4))
        angle_OCA = np.arctan((- A[1, :]) / (self.l1 - A[0, :]))
        angle_OCB = angle_ACB - angle_OCA
        B = np.array([
            self.l1 - self.l4 * np.cos(angle_OCB),
            self.l4 * np.sin(angle_OCB),
            np.zeros(len(self.phi))
        ])
        M = 2 * B - A
        BA = B - A
        n = np.sqrt(BA[0, :] ** 2 + BA[1, :] ** 2)
        self.delta_B_normal = np.array([- BA[1, :] / n, BA[0, :] / n, np.zeros(len(self.phi))])
        self.delta_B_angle = np.array(np.arctan(BA[2, :] / n))
        self.A = A
        self.B = B
        self.M = M

    def draw_lambda_motion_2d(self):
        fig, ax = plt.subplots()
        ax.set_title("Animation of Chebyshev lambda mechanism motion. View 2D \n "
                     "Axis incline normal n = " + str(self.delta_axis_normal) +
                     ", Axis incline angle grad = " + "{:.4f}".format(self.delta_axis_angle))
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.axis('equal')
        camera = Camera(fig)

        colors = ['tab:red', 'tab:blue', 'tab:orange', 'tab:gray']

        for i in range(len(self.phi)):
            ax.grid('on')
            # two fixed motionless hinges
            ax.plot([0, self.l1], [0, 0], "^", color=colors[0], markersize=10)
            # the circle along which the point B moves
            circle = ptc.Circle((self.l1, 0), radius=self.l4, fill=None, color=colors[3], linestyle='--', linewidth=0.5)
            ax.add_patch(circle)
            # trajectories of points A, M, B
            ax.plot(self.A[0, :], self.A[1, :], color=colors[1], linewidth=1)  # circle
            ax.plot(self.M[0, :], self.M[1, :], color=colors[1], linewidth=1)  # "mushroom cap"
            ax.plot(self.B[0, :], self.B[1, :], color=colors[1], linewidth=1)  # part of a circle
            # parts of the lambda mechanism
            ax.plot([self.A[0, i], self.M[0, i]], [self.A[1, i], self.M[1, i]], color=colors[2], linewidth=2)  # rod
            ax.plot([self.A[0, i], 0], [self.A[1, i], 0], colors[2], linewidth=2)  # crank
            ax.plot([self.B[0, i], self.l1], [self.B[1, i], 0], colors[2], linewidth=2)  # rocker
            # points A, M, B
            ax.plot([self.A[0, i], self.B[0, i], self.M[0, i]], [self.A[1, i], self.B[1, i], self.M[1, i]],
                    "o", color=colors[0], linewidth=1)
            camera.snap()

        anim = camera.animate()
        anim.save("lambda_animation_2d.gif", writer="imagemagick")
        plt.close(fig)

    def draw_lambda_motion_3d(self):
        fig = pylab.figure()
        ax = Axes3D(fig)
        ax.set_title("Animation of Chebyshev lambda mechanism motion. View 3D \n "
                     "Axis incline normal n = " + str(self.delta_axis_normal) +
                     ", Axis incline angle grad = " + "{:.4f}".format(self.delta_axis_angle))
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        # ax.view_init(elev=-30, azim=0)
        camera = Camera(fig)

        colors = ['tab:red', 'tab:blue', 'tab:orange', 'tab:gray']

        for i in range(len(self.phi)):
            # two fixed motionless hinges
            ax.plot([0, self.l1], [0, 0], [0, 0], "^", color=colors[0], markersize=10)
            # the circle along which the point B moves
            theta = np.linspace(0, 2 * np.pi, 201)
            ax.plot(self.l4 * np.cos(theta) + self.l1, self.l4 * np.sin(theta), np.zeros(len(theta)), color=colors[3],
                    linestyle='--', linewidth=0.5)
            # trajectories of points A, M, B
            ax.plot(self.A[0, :], self.A[1, :], self.A[2, :], color=colors[1], linewidth=1)  # circle
            ax.plot(self.M[0, :], self.M[1, :], self.M[2, :], color=colors[1], linewidth=1)  # "mushroom cap"
            ax.plot(self.B[0, :], self.B[1, :], self.B[2, :], color=colors[1], linewidth=1)  # part of a circle
            # parts of the lambda mechanism
            ax.plot([self.A[0, i], self.M[0, i]], [self.A[1, i], self.M[1, i]], [self.A[2, i], self.M[2, i]],
                    color=colors[2], linewidth=2)  # rod
            ax.plot([self.A[0, i], 0], [self.A[1, i], 0], [self.A[2, i], 0], colors[2], linewidth=2)  # crank
            ax.plot([self.B[0, i], self.l1], [self.B[1, i], 0], [self.B[2, i], 0], colors[2], linewidth=2)  # rocker
            # points A, M, B
            ax.plot([self.A[0, i], self.B[0, i], self.M[0, i]], [self.A[1, i], self.B[1, i], self.M[1, i]],
                    [self.A[2, i], self.B[2, i], self.M[2, i]], "o", color=colors[0], linewidth=1)

            camera.snap()

        anim = camera.animate()
        anim.save("lambda_animation_3d.gif", writer="imagemagick")
        plt.close(fig)

    def calculate_bearing_motion(self):
        h = np.linspace(- self.Hb / 2, self.Hb / 2, 30)
        self.cylinder = np.array([self.Rb * np.kron(np.cos(self.phi), np.ones((1, len(h))))[0],
                             self.Rb * np.kron(np.sin(self.phi), np.ones((1, len(h))))[0],
                             np.kron(np.ones((1, len(self.phi))), h)[0]])
        self.cylinder_rot = []
        self.cylinder_rot_out = []
        for i in range(len(self.delta_B_angle)):
            rotB = self.__rotation(self.delta_B_normal[:, i], self.delta_B_angle[i])
            curr_rot = rotB.dot(self.cylinder)
            curr_rot_out = [[], [], []]
            for j in range(len(curr_rot[0])):
                if curr_rot[0][j] ** 2 + curr_rot[1][j] ** 2 > self.Rh ** 2:
                    curr_rot_out[0].append(curr_rot[0][j])
                    curr_rot_out[1].append(curr_rot[1][j])
                    curr_rot_out[2].append(curr_rot[2][j])
            self.cylinder_rot_out.append(np.array(curr_rot_out))
            self.cylinder_rot.append(np.array(curr_rot))

    def draw_bearing_animation_3d(self):
        fig = pylab.figure()
        ax = Axes3D(fig)
        ax.set_title("Animation of Bearing mechanism. View 3D \n" + "Rb = 17.5, Hb = 10, Rbh = 0.03 \n" +
                     "Axis incline normal n = " + str(self.delta_axis_normal) +
                     ", Axis incline angle grad = " + "{:.4f}".format(self.delta_axis_angle))
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        camera = Camera(fig)
        colors = ['tab:orange', 'tab:blue', 'tab:red']

        for i in range(len(self.delta_B_angle)):
            BC = np.array([self.l1, 0, 0]).T - self.B[:, i]
            BC = BC * (self.Rb / np.sqrt(np.linalg.norm(BC.T.dot(BC))))
            BA = self.A[:, i] - self.B[:, i]
            BA = self.Rb * BA / np.linalg.norm(BA)

            ax.plot(self.cylinder_rot_out[i][0, :], self.cylinder_rot_out[i][1, :], self.cylinder_rot_out[i][2, :],
                    ".", color=colors[2])
            ax.plot([0, 0], [0, 0], [- self.Hb / 2, self.Hb / 2], color=colors[1])
            ax.plot([0, BC[0]], [0, BC[1]], [- self.Hb / 2, - self.Hb / 2 + BC[2]], color=colors[1])
            ax.plot([0, BA[0]], [0, BA[1]], [- self.Hb / 2, - self.Hb / 2 + BA[2]], color=colors[0])
            ax.plot([0, BC[0]], [0, BC[1]], [self.Hb / 2, self.Hb / 2 + BC[2]], color=colors[1])
            ax.plot([0, BA[0]], [0, BA[1]], [self.Hb / 2, self.Hb / 2 + BA[2]], color=colors[0])
            ax.plot(self.Rb * np.cos(self.phi), self.Rb * np.sin(self.phi), np.zeros(len(self.phi)) - self.Hb / 2,
                    color=colors[1])
            ax.plot(self.Rb * np.cos(self.phi), self.Rb * np.sin(self.phi), np.zeros(len(self.phi)) + self.Hb / 2,
                    color=colors[1])
            # ax.text(-20, -10, -2, "phi = " + "{:.4f}".format(self.phi[i]))
            camera.snap()

        anim = camera.animate()
        anim.save("Bearing_animation_3d.gif", writer="imagemagick")
        plt.close(fig)

    def draw_number_of_jamming_points(self):
        fig, ax = plt.subplots()
        ax.set_title("Dependence of the number of cylinder jamming points on the crank angle \n "
                     "Rb = 17.5, Hb = 10, Rbh = 0.03 \n" +
                     "Axis incline normal n = " + str(self.delta_axis_normal) +
                     ", Axis incline angle grad = " + "{:.4f}".format(self.delta_axis_angle))
        ax.set_xlabel('Crank angle, phi, grad')
        ax.set_ylabel('Area, rel.units')

        x = []
        y = []
        for i in range(len(self.phi)):
            angle = self.phi[i] * 180 / np.pi
            x.append(angle)
            y.append(len(self.cylinder_rot_out[i][0]))

        x = np.array(x)
        y = np.array(y)
        ax.grid('on')
        ax.plot(x, y, '-', color='tab:blue')
        ax.plot([x[0], x[-1]], [y[0], y[-1]], 'o', color='tab:red')
        ax.text(80, 450, 'Start of motion')
        ax.text(390, 450, 'End of motion')
        plt.savefig('number_jamm_points.png')
        plt.close(fig)

    def calculate_area_of_intersection(self):
        h = self.Hb / 2
        area_array = []
        for i in range(0, len(self.phi)):
            print("i = ", i)
            # area = 0
            area = self.do_something(i, h)
            area_array.append(area)
        v = 1
        return area_array

    def draw_area(self, area):
        fig, ax = plt.subplots()
        ax.set_title("Dependence of area on the crank angle \n "
                     "Rb = 17.5, Hb = 10, Rbh = 0.03 \n" +
                     "Axis incline normal n = " + str(self.delta_axis_normal) +
                     ", Axis incline angle grad = " + "{:.4f}".format(self.delta_axis_angle))
        ax.set_xlabel('Crank angle, phi, grad')
        ax.set_ylabel('Area, mm^2')
        ax.grid('on')

        for i in range(0, len(self.phi)):
            angle = self.phi[i] * 180 / np.pi
            ax.plot(angle, area[i][0], '^', color='tab:blue')
            ax.plot(angle, area[i][1], '^', color='tab:blue')

        plt.savefig('area.png')
        plt.close(fig)



    def do_something(self, i, h):
        # radius of the "ideal" circle
        R1 = self.Rb
        # radius of the circle to be transformed
        R2 = self.Rh
        circle1 = np.array([R1 * np.cos(self.phi), R1 * np.sin(self.phi), np.zeros(len(self.phi)) + h])
        circle2 = np.array([R2 * np.cos(self.phi), R2 * np.sin(self.phi), np.zeros(len(self.phi)) + h])
        rotB = self.__rotation(self.delta_B_normal[:, i], self.delta_B_angle[i])
        circle2_rot = rotB.dot(circle2)
        a, b = canonical_form(circle2_rot, rotB)
        print("a = ", a, ", b = ", b)
        if a > self.Rb and b > self.Rb:
            area = iv.matrix([mpf(0.0), mpf(0.0)])
            return [mpf(area[0]), mpf(area[1])]
        area = calculate_area_ellipse_circle_intersect(R1, a, b, 0.0001)
        return area





    def __rotation(self, n, alpha):
        c = np.cos(alpha)
        s = np.sin(alpha)

        rot = np.array([[
            c + (n[0] ** 2) * (1 - c), n[0] * n[1] * (1 - c) + n[2] * s, n[0] * n[2] * (1 - c) - n[1] * s
        ], [
            n[0] * n[1] * (1 - c) - n[2] * s, c + (n[1] ** 2) * (1 - c), n[0] * n[2] * (1 - c) + n[0] * s
        ], [
            n[0] * n[2] * (1 - c) + n[1] * s, n[1] * n[2] * (1 - c) - n[0] * s, c + (n[2] ** 2) * (1 - c)
        ]])

        return rot.T
