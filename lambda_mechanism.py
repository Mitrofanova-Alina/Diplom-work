from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pylab
import matplotlib.pyplot as plt
import matplotlib.patches as ptc
import matplotlib.ticker as ticker
from celluloid import Camera


class LambdaMechanism:
    def __init__(self, size, phi, delta_axis_normal, delta_axis_angle):
        self.l1 = 2 * size
        self.l2 = size
        self.l3 = 2.5 * size
        self.l4 = 2.5 * size
        self.phi = phi
        self.delta_axis_normal = delta_axis_normal
        self.delta_axis_angle = delta_axis_angle

    def motion(self):
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
        delta_B_normal = np.array([- BA[1, :] / n, BA[0, :] / n, np.zeros(len(self.phi))])
        delta_B_angle = np.array(np.arctan(BA[2, :] / n))
        return A, B, M, delta_B_normal, delta_B_angle

    def draw_lambda_mechanism_motion_2d(self, A, M, B):
        fig, ax = plt.subplots()
        ax.set_title("Lambda mechanism 2d animation")
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.axis('equal')
        camera = Camera(fig)

        for i in range(len(self.phi)):
            ax.grid('on')
            ax.plot([0], [0], "o", color='tomato', linewidth=1)
            ax.plot([self.l1], [0], "o", color='tomato', linewidth=1)
            circle = ptc.Circle((self.l1, 0), radius=self.l4, fill=None, color='grey', linestyle='--', linewidth=0.5)
            ax.add_patch(circle)
            ax.plot(A[0, :], A[1, :], color="lightseagreen", linewidth=0.5)
            ax.plot(M[0, :], M[1, :], color="lightseagreen", linewidth=0.5)
            ax.plot(B[0, :], B[1, :], color="lightseagreen", linewidth=0.5)
            ax.plot([A[0, i], M[0, i]], [A[1, i], M[1, i]], color="royalblue", linewidth=2)
            ax.plot([A[0, i], 0], [A[1, i], 0], "royalblue", linewidth=2)
            ax.plot([B[0, i], self.l1], [B[1, i], 0], "royalblue", linewidth=2)
            ax.plot([A[0, i]], [A[1, i]], "o", color='grey', linewidth=1)
            ax.plot([B[0, i]], [B[1, i]], "o", color='grey', linewidth=1)
            ax.plot([M[0, i]], [M[1, i]], "o", color='grey', linewidth=1)
            camera.snap()

        anim = camera.animate()
        anim.save("lambda animation 2d.gif", writer="imagemagick")
        plt.close(fig)

    def draw_lambda_mechanism_motion_3d(self, A, M, B):
        fig = pylab.figure()
        ax = Axes3D(fig)
        ax.set_title("Lambda mechanism 3d animation")
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        camera = Camera(fig)

        for i in range(len(self.phi)):
            ax.plot([0], [0], [0], "o", color='tomato')
            ax.plot([self.l1], [0], [0], "o", color='tomato')
            theta = np.linspace(0, 2 * np.pi, 201)
            ax.plot(self.l4 * np.cos(theta) + self.l1, self.l4 * np.sin(theta), np.zeros(len(theta)), color='grey', linestyle='--', linewidth=0.5)
            ax.plot(A[0, :], A[1, :], A[2, :], color="lightseagreen", linewidth=0.5)
            ax.plot(M[0, :], M[1, :], M[2, :], color="lightseagreen", linewidth=0.5)
            ax.plot(B[0, :], B[1, :], B[2, :], color="lightseagreen", linewidth=0.5)
            ax.plot([A[0, i], M[0, i]], [A[1, i], M[1, i]], [A[2, i], M[2, i]], color="royalblue", linewidth=2)
            ax.plot([A[0, i], 0], [A[1, i], 0], [A[2, i], 0], "royalblue", linewidth=2)
            ax.plot([B[0, i], self.l1], [B[1, i], 0], [B[2, i], 0], "royalblue", linewidth=2)
            ax.plot([A[0, i]], [A[1, i]], [A[2, i]], "o", color='grey', linewidth=1)
            ax.plot([B[0, i]], [B[1, i]], [B[2, i]], "o", color='grey', linewidth=1)
            ax.plot([M[0, i]], [M[1, i]], [M[2, i]], "o", color='grey', linewidth=1)

            camera.snap()

        anim = camera.animate()
        anim.save("lambda animation 3d.gif", writer="imagemagick")
        plt.close(fig)

    def __calculate_bearing_animation(self, i, A, B, delta_B_normal, delta_B_angle, cylinder, Hb, Rb, Rh, ax):
        rotB = self.__rotation(delta_B_normal[:, i], delta_B_angle[i])
        cylinder_rotation = rotB.dot(cylinder)
        cylinder_rotation_out = [[], [], []]
        for j in range(len(cylinder_rotation[0])):
            if cylinder_rotation[0][j] ** 2 + cylinder_rotation[1][j] ** 2 > Rh ** 2:
                cylinder_rotation_out[0].append(cylinder_rotation[0][j])
                cylinder_rotation_out[1].append(cylinder_rotation[1][j])
                cylinder_rotation_out[2].append(cylinder_rotation[2][j])
        cylinder_rotation_out = np.array(cylinder_rotation_out)

        BC = np.array([self.l1, 0, 0]).T - B[:, i]
        BC = BC * (Rb / np.sqrt(np.linalg.norm(BC.T.dot(BC))))
        BA = A[:, i] - B[:, i]
        BA = Rb * BA / np.linalg.norm(BA)

        ax.plot(cylinder_rotation_out[0, :], cylinder_rotation_out[1, :], cylinder_rotation_out[2, :], ".", color='red')
        ax.plot([0, 0], [0, 0], [- Hb / 2, Hb / 2], 'b')
        ax.plot([0, BC[0]], [0, BC[1]], [- Hb / 2, - Hb / 2 + BC[2]], 'blue')
        ax.plot([0, BA[0]], [0, BA[1]], [- Hb / 2, - Hb / 2 + BA[2]], 'red')
        ax.plot([0, BC[0]], [0, BC[1]], [Hb / 2, Hb / 2 + BC[2]], 'blue')
        ax.plot([0, BA[0]], [0, BA[1]], [Hb / 2, Hb / 2 + BA[2]], 'red')
        ax.plot(Rb * np.cos(self.phi), Rb * np.sin(self.phi), np.zeros(len(self.phi)) - Hb / 2, 'blue')
        ax.plot(Rb * np.cos(self.phi), Rb * np.sin(self.phi), np.zeros(len(self.phi)) + Hb / 2, 'blue')

    def draw_bearing_animation_3d(self, A, M, B, delta_B_normal, delta_B_angle):
        fig = pylab.figure()
        ax = Axes3D(fig)
        ax.set_title("Bearing mechanism 3d animation")
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        camera = Camera(fig)

        Rb = 17.5
        Hb = 10
        Rbh = 0.03
        Rh = Rb + Rbh
        h = np.linspace(- Hb / 2, Hb / 2, 30)
        cylinder = np.array([Rb * np.kron(np.cos(self.phi), np.ones((1, len(h))))[0],
             Rb * np.kron(np.sin(self.phi), np.ones((1, len(h))))[0],
             np.kron(np.ones((1, len(self.phi))), h)[0]])

        for i in range(len(delta_B_angle)):
            self.__calculate_bearing_animation(i, A, B, delta_B_normal, delta_B_angle, cylinder, Hb, Rb, Rh, ax)
            camera.snap()

        anim = camera.animate()
        anim.save("Bearing mechanism animation 3d.gif", writer="imagemagick")
        plt.close(fig)

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
