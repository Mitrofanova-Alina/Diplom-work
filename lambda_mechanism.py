from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pylab
import matplotlib.pyplot as plt
import matplotlib.patches as ptc
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
        A = self.__rotation(self.delta_axis_normal, self.delta_axis_angle, self.phi)
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
        delta_B_normal = [- BA[1, :] / n, BA[0, :] / n, np.zeros(len(self.phi))]
        delta_B_angle = np.arctan(BA[2, :] / n)
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

    def __rotation(self, n, alpha, phi):
        c = np.cos(alpha)
        s = np.sin(alpha)

        rot = np.array([[
            c + (n[0] ** 2) * (1 - c), n[0] * n[1] * (1 - c) + n[2] * s, n[0] * n[2] * (1 - c) - n[1] * s
        ], [
            n[0] * n[1] * (1 - c) - n[2] * s, c + (n[1] ** 2) * (1 - c), n[0] * n[2] * (1 - c) + n[0] * s
        ], [
            n[0] * n[2] * (1 - c) + n[1] * s, n[1] * n[2] * (1 - c) - n[0] * s, c + (n[2] ** 2) * (1 - c)
        ]])

        A = (rot * self.l2).T.dot(np.array([np.cos(phi), np.sin(phi), np.zeros(len(phi))]))
        return A
