import numpy as np

from lambda_mechanism import LambdaMechanism

if __name__ == "__main__":
    print("Hello world!")

    size = 1
    delta_axis_normal = [1, 0, 0]
    delta_axis_angle = 20 / 60
    phi0 = np.pi / 2
    phi = np.linspace(phi0, phi0 + 4 * np.pi, 100)
    lambda_mechanism = LambdaMechanism(size, phi, delta_axis_normal, delta_axis_angle)
    A, B, M, delta_B_normal, delta_B_angle = lambda_mechanism.motion()
    lambda_mechanism.draw_lambda_mechanism_motion_2d(A, M, B)
    lambda_mechanism.draw_lambda_mechanism_motion_3d(A, M, B)

    print("Bye, world!")
