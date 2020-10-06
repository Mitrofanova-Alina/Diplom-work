import numpy as np
from mpmath import *

from lambda_mechanism import LambdaMechanism
from interval_analysis import calculate_area_ellipse_circle_intersect, shoelace_formula

if __name__ == "__main__":
    print("Hello world!")

    size = 1
    delta_axis_normal = [1, 0, 0]
    delta_axis_angle = 1 / 30
    phi0 = np.pi / 2
    phi = np.linspace(phi0, phi0 + 2 * np.pi, 100)
    Rb = 17.5
    Hb = 10
    Rbh = 0.03
    lambda_mechanism = LambdaMechanism(size, phi, delta_axis_normal, delta_axis_angle, Rb, Hb, Rbh)
    lambda_mechanism.calculate_lambda_motion()
    # lambda_mechanism.draw_lambda_motion_2d()
    # lambda_mechanism.draw_lambda_motion_3d()
    lambda_mechanism.calculate_bearing_motion()
    # lambda_mechanism.draw_bearing_animation_3d()
    lambda_mechanism.draw_number_of_jamming_points()
    # lambda_mechanism.do_something(0, Hb / 2)
    # area = lambda_mechanism.calculate_area_of_intersection()
    # print(area)
    # lambda_mechanism.draw_area(area)

    # R1 = 5
    # a, b = 6, 4
    # # radius of the circle to be transformed
    # circle = np.array([R1 * np.cos(phi), R1 * np.sin(phi)])
    # ellipse = np.array([a * np.cos(phi), b * np.sin(phi)])
    # # rot = np.array([[0.5, -0.7], [-0.8, 0.6]])
    # # ellipse = rot.dot(ellipse)
    # calculate_area_ellipse_circle_intersect(R1, a, b, 0.01)

    # points = [[0, 0], [4, 0], [4, 4], [0, 4], [0, 0]]
    # eps = 0.5
    # for i in range(len(points)):
    #     points[i] = iv.matrix([[mpi(points[i][0] - eps, points[i][0] + eps)], [mpi(points[i][1] - eps, points[i][1] + eps)]])
    # print(shoelace_formula(points))
    print("Bye, world!")
