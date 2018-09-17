"""
Math functions and constants for moldraw
"""

import numpy as np

def rotate_matrix(angle):
    """ Returns the 2D rotation matrix
        of angle """
    c = np.cos(angle)
    s = np.sin(angle)
    return np.array([[c, s],
                     [-s, c]])

""" Specific commonly used rotations """
R45 = rotate_matrix(np.pi/4)
R60 = rotate_matrix(np.pi/3)
R90 = rotate_matrix(np.pi/2)
R180 = rotate_matrix(np.pi)
R270 = rotate_matrix(3*np.pi/4)

def scale_matrix(N):
    """ Returns the 2x2 N-scale matrix """
    return N * np.identity(2)

def mirror_matrix(theta):
    """ Returns a matrix to flip around
        a line going through the origin
        with angle theta """
    c = np.cos(theta)
    s = np.sin(theta)
    return np.array([[c**2-s**2, 2*c*s],
                     [2*c*s, s**2-c**2]])

def distance(p1, p2):
    """ Distance between two points """
    return np.linalg.norm(p2-p1)

def angle(*args):
    """ Returns the angle between two vectors.
        If only one vector is given, the angle
        returned is between it and the x-axis,
        if more than two vectors are given,
        only the angle between the first two
        will be returned. """
    if len(args) < 1:
        return 0.0
    elif len(args) == 1:
        return np.arctan2(args[0][1], args[0][0])
    else:
        v1 = args[0].flatten()
        v2 = args[1].flatten()
        return np.arccos(np.dot(v1, v2) / (norm(v1) * norm(v2)))

def norm(vec):
    """ Calculates the norm of a vector """
    return np.linalg.norm(vec)

def normalized(vec):
    """ Returns the normalized form of a vector """
    l = norm(vec)
    if l != 0.0:
        return vec / l
    else:
        raise ArithmeticError('Zero vector can\'t be normalized!')

def set_length(vec, length):
    """ Returns a vector with the same
        angle as the imput vector, but
        any desired length """
    return normalized(vec) * length

def perp_2d_vec(vec):
    """ Returns a vector perpendicular to vec.
        Works only for 2D vectors! """
    return normalized(np.dot(R90, vec))
