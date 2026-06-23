import numpy as np

def normalize_vector(v):
    return v / np.linalg.norm(v)

def line_vector(p1, p2):
    return np.asarray(p2) - np.asarray(p1)

def plane_normal(p1, p2, p3):
    return np.cross(line_vector(p1, p2), line_vector(p1, p3))

def distance(p1, p2):
    return np.linalg.norm(line_vector(p1, p2))

def angle(p1, p2, p3):
    v1 = line_vector(p2, p1)
    v2 = line_vector(p2, p3)
    cos_theta = np.dot(v1, v2)
    cos_theta /= np.linalg.norm(v1) * np.linalg.norm(v2)
    return np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))

def dihedral(p1, p2, p3, p4):
    b0 = line_vector(p2, p1)
    b1 = line_vector(p2, p3)
    b2 = line_vector(p3, p4)
    b1 = normalize_vector(b1)
    v = np.cross(b0, b1)
    w = np.cross(b2, b1)

    return np.degrees(np.arctan2(np.dot(np.cross(v, w), b1), np.dot(v, w)))