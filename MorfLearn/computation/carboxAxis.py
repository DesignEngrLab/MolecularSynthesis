import numpy as np

def calcMassCenter(masses, positions):
    massCenter = [sum(x) for x in
                  np.matrix.transpose(np.array([i * j for i, j in zip(np.array(masses), np.array(positions))]))] / sum(
        np.array(masses))
    return massCenter


def calcNewCoor(positions):
    nodes1 = [positions[0], positions[2]]
    nodes2 = positions[-2:]
    centers = np.array([sum(np.array(nodes1)) / 2, sum(np.array(nodes2)) / 2])
    x = centers[0] - centers[1]
    x = x / np.dot(x, x) ** 0.5
    nodes2 = np.array(nodes2)
    nodes1 = np.array(nodes1)
    directions = np.array([nodes1[0] - nodes1[1], nodes2[0] - nodes2[1]])
    if np.dot(directions[0], directions[1]) >= 0:
        y = directions[0] + directions[1]
    else:
        y = directions[0] - directions[1]
    z = np.cross(x, y)
    z = z / np.dot(z, z) ** 0.5
    y = np.cross(z, x)
    transMat = np.mat(np.matrix([x, y, z]))
    return transMat


def newPosition(position, coor, origin):
    transPosition = np.array((coor * np.mat(position - origin).T).T).flatten()
    return transPosition