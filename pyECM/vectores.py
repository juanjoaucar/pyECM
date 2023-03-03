import numpy as np


def borrar_nulos(x, y, z=None):
    """Removes zeros values from input arrays

    :param x: First array to be treated
    :type x: array
    :param y: Second array to be treated
    :type y: array
    :param z: Third array to be treated, defaults to None
    :type z: array, optional
    :return: x, y, z without theirs zero values
    :rtype: arrays
    """
    index = 0
    indices_borrar = np.array([])

    if z is None:
        for i in x:
            if (i == 0) and (y[index] == 0):
                indices_borrar = np.append(indices_borrar, index)
            index = index + 1
        indices_borrar = indices_borrar.astype(int)
        resultado = (
            np.delete(x, indices_borrar),
            np.delete(y, indices_borrar),
            indices_borrar,
        )
    else:
        for i in x:
            if (i == 0) and (y[index] == 0) and (z[index] == 0):
                indices_borrar = np.append(indices_borrar, index)
            index = index + 1
        indices_borrar = indices_borrar.astype(int)
        resultado = (
            np.delete(x, indices_borrar),
            np.delete(y, indices_borrar),
            np.delete(z, indices_borrar),
            indices_borrar,
        )

    return resultado


def vector_to_versor(x, y, z=None):
    """Normalize 2D/3D vector (to unity)

    :param x: x component/s
    :type x: float or array
    :param y: y component/s
    :type y: float or array
    :param z: z component/s, defaults to None
    :type z: float or array, optional
    :return: Normalized vector/s
    :rtype: float or array
    """
    # z = np.array(z)
    if np.array(z).all() is None:
        n = np.sqrt(np.power(x, 2) + np.power(y, 2))
        versor = [x / n, y / n]
    else:
        n = np.sqrt(np.power(x, 2) + np.power(y, 2) + np.power(z, 2))
        versor = [x / n, y / n, z / n]
    return versor
