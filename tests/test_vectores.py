import numpy as np
import pytest

from pyECM.vectores import vector_to_versor
from pyECM.vectores import borrar_nulos


# ===============================================================
# vector_to_versor function
# ===============================================================
testdata1 = [(1.232,-2.500,1.544,[0.3866696068624928, -0.7846380009385, 0.4845924293796176]),]

@pytest.mark.parametrize("x, y, z, expected",testdata1, ids=None)


def test_vectortoversor(x,y,z, expected):
    """Testing vector conversion

    Examples
    --------
    >>>Atom('F', 0, 0, 1.97)
    {'element': 'F', 'x': 0, 'y': 0, 'z': 1.97}
    """
    input = str(vector_to_versor(x,y,z))

    assert input == str(expected)

@pytest.mark.parametrize(
    "x, y, expected_result",
    [
        (3,-6,[0.4472135954999579, -0.8944271909999159]),
    ],
)

# ===============================================================
# vector_to_versor (without z) function
# ===============================================================
def test_vectortoversor_withouz_z(x,y, expected_result):
    """Testing vector conversion

    Examples
    --------
    >>>Atom('F', 0, 0, 1.97)
    {'element': 'F', 'x': 0, 'y': 0, 'z': 1.97}
    """
    input = str(vector_to_versor(x,y))

    assert input == str(expected_result)

# ===============================================================
# borrar_nulos function
# ===============================================================
@pytest.mark.parametrize(
    "x, y, expected_result",
    [
        ([5,0,9,4],[5,0,-3,0],"(array([5, 9, 4]), array([ 5, -3,  0]), array([1]))"),
    ],
)
def test_borrarnulos(x,y,expected_result):
    """Test of removing null elements

    Examples
    --------
    >>>Atom('F', 0, 0, 1.97)
    {'element': 'F', 'x': 0, 'y': 0, 'z': 1.97}
    """
    input = str(borrar_nulos(x,y))

    assert input == str(expected_result)


@pytest.mark.parametrize(
    "x, y, z, expected_result",
    [
        ([0,0,1,-56],[0,0,-3,0],[0,53,-6,2],"(array([  0,   1, -56]), array([ 0, -3,  0]), array([53, -6,  2]), array([0]))"),
    ],
)
def test_borrarnulos_withz(x,y,z,expected_result):
    """Test of removing null elements

    Examples
    --------
    >>>Atom('F', 0, 0, 1.97)
    {'element': 'F', 'x': 0, 'y': 0, 'z': 1.97}
    """
    input = str(borrar_nulos(x,y,z))

    assert input == str(expected_result)
