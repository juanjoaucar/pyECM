from pyECM.geometric_figures import plot_vector
import matplotlib.pyplot as plt
import pytest

@pytest.mark.mpl_image_compare
def test_plot_vector_1():
    origin = [0,0,0]
    end = [0.04,0.05,-0.03]
    color = "blue"
    fig=plt.figure()

    return plot_vector(fig,origin,end,color)

@pytest.mark.mpl_image_compare
def test_plot_vector_2():
    origin = [0,0,0]
    end = [0.04,0.05,-0.03]
    color = "blue"
    fig=plt.figure()
    fig.add_subplot(projection="3d")

    return plot_vector(fig,origin,end,color)
