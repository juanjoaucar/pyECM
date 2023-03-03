import numpy as np
import pytest

from pyECM.molecule_class import molecula
from pyECM.data.moleculas import mol_SeIClO_rotated
import matplotlib.pyplot as plt
import os

main_directory=os.path.realpath(os.path.dirname(__file__))+'/../'
# ===============================================================
# Plots
# ===============================================================
@pytest.mark.mpl_image_compare
def test_molecule_plot_1():
    options = {
            'xlabel' : "ax.set_xlabel('$x$', fontsize=10, rotation=0)",
            'ylabel' : "ax.set_ylabel('$y$', fontsize=10, rotation=0)",
            'zlabel' : "ax.set_zlabel('$z$', fontsize=10, rotation=0)"
            }
    fig=plt.figure()
    vector = np.array([-0.1219, -0.7434, 0.6576]) #Normal vector that defines the plane
    origen_achiral_atom = np.array([-1.31002, -0.33894, -0.62598]) # Point that defines the plane. Any atom of the symmetric structure
    c = molecula(figure = fig, XYZ_file = main_directory+'pyECM/data/import/AP1_chiral.xyz', dipolo=vector, origen=origen_achiral_atom,**options)
    c.rotate_to_align_dipole_with_z()
    c.plot_dipolo()
    c.plot_plano()
    c.plot_sphere()
    c.plot_opciones()

    return c.fig


@pytest.mark.mpl_image_compare
def test_molecule_plot_2():
    options = {
            'xlabel' : "ax.set_xlabel('$x$', fontsize=10, rotation=0)",
            'ylabel' : "ax.set_ylabel('$y$', fontsize=10, rotation=0)",
            'zlabel' : "ax.set_zlabel('$z$', fontsize=10, rotation=0)"
            }
    fig=plt.figure()
    vector = np.array([-0.1219, -0.7434, 0.6576]) #Normal vector that defines the plane
    origen_achiral_atom = np.array([-1.31002, -0.33894, -0.62598]) # Point that defines the plane. Any atom of the symmetric structure
    c = molecula(figure = fig, XYZ_file = main_directory+'pyECM/data/import/AP1_chiral.xyz', dipolo=vector, origen=origen_achiral_atom,**options)
    c.rotate_to_align_dipole_with_z()
    c.plot_plano()
    c.plot_dipolo()
    c.plot_sphere()
    c.plot_opciones()

    return c.fig


@pytest.mark.mpl_image_compare
def test_molecule_plot_3():
    options = {
            'xlabel' : "ax.set_xlabel('$x$', fontsize=10, rotation=0)",
            'ylabel' : "ax.set_ylabel('$y$', fontsize=10, rotation=0)",
            'zlabel' : "ax.set_zlabel('$z$', fontsize=10, rotation=0)"
            }
    fig=plt.figure()
    vector = np.array([-0.1219, -0.7434, 0.6576]) #Normal vector that defines the plane
    origen_achiral_atom = np.array([-1.31002, -0.33894, -0.62598]) # Point that defines the plane. Any atom of the symmetric structure
    c = molecula(figure = fig, XYZ_file = main_directory+'pyECM/data/import/AP1_chiral.xyz', dipolo=vector, origen=origen_achiral_atom,**options)
#    c.rotate_to_align_dipole_with_z()
    c.plot_opciones()
    c.plot_plano()
    c.plot_dipolo()
    c.plot_sphere()

    return c.fig


@pytest.mark.mpl_image_compare
def test_S_preloaded_molecule_plot():

    fig=plt.figure()
    options = {
            'xlabel' : "ax.set_xlabel('$x$', fontsize=10, rotation=0)",
            'ylabel' : "ax.set_ylabel('$y$', fontsize=10, rotation=0)",
            'zlabel' : "ax.set_zlabel('$z$', fontsize=10, rotation=0)",
            'titulo' : 'self.fig.suptitle("SeIClO test molecule", fontsize=10)'
            }
    s = molecula(figure=fig, preloaded_molecule=mol_SeIClO_rotated('S'),**options)
    s.plot_sphere()
    s.plot_enlaces()
    s.plot_opciones()

    return s.fig


@pytest.mark.mpl_image_compare
def test_R_preloaded_molecule_plot():

    fig=plt.figure()
    options = {
            'xlabel' : "ax.set_xlabel('$x$', fontsize=10, rotation=0)",
            'ylabel' : "ax.set_ylabel('$y$', fontsize=10, rotation=0)",
            'zlabel' : "ax.set_zlabel('$z$', fontsize=10, rotation=0)",
            'titulo' : 'self.fig.suptitle("SeIClO test molecule", fontsize=10)'
            }
    r = molecula(figure=fig, preloaded_molecule=mol_SeIClO_rotated('R'),**options)
    r.plot_enlaces()
    r.plot_sphere()
    r.plot_opciones()

    return r.fig
