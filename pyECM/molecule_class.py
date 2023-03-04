import os
import sys
import time


import mendeleev as mendeleev
import numpy as np
from numpy import matmul as mm
from numpy import transpose as tp
from pyscf import gto, lib, scf
from pyscf.lib.misc import light_speed
from scipy.linalg import fractional_matrix_power as matrix_power
from pyECM.geometric_figures import define_plane, define_sphere
from pyECM.geometric_figures import plot_vector, rotation_matrix_from_vectors

module_path = os.path.abspath(os.path.join(".."))

if module_path not in sys.path:
    sys.path.append(module_path)


class molecula:
    """It is possible to create a molecule object from a xyz file.
    A vector can be associated to the molecule so that it is possible to align it
    with the z direction, as the virtual mirror path will be created
    by reflections in the z=0 plane.
    The vector should be obtained uploading the xzy file to
    https://csm.ouproj.org.il/molecule.

    :examples: >>> vector = np.array([-0.1807, -0.9725, -0.1469])
        >>> origin_at = np.array([0.0000, 0.200, 0.1000])
        >>> file='pyECM/data/import/CFMAR_chiral.xyz'
        >>> molecule = molecula(XYZ_file = file, dipolo=vector, origen=origin_at)
    :param preloaded_molecule: _description_, defaults to None
    :type preloaded_molecule: _type_, optional
    :param figure: _description_, defaults to None
    :type figure: _type_, optional
    :param campo: _description_, defaults to None
    :type campo: _type_, optional
    :param origen: _description_, defaults to None
    :type origen: _type_, optional
    :param XYZ_file: _description_, defaults to None
    :type XYZ_file: _type_, optional
    :param dipolo: _description_, defaults to None
    :type dipolo: _type_, optional
    """

    def __init__(
        self,
        preloaded_molecule=None,
        figure=None,
        campo=None,
        origen=None,
        XYZ_file=None,
        dipolo=None,
        **kwargs
    ):
        self.fig = figure
        self.opciones = kwargs
        self.origen = origen
        self.bohrtoang = 0.529177249

        if preloaded_molecule is not None:
            self.preloaded_molecule = preloaded_molecule
            self.nro_atoms, self.atoms, self.dipolo, self.bonds = preloaded_molecule

        if XYZ_file is not None:
            self.XYZ_file = XYZ_file
            self.load_from_xyz(filename=self.XYZ_file)
            self.dipolo = dipolo

        # self.validar()
        self.atoms_positions()
        self.coordenadas_central()

        if self.origen is not None:
            self.origin_on_atom()
        else:
            self.origen = np.array(
                [self.positions[0][0], self.positions[0][1], self.positions[0][2]]
            )

    # def validar(self):
    #    if not (type(self.nro_atoms) == int and type(self.atoms) == tuple):
    #        raise TypeError("Type error on nro_atoms and/or atoms variable/s")

    def atoms_positions(self):
        """Defines nuclei positions."""
        positions = []
        for i in range(self.nro_atoms):
            positions.append([self.atoms[0][i], self.atoms[1][i], self.atoms[2][i]])
        self.positions = np.array(positions)

    def coordenadas_central(self):
        """Defines central point."""
        if self.origen is None:
            pass
        elif type(self.origen) == str:
            for i in range(self.nro_atoms):
                if self.atoms[4][i] == self.origen:
                    self.punto_central = np.array(
                        [
                            self.positions[i][0],
                            self.positions[i][1],
                            self.positions[i][2],
                        ]
                    )
        elif type(self.origen) == np.ndarray:
            self.punto_central = self.origen

    def origin_on_atom(self):
        """Gets central point coordinates."""
        if type(self.origen) == str:
            for i in range(self.nro_atoms):
                if self.atoms[4][i] == self.origen:
                    coordenadas_punto_central = np.array(
                        [
                            self.positions[i][0],
                            self.positions[i][1],
                            self.positions[i][2],
                        ]
                    )
        if type(self.origen) == np.ndarray:
            coordenadas_punto_central = self.origen
        new_positions = np.zeros(3)
        for i in range(self.nro_atoms):
            new_line = self.positions[i] - coordenadas_punto_central
            new_positions = np.vstack([new_positions, new_line])
        new_positions = np.delete(new_positions, (0), axis=0)  # Borro fila 1 de ceros
        self.positions = new_positions
        self.coordenadas_punto_central = coordenadas_punto_central

    def rotate_to_align_dipole_with_z(self):
        """Rotates the molecule so that its vector is aligned with z."""
        direction = np.array([0, 0, 1])
        if direction.all() != self.dipolo.all():
            matrix = rotation_matrix_from_vectors(self.dipolo, direction)
            for i in range(self.nro_atoms):
                self.positions[i] = matrix @ self.positions[i]
            self.dipolo = matrix @ self.dipolo

    def plot_dipolo(self):
        """Plots the molecule dipole."""
        point = np.array([0, 0, 0])
        dipolo = np.array([self.dipolo[0], self.dipolo[1], self.dipolo[2]])
        plot_vector(self.fig, point, dipolo)

    def plot_plano(self):
        """Plots the plane normal to the dipole."""
        point = np.array([0.0, 0.0, 0.0])
        dipolo = np.array([self.dipolo[0], self.dipolo[1], self.dipolo[2]])
        if self.fig.get_axes():
            ax = self.fig.gca()
        else:
            ax = self.fig.add_subplot(projection="3d")
        define_plane(ax, point, dipolo, size=1)

    def plot_sphere(self):
        """Plots the nuclei as spheres."""
        r = 0.05
        if self.fig.get_axes():
            ax = self.fig.gca()
        else:
            ax = self.fig.add_subplot(projection="3d")

        for i in range(self.nro_atoms):
            (xs, ys, zs) = define_sphere(
                self.positions[i][0], self.positions[i][1], self.positions[i][2], r
            )
            ax.plot_wireframe(xs, ys, zs, color=self.atoms[3][i])

    def plot_enlaces(self):
        """Plot the molecule bonds."""
        for i in range(len(self.bonds)):
            Ax = float(np.array(self.positions)[self.bonds[i, 0], 0])
            Ay = float(np.array(self.positions)[self.bonds[i, 0], 1])
            Az = float(np.array(self.positions)[self.bonds[i, 0], 2])
            Bx = float(np.array(self.positions)[self.bonds[i, 1], 0])
            By = float(np.array(self.positions)[self.bonds[i, 1], 1])
            Bz = float(np.array(self.positions)[self.bonds[i, 1], 2])
            # draw diagonal line from (Ax,Ay,Az) to (Bx,By,Bz)
            # ax = self.fig.gca(projection="3d")
            if self.fig.get_axes():
                ax = self.fig.gca()
            else:
                ax = self.fig.add_subplot(projection="3d")

            ax.plot([Ax, Bx], [Ay, By], zs=[Az, Bz])

    def plot_opciones(self):
        """Options plots."""
        if self.fig.get_axes():
            ax = self.fig.gca()
        else:
            ax = self.fig.add_subplot(projection="3d")
        ax.set_axis_on()
        for i in self.opciones:
            exec(self.opciones[i])

    def save_xyz(self, filename="MOL"):
        """Saves the molecule in xyz format

        :param filename: xyz file name, defaults to "MOL"
        :type filename: str, optional
        """
        atoms_name_xyz = []
        atoms_Z = []
        for j in range(self.nro_atoms):
            atom_name_xyz = "".join(
                [i for i in str(self.atoms[4][j]) if not i.isdigit()]
            )
            atoms_name_xyz.append(atom_name_xyz)
            atoms_Z.append(eval("mendeleev." + atom_name_xyz + ".atomic_number"))

        # filas =  np.zeros(4)
        fila_0 = [self.nro_atoms, "", "", ""]
        fila_1 = ["XYZ file", "", "", ""]
        filas = np.vstack([fila_0, fila_1])
        for j in range(self.nro_atoms):
            new_line = np.array(
                [
                    self.atoms[4][j],
                    self.positions[j][0],
                    self.positions[j][1],
                    self.positions[j][2],
                ]
            )
            filas = np.vstack([filas, new_line])
        # with open(filename, "ab") as f:
        #    np.savetxt(f, filas, fmt="%s")

        with open(filename, "w") as f:
            f.write(str(self.nro_atoms) + "\n")
            f.write("XYZ file created by pyECM\n")
            for j in range(self.nro_atoms):
                f.write(
                    atoms_name_xyz[j]
                    + "   "
                    + "{:.6f}".format(self.positions[j][0])
                    + " "
                    + "{:.6f}".format(self.positions[j][1])
                    + " "
                    + "{:.6f}".format(self.positions[j][2])
                    + "\n"
                )

    def load_from_xyz(self, filename="MOL"):
        """Loads the molecule from a xyz file.

        :param filename: xyz file name, defaults to "MOL"
        :type filename: str, optional
        """

        xyz = open(filename)
        n_atoms = int(xyz.readline())
        xyz.readline()  # title

        coord_x = np.zeros(n_atoms)
        coord_y = np.zeros(n_atoms)
        coord_z = np.zeros(n_atoms)
        colores = ["black"] * n_atoms
        nombres = [None] * n_atoms

        i = 0
        for line in xyz:
            atom, x, y, z = line.split()
            coord_x[i] = x
            coord_y[i] = y
            coord_z[i] = z
            nombres[i] = atom
            i += 1
        xyz.close()
        self.nro_atoms = n_atoms
        self.atoms = (coord_x, coord_y, coord_z, colores, nombres)

    def xyz_mirror_path(
        self,
        prefix_name="MOL",
        DIRAC=False,
        folder=None,
        lim_inf=0.5,
        lim_sup=1.0,
        points=11,
    ):
        """Generates the molecules in the virtual mirror path, which is defined
        by reflecting the molecule in the plane z=0.

        :param prefix_name: prefix name for xyz files, defaults to "MOL"
        :type prefix_name: str, optional
        :param DIRAC: create DIRAC mol files, defaults to False
        :type DIRAC: bool, optional
        :param folder: directory where saving the files, defaults to None
        :type folder: str, optional
        :param lim_inf: minimum z-rate, defaults to 0.0
        :type lim_inf: float, optional
        :param lim_sup: maximum z-rate, defaults to 1.0
        :type lim_sup: float, optional
        :param points: number of z-rates defining the path, defaults to 10.0
        :type points: float, optional
        """
        atoms_name_xyz = []
        atoms_name_dirac = []
        atoms_Z = []
        for j in range(self.nro_atoms):
            atom_name_xyz = "".join(
                [i for i in str(self.atoms[4][j]) if not i.isdigit()]
            )
            atom_name_dirac = str(self.atoms[4][j])

            atoms_name_xyz.append(atom_name_xyz)
            atoms_name_dirac.append(atom_name_dirac + str(j + 1))
            atoms_Z.append(eval("mendeleev." + atom_name_xyz + ".atomic_number"))

        if points is None:
            points = int((lim_sup - lim_inf) / 0.05 + 1)

        # We always need the nearest assymetric structure
        if lim_inf != 0:
            filename = folder + prefix_name + "_" + "0.00.xyz"
            fila_0 = [self.nro_atoms, "", "", ""]
            fila_1 = ["XYZ file", "", "", ""]
            filas = np.vstack([fila_0, fila_1])
            for j in range(self.nro_atoms):
                new_line = np.array(
                    [
                        atoms_name_xyz[j],
                        self.positions[j][0],
                        self.positions[j][1],
                        self.positions[j][2] * 0,
                    ]
                )
                filas = np.vstack([filas, new_line])
            # Remove xyz files if they exist
            try:
                os.remove(filename)
            except OSError:
                pass

            with open(filename, "ab") as f:
                np.savetxt(f, filas, fmt="%s")

        # z rate is i divided by 100
        for i in np.linspace(lim_inf, lim_sup, points):
            # filename = folder+prefix_name+'_'+str(round(i,2))+'.xyz'
            filename = folder + prefix_name + "_" + "{:.2f}".format(i) + ".xyz"
            fila_0 = [self.nro_atoms, "", "", ""]
            fila_1 = ["XYZ file", "", "", ""]
            filas = np.vstack([fila_0, fila_1])
            for j in range(self.nro_atoms):
                new_line = np.array(
                    [
                        atoms_name_xyz[j],
                        self.positions[j][0],
                        self.positions[j][1],
                        self.positions[j][2] * i,
                    ]
                )
                filas = np.vstack([filas, new_line])

            # Remove xyz files if they exist
            try:
                os.remove(filename)
            except OSError:
                pass

            with open(filename, "w") as f:
                f.write(str(self.nro_atoms) + "\n")
                f.write("XYZ file created by pyECM\n")
                for j in range(self.nro_atoms):
                    f.write(
                        atoms_name_xyz[j]
                        + "   "
                        + "{:.6f}".format(self.positions[j][0])
                        + " "
                        + "{:.6f}".format(self.positions[j][1])
                        + " "
                        + "{:.6f}".format(self.positions[j][2] * i)
                        + "\n"
                    )

            if DIRAC:
                filename_DIRAC = (
                    folder + "D22_" + prefix_name + "_" + "{:.2f}".format(i) + ".mol"
                )
                with open(filename_DIRAC, "w") as dirac:
                    dirac.write("DIRAC\n")
                    dirac.write("\n")
                    dirac.write("\n")
                    dirac.write(
                        "C " + "{:3d}".format(self.nro_atoms) + "  0 0         A\n"
                    )
                    for j in range(self.nro_atoms):
                        dirac.write("     " + "{:3d}".format(atoms_Z[j]) + ".     1\n")
                        dirac.write(
                            atoms_name_dirac[j]
                            + "   "
                            + "{:.6f}".format(self.positions[j][0])
                            + " "
                            + "{:.6f}".format(self.positions[j][1])
                            + " "
                            + "{:.6f}".format(self.positions[j][2] * i)
                            + "\n"
                        )
                        dirac.write("LARGE BASIS base\n")
                    dirac.write("FINISH")

    def CCM_on_path(self, lim_inf=0.20, lim_sup=1.00, points=10):
        """Calculates the CCM for the molecules in the virtual mirror path.

        :param lim_inf: minimum z-rate, defaults to 0.20
        :type lim_inf: float, optional
        :param lim_sup: maximum z-rate, defaults to 1.00
        :type lim_sup: float, optional
        :param points: number of z-rates where calculating the CCM, defaults to 10
        :type points: int, optional
        :return: z-rates, NORMs1, CCMs1, NORMs2, CCMs2
        :rtype: numpy.ndarray(s)
        """
        print(
            "Warning: it is assumed that the reflection plane (and symm conf) is at z=0"
        )
        z_rate = np.zeros(points)
        CCMs_1 = np.zeros(points)
        Norms_1 = np.zeros(points)
        CCMs_2 = np.zeros(points)
        Norms_2 = np.zeros(points)
        index = 0
        for j in np.linspace(lim_inf, lim_sup, points):
            z_rate[index] = j
            coordenadas_x = self.atoms[0]
            coordenadas_y = self.atoms[1]
            coordenadas_z = self.atoms[2] * j

            mean_x = np.average(coordenadas_x)
            mean_y = np.average(coordenadas_y)
            mean_z = np.average(coordenadas_z)

            # Method 1
            # Add citation

            mendeleev.Fe.atomic_weight
            atomic_weights = np.zeros(self.nro_atoms)
            i = 0
            for name in self.atoms[4]:
                atomic_weights[i] = eval("mendeleev." + name + ".atomic_weight")
                i = i + 1

            total_mass = np.sum(atomic_weights)
            center_of_mass_x = np.sum(coordenadas_x * atomic_weights) / total_mass
            center_of_mass_y = np.sum(coordenadas_y * atomic_weights) / total_mass
            center_of_mass_z = np.sum(coordenadas_z * atomic_weights) / total_mass
            distances_to_CM_x = coordenadas_x - center_of_mass_x
            distances_to_CM_y = coordenadas_y - center_of_mass_y
            distances_to_CM_z = coordenadas_z - center_of_mass_z
            distances_to_CM = np.sqrt(
                distances_to_CM_x**2 + distances_to_CM_y**2 + distances_to_CM_z**2
            )

            D = np.max(distances_to_CM)
            Norm_1 = D**2 * self.nro_atoms
            CCM_1 = (1 / Norm_1) * np.sum(coordenadas_z**2) * 100

            # Method 2
            # Add citation
            Norm_2 = np.sum(
                (coordenadas_x - mean_x) ** 2
                + (coordenadas_y - mean_y) ** 2
                + (coordenadas_z - mean_z) ** 2
            )
            suma_2 = np.sum(coordenadas_z**2)
            CCM_2 = suma_2 / Norm_2 * 100

            CCMs_1[index] = CCM_1
            Norms_1[index] = Norm_1
            CCMs_2[index] = CCM_2
            Norms_2[index] = Norm_2
            index = index + 1

        return z_rate, Norms_1, CCMs_1, Norms_2, CCMs_2

    def ECM_on_path(
        self,
        debug=0,
        tracking=False,
        basis_set="sto-6g",
        name=None,
        lim_inf=0.20,
        lim_sup=1.00,
        points=10,
        NR=True,
        fourcomp=False,
        cvalue=137.03599967994,
        cartesian=False,
    ):
        """Calculate ECM in the virtual mirror path

        :param debug: Debug printing level, defaults to 0
        :type debug: int, optional
        :param tracking: Print results for each z-rate
        :type tracking: bool, optional
        :param basis_set: basis set for describing
            the electronic wave function, defaults to "sto-6g"
        :type basis_set: str, optional
        :param name: name of the xyz molecule file,
            including its directory, defaults to None
        :type name: str, optional
        :param lim_inf: minimum z-rate, defaults to 0.20
        :type lim_inf: float, optional
        :param lim_sup: maximum z-rate, defaults to 1.00
        :type lim_sup: float, optional
        :param points: number of z-rates where calculating the ECM, defaults to 10
        :type points: int, optional
        :param NR: calculates ECM at NR level, defaults to True
        :type NR: bool, optional
        :param fourcomp: calculates ECM at 4c level, defaults to False
        :type fourcomp: bool, optional
        :param cvalue: speed-ligh velocity, defaults to 137.03599967994
        :type cvalue: float, optional
        :param cartesian: use cartesian basis set, defaults to False
        :type cartesian: bool, optional
        :return: z-rates, ECMs(NR), molecular orbital contributions
            to ECMs(NR), ECMs(4c)
        :rtype: numpy.ndarray(s)
        """

        z_rate = np.zeros(points)
        ECMs_NR = np.zeros(points)
        ovlp_NR = np.zeros(points)
        N_ovlp_NR = np.zeros(points)  # normalized
        ECMs_4c = np.zeros(points)
        ovlp_4c = np.zeros(points)
        N_ovlp_4c = np.zeros(points)  # normalized
        ECMs_molcontr = []
        index = 0
        if tracking:
            print(name)

        for j in np.linspace(lim_inf, lim_sup, points):
            z_rate[index] = "{:.2f}".format(j)
            mol_chiral = gto.M(
                atom=name + "_" + "{:.2f}".format(j) + ".xyz",
                basis=basis_set,
                verbose=0,
                max_memory=5000.0,
            )

            mol_achiral = gto.M(
                atom=name + "_0.00.xyz", basis=basis_set, verbose=0, max_memory=5000.0
            )

            if cartesian:
                mol_chiral, ctr_coeff1 = mol_chiral.to_uncontracted_cartesian_basis()
                mol_achiral, ctr_coeff2 = mol_achiral.to_uncontracted_cartesian_basis()

            mol_super = mol_chiral + mol_achiral

            if NR:
                start_NRtime = time.time()
                mf_chiral = scf.RHF(mol_chiral)
                mf_chiral.kernel()

                # from pyscf import dft
                # mf_chiral = dft.RKS(mol_chiral)
                # mf_chiral.xc = 'b3lyp'
                # mf_chiral.kernel()

                naos_sph = mol_chiral.intor("int1e_ovlp_sph").shape[0]
                occupied_MO = mf_chiral.mol.nelec[0]
                AO_number = mol_chiral.nao
                AO_number_supermol = np.array([mol_super.nao])[0]

                overlap_mixed_fullspace = mol_super.intor("int1e_ovlp")
                overlap_mixed = overlap_mixed_fullspace[
                    0 : int(AO_number_supermol / 2),
                    int(AO_number_supermol / 2) : AO_number_supermol,
                ]

                overlap_chiral = mol_chiral.intor("int1e_ovlp")
                overlap_achiral = mol_achiral.intor("int1e_ovlp")

                overlap_pot_chiral = matrix_power(overlap_chiral, 0.5)
                overlap_pot_achiral = matrix_power(overlap_achiral, -0.5)

                ocupp_mo_coeff = mf_chiral.mo_coeff[0:AO_number, 0:occupied_MO]
                norma_chiral = np.trace(
                    mm(mm(tp(ocupp_mo_coeff), overlap_chiral), ocupp_mo_coeff)
                )

                # old basis: chiral basis set
                # new basis: achiral basis set
                C_achiral_newbasis = mm(
                    mm(overlap_pot_achiral, overlap_pot_chiral), mf_chiral.mo_coeff
                )

                achiral_norm = 0
                solapamiento_NR = 0
                for k in range(occupied_MO):
                    achiral_norm = (
                        achiral_norm
                        + mm(
                            mm(tp(C_achiral_newbasis).conjugate(), overlap_achiral),
                            C_achiral_newbasis,
                        )[k, k]
                    )
                    solapamiento_NR = (
                        solapamiento_NR
                        + mm(
                            mm(tp(mf_chiral.mo_coeff).conjugate(), overlap_mixed),
                            C_achiral_newbasis,
                        )[k, k]
                    )

                    ECMs_molcontr.append(
                        100
                        * (
                            1
                            - np.abs(
                                mm(
                                    mm(tp(mf_chiral.mo_coeff), overlap_mixed),
                                    C_achiral_newbasis,
                                )[k, k]
                            )
                        )
                    )
                #     print("Overlap contribution from",k,"orbital:"
                # ,mm(mm(tp(mf_chiral.mo_coeff),overlap_mixed),C_achiral_newbasis)[k,k])
                #     print(" ")

                ECMs_NR[index] = 100 * (1 - np.abs(solapamiento_NR) / norma_chiral)
                ovlp_NR[index] = solapamiento_NR
                N_ovlp_NR[index] = solapamiento_NR / occupied_MO

                end_NRtime = time.time()
                if debug > 0:
                    print("naos_cart:", mol_chiral.nao)
                    print("naos_sph:", naos_sph)
                    print("MO (doubly) occupied:", occupied_MO)
                    print("electrones:", 2 * occupied_MO)
                    print("norma WF chiral (chiral basis):", norma_chiral)
                    print("norma WF achiral (achiral basis):", achiral_norm)
                    print("overlap:", solapamiento_NR / norma_chiral)
                    print("NR time (min):", (end_NRtime - start_NRtime) / 60)

            # nuevo_NR = np.transpose(np.reshape(
            # np.ravel(ECMs_molcontr),(points,occupied_MO)))
            # return z_rate, ECMs, nuevo_NR

            if fourcomp:
                start_4ctime = time.time()
                # nuevo_4c = None
                with light_speed(cvalue):
                    c = lib.param.LIGHT_SPEED

                    mf_chiral_rel = scf.DHF(mol_chiral)
                    mf_chiral_rel.kernel()

                    n4c, nmo = mf_chiral_rel.mo_coeff.shape
                    n2c = n4c // 2
                    nNeg = nmo // 2  # molecular orbitals of negative energy
                    nocc = mol_chiral.nelectron
                    # nvir = nmo // 2 - nocc  # virtual orbitals
                    mo_pos_l = mf_chiral_rel.mo_coeff[:n2c, nNeg:]
                    mo_pos_s = mf_chiral_rel.mo_coeff[n2c:, nNeg:]  # * (.5/c)
                    Lo = mo_pos_l[:, :nocc]
                    So = mo_pos_s[:, :nocc]
                    # Lv = mo_pos_l[:,nocc:]
                    # Sv = mo_pos_s[:,nocc:]

                    from pyscf.scf.dhf import get_ovlp

                    overlap_chiral_4c = get_ovlp(mol_chiral)
                    overlap_achiral_4c = get_ovlp(mol_achiral)

                    overlap_chiral_large = overlap_chiral_4c[:n2c, :n2c]
                    overlap_chiral_small = overlap_chiral_4c[n2c:, n2c:]

                    overlap_achiral_large = overlap_achiral_4c[:n2c, :n2c]
                    overlap_achiral_small = overlap_achiral_4c[n2c:, n2c:]

                    mol_super = mol_chiral + mol_achiral
                    overlap_supermol_4c = get_ovlp(mol_super)

                    # Its ok:
                    # overlap_chiral_large -
                    # overlap_supermol_4c[:n2c,0:n2c]
                    # overlap_achiral_large -
                    # overlap_supermol_4c[n2c:2*n2c,n2c:2*n2c]
                    # overlap_chiral_small -
                    # overlap_supermol_4c[2*n2c:3*n2c,2*n2c:3*n2c]
                    # overlap_achiral_small -
                    # overlap_supermol_4c[3*n2c:4*n2c,3*n2c:4*n2c]

                    overlap_mixed_SchiralSachiral = overlap_supermol_4c[
                        2 * n2c : 3 * n2c, 2 * n2c : 3 * n2c
                    ]
                    overlap_mixed_LchiralLachiral = overlap_supermol_4c[
                        :n2c, n2c : 2 * n2c
                    ]

                    AO_number = mol_chiral.nao
                    AO_number_supermol = np.array([mol_super.nao])[0]

                    overlap_mixed_fullspace = mol_super.intor("int1e_ovlp_spinor")
                    overlap_mixed = overlap_mixed_fullspace[
                        0 : int(AO_number_supermol),
                        int(AO_number_supermol) : 2 * AO_number_supermol,
                    ]

                    overlap_chiral = mol_chiral.intor("int1e_ovlp_spinor")
                    overlap_achiral = mol_achiral.intor("int1e_ovlp_spinor")

                    overlap_pot_chiral = matrix_power(overlap_chiral, 0.5)
                    overlap_pot_achiral = matrix_power(overlap_achiral, -0.5)

                    # overlap_pot_chiral_large = matrix_power(overlap_chiral_large, 0.5)
                    overlap_pot_achiral_large = matrix_power(
                        overlap_achiral_large, -0.5
                    )

                    # overlap_pot_chiral_small = matrix_power(overlap_chiral_small, 0.5)
                    overlap_pot_achiral_small = matrix_power(
                        overlap_achiral_small, -0.5
                    )

                    overlap_ll_pot_chiral = matrix_power(overlap_chiral_large, 0.5)
                    # overlap_ll_pot_achiral=matrix_power(overlap_achiral_large,0.5)
                    overlap_ss_pot_chiral = matrix_power(overlap_chiral_small, 0.5)
                    # overlap_ss_pot_achiral=matrix_power(overlap_achiral_small,0.5)

                    norma_LoLo_chiral = np.trace(
                        mm(mm(tp(Lo).conjugate(), overlap_chiral_large), Lo)
                    )
                    norma_SoSo_chiral = np.trace(
                        mm(mm(tp(So).conjugate(), overlap_chiral_small), So)
                    )
                    norma_total_chiral = norma_LoLo_chiral + norma_SoSo_chiral

                    # MCOEFF new basis:
                    # C_Lo_achiral_newbasis =
                    # mm(mm(overlap_pot_achiral,overlap_ll_pot_chiral),Lo)
                    # C_So_achiral_newbasis =
                    # mm(mm(overlap_pot_achiral,overlap_ss_pot_chiral),So)
                    C_Lo_achiral_newbasis = mm(
                        mm(overlap_pot_achiral_large, overlap_ll_pot_chiral), Lo
                    )
                    C_So_achiral_newbasis = mm(
                        mm(overlap_pot_achiral_small, overlap_ss_pot_chiral), So
                    )

                    achiral_norm_So = 0
                    achiral_norm_Lo = 0
                    solapamiento_LoLo = 0
                    solapamiento_SoSo = 0
                    # test2=0
                    for k in range(nocc):
                        achiral_norm_Lo = (
                            achiral_norm_Lo
                            + mm(
                                mm(
                                    tp(C_Lo_achiral_newbasis).conjugate(),
                                    overlap_achiral_large,
                                ),
                                C_Lo_achiral_newbasis,
                            )[k, k]
                        )
                        achiral_norm_So = (
                            achiral_norm_So
                            + mm(
                                mm(
                                    tp(C_So_achiral_newbasis).conjugate(),
                                    overlap_achiral_small,
                                ),
                                C_So_achiral_newbasis,
                            )[k, k]
                        )
                        solapamiento_LoLo = (
                            solapamiento_LoLo
                            + mm(
                                mm(tp(Lo).conjugate(), overlap_mixed_LchiralLachiral),
                                C_Lo_achiral_newbasis,
                            )[k, k]
                        )
                        solapamiento_SoSo = (
                            solapamiento_SoSo
                            + mm(
                                mm(tp(So).conjugate(), overlap_mixed_SchiralSachiral),
                                C_So_achiral_newbasis,
                            )[k, k]
                        )

                    # test2 = test1 + mm(mm(tp(Lo).conjugate()
                    # ,overlap_chiral_large),Lo)[k,k]

                    solapamiento_total = (solapamiento_LoLo + solapamiento_SoSo).real
                    ECMs_4c[index] = 100 * (
                        1 - np.abs(solapamiento_total) / np.abs(norma_total_chiral)
                    )
                    ovlp_4c[index] = solapamiento_total.real
                    N_ovlp_4c[index] = solapamiento_total.real / nocc

                    end_4ctime = time.time()
                    if debug > 0:
                        print("LoLo Norm:", norma_LoLo_chiral)
                        print("SoSo Norm:", norma_SoSo_chiral)
                        print("Total (chiral) Norm:", norma_total_chiral)
                        print("Achiral LoLo Norm:", achiral_norm_Lo)
                        print("Achiral SoSo Norm:", achiral_norm_So)
                        print("LoLo chiral/achiral overlap:", solapamiento_LoLo / nocc)
                        print("SoSo chiral/achiral overlap:", solapamiento_SoSo / nocc)
                        print("suma solapamientos:", solapamiento_total / nocc)
                        print("cvalue", c)
                        print("4c energy", mf_chiral_rel.e_tot)
                        print("ECM LL+SS:", ECMs_4c[index])
                        print("4c time (min):", (end_4ctime - start_4ctime) / 60)

            if tracking:
                print(
                    "{:.2f}".format(z_rate[index]),
                    "{:.2f}".format(ECMs_NR[index]),
                    "{:.2f}".format(ECMs_4c[index]),
                )
            index = index + 1

        nuevo_NR = np.transpose(
            np.reshape(np.ravel(ECMs_molcontr), (points, occupied_MO))
        )
        return z_rate, ECMs_NR, nuevo_NR, ECMs_4c
