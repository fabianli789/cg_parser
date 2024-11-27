from typing import (
    TYPE_CHECKING,
)

if TYPE_CHECKING:
    from nomad.datamodel.datamodel import (
        EntryArchive,
    )
    from structlog.stdlib import (
        BoundLogger,
    )

import numpy as np
from nomad.config import config
from nomad.datamodel.data import Schema
from nomad.metainfo import Quantity, SchemaPackage, Section, MSection, SubSection
from runschema.run import Run
from runschema.calculation import Calculation

configuration = config.get_plugin_entry_point(
    'cg_parser.schema_packages:schema_package_entry_point'
)

m_package = SchemaPackage()

class CGIntegration(MSection):
    m_def = Section(validate=False)
    name=Quantity(type=str, description='Describes what is being integrated')
    parameters = Quantity(type=str, shape=['*'], description='Parameters necessary for integration.')
class CGBox_size(MSection):
    m_def = Section(validate=False)
    xlo = Quantity(type=float, description='lower boundery for x-direction of simulation box, in Angstrom.')
    xhi = Quantity(type=float, description='upper boundery for x-direction of simulation box, in Angstrom.')
    ylo = Quantity(type=float, description='lower boundery for y-direction of simulation box, in Angstrom.')
    yhi = Quantity(type=float, description='upper boundery for y-direction of simulation box, in Angstrom.')
    zlo = Quantity(type=float, description='lower boundery for z-direction of simulation box, in Angstrom.')
    zhi = Quantity(type=float, description='upper boundery for z-direction of simulation box, in Angstrom.')

class CGMolecule_velocities(MSection):
    m_def = Section(validate=False)
    v_x = Quantity(type=np.float64, shape=['*'], description='x-velocity of molecules in Angstrom/femtosecond.')
    v_y = Quantity(type=np.float64, shape=['*'], description='y-velocity of molecules in Angstrom/femtosecond.')
    v_z = Quantity(type=np.float64, shape=['*'], description='z-velocity of molecules in Angstrom/femtosecond.')
    w_x = Quantity(type=np.float64, shape=['*'], description='Angular x-velocity of molecules in radians/femtosecond.')
    w_y = Quantity(type=np.float64, shape=['*'], description='Angular y-velocity of molecules in radians/femtosecond.')
    w_z = Quantity(type=np.float64, shape=['*'], description='Angular z-velocity of molecules in radians/femtosecond.')
class CGEllipsoid_properties(MSection):
    m_def = Section(validate=False)
    length = Quantity(type=np.float64, shape=['*'], description='Length of ellipsoids in Angstrom. Length of this array must fit n_ellipsoids.')
    width = Quantity(type=np.float64, shape=['*'], description='Width of ellipsoids in Angstrom.')
    height = Quantity(type=np.float64, shape=['*'], description='Height of ellipsoids in Angstrom.')
    w_orientation = Quantity(type=np.float64, shape=['*'], description='W-orientation of ellipsoids.')
    x_orientation = Quantity(type=np.float64, shape=['*'], description='X-orientation of ellipsoids.')
    y_orientation = Quantity(type=np.float64, shape=['*'], description='Y-orientation of ellipsoids.')
    z_orientation = Quantity(type=np.float64, shape=['*'], description='Z-orientation of ellipsoids.')
class CGInput(MSection):
    m_def  = Section(validate=False)
    n_atoms = Quantity(type=int, description='Number of atoms')
    n_ellipsoids = Quantity(type=int, description='Number of ellipsoids')
    n_bonds = Quantity(type=int, description='Number of bonds')
    n_angles = Quantity(type=int, description='Number of angles')
    n_dihedrals = Quantity(type=int, description='Number of dihedrals')
    n_impropers = Quantity(type=int, description='Number of impropers')
    n_atom_types = Quantity(type=int, description='Number of atom types')
    box_size  = SubSection(sub_section=CGBox_size.m_def)
    atom_type = Quantity(type=np.int32, shape=['*'], description='vector containing atom type')
    ellipsoid_flag = Quantity(type=np.int32, shape=['*'], description = '1 == ellipsoid, 0 == spherical')
    density_of_particle = Quantity(type=np.float64, shape=['*'], description='density of particle in g/cm3')
    ellipsoid_properties = SubSection(sub_section=CGEllipsoid_properties.m_def)
    molecule_velocities = SubSection(sub_section=CGMolecule_velocities.m_def)
    atom_coordinates = Quantity(type=np.float64, shape=['*', 3], description='Cartesian coordinates of atoms in Angstrom.')
    
class CGVacf(MSection):
    m_def = Section(validate=False)
    
    vacf = Quantity(type=np.float64, shape=['*'], description='velocity auto-correlation function')
    vacf_time = Quantity(type=np.float64, shape=['*'], description='time steps in femtoseconds corresponding to vacf-array.')
class CGSettings(MSection):
    m_def = Section(validate=False)
    ndim = Quantity(type=int, description='number of dimensions of simulation')
    seed = Quantity(type=float, description='seed for random number generator')
    rc_global = Quantity(type=float, description='global cutoff radius in Angstrom')
    skin = Quantity(type=str, description='surround skin for neighbor list in Angstrom')
    dw_c = Quantity(type=float, description='softness of conservative force')
    l_c = Quantity(type=float, description='length of particle in conservative force in Angstrom')
    d_c = Quantity(type=float, description='width of particle in conservative force in Angstrom')
    esp0_c = Quantity(type=float, description='well-depth of cross configuration in conservative force in kcal/mol')
    esp_rat_c = Quantity(type=float, description='well-depth ratio of face-to-face and side-to-side in conservative force')
    dw_f = Quantity(type=float, description='softness of the frictional force')
    l_f = Quantity(type=float, description='length of particle in frictional force in Angstrom')
    d_f = Quantity(type=float, description='width of particle in frictional force in Angstrom')
    ld_f = Quantity(type=float, description='geometry ratio of particle in frictional force')
    esp0_f = Quantity(type=float, description='well-depth of cross configuration in frictional force in kcal/mol')
    esp_rat_f = Quantity(type=float, description='well-depth ratio of face-to-face and side-to-side in frictional force')
    relax = Quantity(type=float, description='timestep for relax of the system')
    ntime = Quantity(type=int, description='total number of timestep for the simulation')
    thermo = Quantity(type=float, description='frequency for screen print')
    dt = Quantity(type=float, description='timestep in femto-seconds')
    verlet = Quantity(type=float, description='factor for integration')
    atom_style =Quantity(type=str)
    pair_style = Quantity(type=str, shape=['*'], description='parameters needed for this force-field pair_style')
    pair_coeff = Quantity(type=str, shape=['*'], description='parameters needed for pair-wise interactions. First parameter is itype, second one is jtype.')
    integration = SubSection(sub_section=CGIntegration.m_def)
class CGCalculation(Calculation):
    m_def = Section(validate=False, extends_base_section=False)    

    settings = SubSection(sub_section=CGSettings.m_def)
    input = SubSection(sub_section=CGInput.m_def)
    vacf = SubSection(sub_section=CGVacf.m_def)


m_package.__init_metainfo__()
