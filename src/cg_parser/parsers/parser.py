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

import yaml
import os
import re
import datetime
import numpy as np
import filecmp
from pathlib import Path

from nomad.config import config
from nomad.datamodel.metainfo.workflow import Workflow
from nomad.parsing.parser import MatchingParser
from runschema.run import Run, Program
from runschema.calculation import Calculation
from cg_parser.schema_packages.schema_package import CGCalculation, CGSettings, CGIntegration, CGInput, CGBox_size, CGEllipsoid_properties, CGMolecule_velocities, CGVacf

configuration = config.get_plugin_entry_point(
    'cg_parser.parsers:parser_entry_point'
)




def DetailedCGParser(filepath, archive):
    run = Run()
    archive.run.append(run)
    

    calculation = CGCalculation()
    run.calculation.append(calculation)
    
    with open(str(filepath.parent) + r'/in.*cg.lmp') as in_file:
        cgsettings = CGSettings()
        calculation.settings = cgsettings
        
        for i, line in enumerate(in_file):
            line = line.replace('#', '').strip('\n').strip(' ').strip('\t')
            if i==4:
                run.program = Program(name="Ka Chun " + line)
#            if re.search(r'case', line.lower()):
#                parts = line.split(': ')
#                sec_material = archive.m_setdefault("results.material")
#                sec_material.material_name = parts[1]
            if re.search(r'ndim\s+equal', line):
                parts = line.replace(' ', '').split('equal')
                
                cgsettings.ndim = int(parts[1])
            if re.search(r'temp\s+equal', line.lower()):
                parts = line.replace(' ', '').split('equal')
                calculation.temperature = float(parts[1])
            if re.search(r'seed\s+equal', line.lower()):
                parts = line.replace(' ', '').split('equal')
                cgsettings.seed = float(parts[1])
            if re.search(r'rc_global\s+equal', line.lower()):
                parts = line.replace(' ', '').split('equal')
                cgsettings.rc_global = float(parts[1])
            if re.search(r'skin\s+equal', line.lower()):
                parts = line.replace(' ', '').split('equal')
                cgsettings.skin = str(parts[1])
            if re.search(r'dw_c\s+equal', line.lower()):
                parts = line.replace(' ', '').split('equal')
                cgsettings.dw_c = float(parts[1])
            if re.search(r'l_c\s+equal', line.lower()):
                parts = line.replace(' ', '').split('equal')
                cgsettings.l_c = float(parts[1])
            if re.search(r'd_c\s+equal', line.lower()):
                parts = line.replace(' ', '').split('equal')
                cgsettings.d_c = float(parts[1])
            if re.search(r'esp0_c\s+equal', line.lower()):
                parts = line.replace(' ', '').split('equal')
                cgsettings.esp0_c = float(parts[1])
            if re.search(r'esp_rat_c\s+equal', line.lower()):
                parts = line.replace(' ', '').split('equal')
                cgsettings.esp_rat_c = float(parts[1])
            if re.search(r'dw_f\s+equal', line.lower()):
                parts = line.replace(' ', '').split('equal')
                cgsettings.dw_f = float(parts[1])
            if re.search(r'l_f\s+equal', line.lower()):
                parts = line.replace(' ', '').split('equal')
                cgsettings.l_f = float(parts[1])
            if re.search(r'd_f\s+equal', line.lower()):
                parts = line.replace(' ', '').split('equal')
                cgsettings.d_f = float(parts[1])
            if re.search(r'ld_f\s+equal', line.lower()):
                parts = line.replace(' ', '').split('equal')
                cgsettings.ld_f = float(parts[1])
            if re.search(r'esp0_f\s+equal', line.lower()):
                parts = line.replace(' ', '').split('equal')
                cgsettings.esp0_f = float(parts[1])
            if re.search(r'esp_rat_f\s+equal', line.lower()):
                parts = line.replace(' ', '').split('equal')
                cgsettings.esp_rat_f = float(parts[1])
            if re.search(r'relax\s+equal', line.lower()):
                parts = line.replace(' ', '').split('equal')
                cgsettings.relax = float(parts[1])
            if re.search(r'ntime\s+equal', line.lower()):
                parts = line.replace(' ', '').split('equal')
                calculation.step = int(parts[1])
                
            if re.search(r'thermo\s+equal', line.lower()):
                parts = line.replace(' ', '').split('equal')
                cgsettings.thermo = float(parts[1])
            if re.search(r'dt\s+equal', line.lower()):
                parts = line.replace(' ', '').split('equal')
                
                cgsettings.dt = float(parts[1])
            if re.search(r'verlet\s+equal', line.lower()):
                parts = line.replace(' ', '').split('equal')
                cgsettings.verlet = float(parts[1])
            if line.startswith(r'atom_style'):
                
                parts =  line.replace(' ', '').replace('\t', '').split('style')
                
                cgsettings.atom_style = parts[1]
            if line.startswith(r'pair_style'):
                parts1, parts2 = line.split('style')
                parts = parts2.replace('{', '').replace('}', '').replace(' ', '').split('$')
                cgsettings.pair_style = []
                for i in range(len(parts)):
                    cgsettings.pair_style.append(parts[i])
            if line.startswith(r'pair_coeff'):
                cgsettings.pair_coeff = []
                line = ' '.join(line.split())
                parts = line.replace('$', '').replace('{', '').replace('}', '').split(' ')
                for i in range(len(parts)-1):
                    cgsettings.pair_coeff.append(parts[i+1])
            if re.search(r'integrate', line):
                line  = line.split()
                _integration_name = []
                
                cgintegration = CGIntegration()
                settings.integration = cgintegration
                cgintegration.parameters = []
                j = 0
                
                for word in line:
                    if 'integrate' in word:
                        j = 1
                        continue
                    if '$' in word:
                        cgintegration.parameters.append(word.replace('$', '').replace('{',  '').replace('}', ''))
                        
                    if j>0 and '$' not in word:
                        _integration_name.append(word)
                cgintegration.name = ' '.join(_integration_name)
                                
    with open(str(filepath.parent) + r'/data_quat.in') as data_file:
        cginput = CGInput()
        calculation.input = cginput
        
        cgbox_size = CGBox_size()
        cginput.box_size = cgbox_size

        cgellipsoid_properties = CGEllipsoid_properties()
        cginput.ellipsoid_properties = cgellipsoid_properties

        cgmolecule_velocities = CGMolecule_velocities()
        cginput.molecule_velocities = cgmolecule_velocities

        j=0
        atom_type = []
        ellipsoid_flag = []
        density_of_particle = []
        
        atom_x_positions = []
        atom_y_positions = []
        atom_z_positions = []
        length = []
        width = []
        height = []
        w_orientation = []
        x_orientation = []
        y_orientation = []
        z_orientation = []
        v_x = []
        v_y = []
        v_z = []
        w_x = []
        w_y = []
        w_z = []
        for i, line in enumerate(data_file):
            if re.search(r'\s+atoms', line):
                parts = line.split()
                cginput.n_atoms = int(parts[0])
                atom_positions_array = np.zeros((int(parts[0]), 3))
            if re.search(r'\s+ellipsoids', line):
                parts = line.split()
                cginput.n_ellipsoids = int(parts[0])
            if re.search(r'\s+bonds', line):
                parts = line.split()
                cginput.n_bonds = int(parts[0])    
            if re.search(r'\s+angles', line):
                parts = line.split()
                cginput.n_angles = int(parts[0])
            if re.search(r'\s+dihedrals', line):
                parts = line.split()
                cginput.n_dihedrals = int(parts[0])
            if re.search(r'\s+impropers', line):
                parts = line.split()
                cginput.n_impropers = int(parts[0])
            if re.search(r'\s+atom\s+types', line):
                parts = line.split()
                cginput.n_atom_types = int(parts[0])
            if re.search(r'xlo', line):
                parts = line.split()
                cgbox_size.xlo = float(parts[0])
            if re.search(r'xhi', line):
                parts  = line.split()
                cgbox_size.xhi = float(parts[1])
            if re.search(r'ylo', line):
                parts = line.split()
                cgbox_size.ylo = float(parts[0])
            if re.search(r'yhi', line):
                parts = line.split()
                cgbox_size.yhi  = float(parts[1])
            if re.search(r'zlo', line):
                parts = line.split()
                cgbox_size.zlo = float(parts[0])
            if re.search(r'zhi', line):
                parts = line.split()
                cgbox_size.zhi = float(parts[1])           
            if re.search(r'Atoms', line):
                j = 1
                continue

            if j == 1 and len(line.split(' ')) > 1:
                parts = line.split()

                atom_type.append(int(parts[1]))
                atom_type_array = np.array(atom_type)
                cginput.atom_type = atom_type_array
                
                ellipsoid_flag.append(int(parts[2]))
                ellipsoid_flag_array = np.array(ellipsoid_flag)
                cginput.ellipsoid_flag = ellipsoid_flag_array

                density_of_particle.append(float(parts[3]))
                density_of_particle_array = np.array(density_of_particle)
                cginput.density_of_particle = density_of_particle_array
                
                atom_x_positions.append(np.float64(parts[4]))
                atom_y_positions.append(np.float64(parts[5]))
                atom_z_positions.append(np.float64(parts[6]))

            if re.search(r'Ellipsoids', line):
                j =  2
                continue
            if j == 2 and len(line.split(' ')) > 1:
                parts = line.split()
                
                length.append(float(parts[1]))
                length_array = np.array(length)
                cgellipsoid_properties.length = length_array

                width.append(float(parts[2]))
                width_array = np.array(width)
                cgellipsoid_properties.width = width_array

                height.append(float(parts[3]))
                height_array = np.array(height)
                cgellipsoid_properties.height = height_array

                w_orientation.append(float(parts[4]))
                w_orientation_array = np.array(w_orientation)
                cgellipsoid_properties.w_orientation = w_orientation_array
                
                x_orientation.append(float(parts[5]))
                x_orientation_array = np.array(x_orientation)
                cgellipsoid_properties.x_orientation = x_orientation_array
                
                y_orientation.append(float(parts[6]))
                y_orientation_array = np.array(y_orientation)
                cgellipsoid_properties.y_orientation =  y_orientation_array
                
                z_orientation.append(float(parts[7]))
                z_orientation_array = np.array(z_orientation)
                cgellipsoid_properties.z_orientation = z_orientation_array

            if re.search(r'Velocities', line):
                j = 3
                continue
            if j == 3 and len(line.split(' ')) > 1:
                parts = line.split()
                v_x.append(float(parts[1]))
                v_x_array = np.array(v_x)
                cgmolecule_velocities.v_x = v_x_array
                
                v_y.append(float(parts[2]))
                v_y_array = np.array(v_y)
                cgmolecule_velocities.v_y = v_y_array

                v_z.append(float(parts[3]))
                v_z_array = np.array(v_z)
                cgmolecule_velocities.v_z = v_z_array

                w_x.append(float(parts[4]))
                w_x_array = np.array(w_x)
                cgmolecule_velocities.w_x = w_x_array
                
                w_y.append(float(parts[5]))
                w_y_array = np.array(w_y)
                cgmolecule_velocities.w_y = w_y_array

                w_z.append(float(parts[6]))
                w_z_array = np.array(w_z)
                cgmolecule_velocities.w_z = w_z_array

        atom_positions_array[:, 0] = np.array(atom_x_positions)
        atom_positions_array[:, 1] = np.array(atom_y_positions)
        atom_positions_array[:, 2] = np.array(atom_z_positions)
        cginput.atom_coordinates = atom_positions_array
        
    with open(str(filepath.parent) + r'/vacf_atom.dat') as vacf_file:
        cgvacf = CGVacf()
        calculation.vacf = cgvacf
        
        _vacf_time = []
        _vacf= []

        for i, line in enumerate(vacf_file):
            if re.search(r'variables', line.lower()):
                continue
            parts = line.split()
            _vacf_time.append(float(parts[0]))
            _vacf_time_array = np.array(_vacf_time)
                        
            _vacf.append(float(parts[1]))        
            _vacf_array = np.array(_vacf)
        
        cgvacf.vacf = _vacf_array
        cgvacf.vacf_time = _vacf_time_array


class NewParser(MatchingParser):
    def parse(
        self,
        mainfile: str,
        archive: 'EntryArchive',
        logger: 'BoundLogger',
        child_archives: dict[str, 'EntryArchive'] = None,
    ) -> None:
        logger.info('NewParser.parse', parameter=configuration.parameter)

        archive.workflow2 = Workflow(name='test')
        mainfile = Path(mainfile)
        DetailedCGParser(mainfile, archive)