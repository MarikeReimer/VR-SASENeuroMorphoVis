####################################################################################################
# Copyright (c) 2016 - 2024, EPFL / Blue Brain Project
# Author(s): Marwan Abdellah <marwan.abdellah@epfl.ch>
#
# This file is part of NeuroMorphoVis <https://github.com/BlueBrain/NeuroMorphoVis>
#
# This program is free software: you can redistribute it and/or modify it under the terms of the
# GNU General Public License as published by the Free Software Foundation, version 3 of the License.
#
# This Blender-based tool is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.
# If not, see <http://www.gnu.org/licenses/>.
####################################################################################################

# System imports
import argparse

# Blender imports
from mathutils import Vector, Matrix
import neuron


####################################################################################################
# @parse_command_line_arguments
####################################################################################################
def parse_command_line_arguments(arguments=None):
    """Parses the input arguments.

    :param arguments:
        Command line arguments.
    :return:
        Arguments list.
    """

    # add all the options
    description = 'NeuroRender: a simple tool for rendering neurons generated by NeuroMorphoVis'
    parser = argparse.ArgumentParser(description=description)

    arg_help = 'Rendering configuration file'
    parser.add_argument('--config',
                        action='store', dest='config', help=arg_help)

    arg_help = 'Sytle file'
    parser.add_argument('--style-file',
                        action='store', dest='style_file', help=arg_help)

    arg_help = 'Input directory, where the meshes are stored'
    parser.add_argument('--input-directory',
                        action='store', dest='input_directory', help=arg_help)

    arg_help = 'Input data type: blend, ply, obj'
    parser.add_argument('--input-type',
                        action='store', dest='input_type', help=arg_help)

    arg_help = 'Base image resolution'
    parser.add_argument('--resolution',
                        action='store', default=512, dest='resolution', help=arg_help)

    arg_help = 'Projection, orthographic or perspective'
    parser.add_argument('--projection',
                        action='store', default='orthographic', dest='projection', help=arg_help)

    arg_help = 'Number of samples per pixel SPP'
    parser.add_argument('--spp',
                        action='store', default=32, dest='num_samples', help=arg_help)

    arg_help = 'Output directory, where the final image and scene will be stored'
    parser.add_argument('--output-directory',
                        action='store', dest='output_directory', help=arg_help)

    arg_help = 'Load positions and draw spheres'
    parser.add_argument('--use-spheres',
                        action='store_true', default=False, dest='use_spheres',help=arg_help)

    arg_help = 'Transform neurons to local positions'
    parser.add_argument('--transform',
                        action='store_true', default=False, dest='transform', help=arg_help)

    arg_help = 'Image and scene prefix'
    parser.add_argument('--prefix',
                        action='store', default='image', dest='prefix', help=arg_help)

    # Parse the arguments
    return parser.parse_args()


####################################################################################################
# @parse_rendering_configuration
####################################################################################################
def parse_rendering_configuration(configuration_file):
    """Parse a rendering configuration file.

    :param configuration_file:
        A given configuration file to load the neurons and rendering the scene.
    :return:
        A list of neurons.
    """

    # File contents
    config_data = list()

    # Read the configuration file
    config_file_handle = open(configuration_file, 'r')
    for line in config_file_handle:
        config_data.append(line)
    config_file_handle.close()

    # A list of all the parsed neurons
    neurons = list()

    # Parse the configuration
    index = 0
    while True:

        # We have reached the end of the file
        if index > len(config_data) - 1:
            break

        # Get a new line
        line = config_data[index]

        # If neuron is found
        if 'NEURON' in line:

            # Construct a new neuron object
            neuron_object = neuron.Neuron()

            # Parse neuron properties
            while True:

                # Get a new line
                index += 1

                if index > len(config_data) - 1:
                    break
                line = config_data[index]

                # GID
                if 'GID' in line:
                    line = line.replace('    ', '')
                    line = line.replace('\n', '')
                    line = line.split('GID: ')
                    neuron_object.gid = line[1]
                elif 'TAG' in line:
                    line = line.replace('    ', '')
                    line = line.replace('\n', '')
                    line = line.split('TAG: ')
                    neuron_object.tag = int(line[1])
                elif 'MORPHOLOGY_TYPE' in line:
                    line = line.replace('    ', '')
                    line = line.replace('\n', '')
                    line = line.split('MORPHOLOGY_TYPE: ')
                    neuron_object.mtype = line[1]
                elif 'MORPHOLOGY_LABEL' in line:
                    line = line.replace('    ', '')
                    line = line.replace('\n', '')
                    line = line.split('MORPHOLOGY_LABEL: ')
                    neuron_object.mlabel = line[1]
                elif 'POSITION' in line:
                    line = line.replace('    ', '')
                    line = line.replace('\n', '')
                    line = line.split('POSITION: ')[1]
                    line = line.split(' ')
                    x = float(line[0])
                    y = float(line[1])
                    z = float(line[2])
                    neuron_object.position = Vector((x, y, z))
                elif 'ORIENTATION' in line:
                    line = line.replace('    ', '')
                    line = line.replace('\n', '')
                    line = line.split('ORIENTATION: ')
                    neuron_object.orientation = float(line[1])
                elif 'COLUMN' in line:
                    line = line.replace('    ', '')
                    line = line.replace('\n', '')
                    line = line.split('COLUMN: ')
                    neuron_object.column = line[1]
                elif 'LAYER' in line:
                    line = line.replace('    ', '')
                    line = line.replace('\n', '')
                    line = line.split('LAYER: ')
                    neuron_object.layer = line[1]
                elif 'MEAN_RADIUS' in line:
                    line = line.replace('    ', '')
                    line = line.replace('\n', '')
                    line = line.split('MEAN_RADIUS: ')
                    neuron_object.soma_mean_radius = float(line[1])
                elif 'MIN_RADIUS' in line:
                    line = line.replace('    ', '')
                    line = line.replace('\n', '')
                    line = line.split('MIN_RADIUS: ')
                    neuron_object.soma_min_radius = float(line[1])
                elif 'MAX_RADIUS' in line:
                    line = line.replace('    ', '')
                    line = line.replace('\n', '')
                    line = line.split('MAX_RADIUS: ')
                    neuron_object.soma_max_radius = float(line[1])
                elif 'TRANSFORM' in line:
                    line = line.replace('    ', '')
                    line = line.replace('\n', '')
                    line = line.split('TRANSFORM: ')[1]
                    line = line.split(' ')
                    t = Matrix()
                    t[0][0:4] = float(line[0]), float(line[1]), float(line[2]), float(line[3])
                    t[1][0:4] = float(line[4]), float(line[5]), float(line[6]), float(line[7])
                    t[2][0:4] = float(line[8]), float(line[9]), float(line[10]), float(line[11])
                    t[3][0:4] = float(line[12]), float(line[13]), float(line[14]), float(line[15])
                    neuron_object.transform = t
                elif not line.strip():
                    # Neuron is done
                    break
                else:
                    # Unrecognized parameter
                    continue

            # Add the neurons to the list
            neurons.append(neuron_object)

            # Increment the index
            index += 1

    return neurons


####################################################################################################
# @parse_style_file
####################################################################################################
def parse_style_file(style_file):
    """Parse a style file for coloring the neurons.

    :param style_file:
        Input style file.
    :return:
        A style list.
    """

    style_map = list()

    # Read the configuration file
    style_file_handle = open(style_file, 'r')
    for line in style_file_handle:
        if '#' in line:
            continue
        elif not line.strip():
            continue
        else:
            line = line.replace('\n', '')
            line = line.split(' ')
            tag = int(line[0])
            r = float(line[1])
            g = float(line[2])
            b = float(line[3])
            if r > 1.0:
                r = r / 256.0
            if g > 1.0:
                g = g / 256.0
            if b > 1.0:
                b = b / 256.0
            alpha = float(line[4])
            color = Vector((r, g, b))
            shader = line[5]
            style_map.append([tag, color, alpha, shader])
    style_file_handle.close()

    return style_map
