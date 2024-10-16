# -*- coding: utf-8 -*-
####################################################################################################
# Copyright (c) 2020, EPFL / Blue Brain Project
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
import sys
import os
import argparse
import subprocess

# NeuroMorphoVis imports
import nmv.interface

sys.path.append(('%s/.' % (os.path.dirname(os.path.realpath(__file__)))))
sys.path.append(('%s/../' % (os.path.dirname(os.path.realpath(__file__)))))

# Blender imports
import analysis_input_vs_optimized


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
    description = 'This script takes the path to an input mesh and creates the corresponding ' \
                  'watertight mesh and the stats. of both meshes and creates a comparative result'
    parser = argparse.ArgumentParser(description=description)

    arg_help = 'Blender executable'
    parser.add_argument('--blender-executable',
                        action='store', dest='blender_executable', help=arg_help)

    arg_help = 'The input directory that contains all the meshes'
    parser.add_argument('--input-directory', action='store', help=arg_help)

    arg_help = 'Output directory, where the final image/movies and scene will be stored'
    parser.add_argument('--output-directory', action='store', help=arg_help)

    arg_help = 'Quality checker executable from Ultraliser'
    parser.add_argument('--quality-checker-executable', action='store', help=arg_help)

    arg_help = 'Number of parallel cores'
    parser.add_argument('--num-cores',
                        action='store', default=1, type=int, help=arg_help)

    # Parse the arguments
    return parser.parse_args()

####################################################################################################
# @construct_per_mesh_command
####################################################################################################
def construct_per_mesh_command(args, input_mesh):

    command = ''
    command += ' %s ' % args.blender_executable
    command += ' -b --verbose 0 --python create-original-vs-optimized-histograms-for-mesh.py -- '
    command += ' --input-directory %s ' % args.input_directory
    command += ' --input-mesh %s ' % input_mesh
    command += ' --quality-checker-executable %s ' % args.quality_checker_executable
    command += ' --output-directory %s ' % args.output_directory
    return command


####################################################################################################
# @execute_command
####################################################################################################
def execute_command(command):

    print(command)
    subprocess.call(command, shell=True)


####################################################################################################
# @execute_commands_parallel
####################################################################################################
def execute_commands_parallel(shell_commands, num_cores):

    from joblib import Parallel, delayed
    Parallel(n_jobs=num_cores)(delayed(execute_command)(i) for i in shell_commands)

####################################################################################################
# @ Main
####################################################################################################
if __name__ == "__main__":

    # Get all arguments after the '--'
    args = sys.argv
    sys.argv = args[args.index("--") + 0:]

    # Parse the command line arguments
    args = parse_command_line_arguments()

    # Create the output hierarchy
    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)
    intermediate_directory = '%s/intermediate-images' % args.output_directory
    if not os.path.exists(intermediate_directory):
        os.makedirs(intermediate_directory)
    images_directory = '%s/images' % args.output_directory
    if not os.path.exists(images_directory):
        os.makedirs(images_directory)
    scenes_directory = '%s/scenes' % args.output_directory
    if not os.path.exists(scenes_directory):
        os.makedirs(scenes_directory)

    original_meshes_directory = '%s/meshes' % args.input_directory

    # Get all the meshes in the path, either obj or ply
    list_meshes = nmv.file.get_files_in_directory(original_meshes_directory, file_extension='.obj')

    commands = list()
    for input_mesh in list_meshes:
        commands.append(construct_per_mesh_command(args=args, input_mesh=input_mesh))

    execute_commands_parallel(commands, num_cores=args.num_cores)
