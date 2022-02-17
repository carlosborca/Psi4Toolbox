#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
This is code is designed to extract molecular geometries out of a Psi4 
geometry optimization output file to create and XYZ trajectory file.
"""

#
# @BEGIN LICENSE
#
# Psi4Toolbox: Automated scripts to pre- and post-process Psi4 I/O files
#
# Copyright (c) 2022
# Carlos H. Borca
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4Toolbox.
#
# Psi4Toolbox is free software; you can redistribute it and/or modify
# it under the tesms of the GNU Lesser General Public License as
# published by the Free Software Foundation, version 3.
#
# Psi4Toolbox is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with Psi4Toolbox; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301 USA.
#
# @END LICENSE
#

__author__  = "Carlos H. Borca"
__version__ = "0.0.1"
__credits__ = ["Carlos H. Borca"]
__email__   = "carlosborca@gmail.com"
__license__ = "LGPLv3"

# Standard Python imports.
import argparse
import time
import sys

start = time.time() # Start the timer.

#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#
def createParser():
    """
    Creates a parser object from system argument variables.
                
    Returns
    -------
    parser : ArgumentPaser
        Parser object of argparse module
    """

    # For the argument parser. An example of usage and a description of what this script does.
    exas =  "python {} 'h2o.out' -t h2o.xyz".format(sys.argv[0])
    docs =  "Psi4Toolbox: Automated scripts to post process Psi4 output files"
    docs += "This is code is designed to extract molecular geometries out of a Psi4"
    docs += "geometry optimization output file to create and XYZ trajectory file."
    epis =  "Thank you for using the Psi4Toolbox."

    # Help strings.
    hg = "Psi4 geometry optimization output file name, i.e. 'h2o.out'"
    ht = "XYZ trjectory file name, i.e. 'h2o' (default: 'traj.xyz')"

    # Create parser.
    parser = argparse.ArgumentParser(description=docs, epilog=epis, usage=exas)

    # Positional Arguments: Must be specified.
    parser.add_argument("g", metavar="Input", type=str, help=hg)

    # Optional arguments the script might expect.
    parser.add_argument("-t", "--traj", dest="t", default="traj.xyz", help=ht)

    return parser
#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#

#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#
def processArgs(args):
    """
    Take the parsed arguments and return usable objects in a dictionary.
    
    Parameters
    ----------
    args : class 'argparse.Namespace'
        An object containing the runtime arguments
        
    Returns
    -------
    keywords : dict
        Set of runtime variables 
    """

    keywords = {} # A dictionary for variable objects.

    keywords['Psi4GeomOptOutFile'] = args.g
    keywords['XYZTrajectoryFile']  = args.t

    return keywords
#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#

#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#
def getLineMarksNumAtoms(psi4_go_out_file):
    """
    Extraction of beginning and end lines indices enclosing geometries.
        
    Parameters
    ----------
    psi4_go_out_file : str
        Name of a Psi4 geometry optimization output file
        
    Returns
    -------
    number_of_atoms : int
        Number of atoms in the molecular geometry
    cart_geom_lines : list
        List of line numbers where geometries start in the Psi4 output
    """
    
    f = open(psi4_go_out_file, "r")
    
    cart_geom_lines = []
    optk_fini_lines = []
    
    i = 0
    
    for line in f:
        
        if "Cartesian Geometry (in Angstrom)" in line:
            
            cart_geom_lines.append(i)
            
        if "OPTKING Finished Execution" in line:
            
            optk_fini_lines.append(i)
        
        i += 1
    
    number_of_atoms = optk_fini_lines[0] - cart_geom_lines[0] - 2
    
    return number_of_atoms, cart_geom_lines
#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#

#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#
def formatTraj(psi4_go_out_file, number_of_atoms, cart_geom_lines):
    """
    Format the trajectory file scrapping geometries from the Psi4 output.
    
    Parameters
    ----------
    psi4_go_out_file : str
        Name of a Psi4 geometry optimization output file
    number_of_atoms : int
        Number of atoms in the molecular geometry
    cart_geom_lines : list
        List of line numbers where geometries start in the Psi4 output

    Returns
    -------
    traj_lin_list : list
        List of lines to write the XYZ trajectory
    """
    
    list_of_lines = [l for l in open(psi4_go_out_file, "r")]
    traj_lin_list = [] # List to create the output file from.
    
    frame = 0  # Trajectory frames counter. 
    
    for geom in cart_geom_lines:
           
        traj_lin_list.append("{}".format(number_of_atoms)) # Number atoms line     
        traj_lin_list.append("File: {}, Frame: {}".format(psi4_go_out_file, frame)) # Comment line
    
        for atom in range(number_of_atoms): # Coordinates lines.
            traj_lin_list.append(list_of_lines[geom + 1 + atom].strip())
            
        frame += 1
        
    return traj_lin_list
#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#

#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#
def writeTrajFromList(trajec_file_name, traj_lin_list):
    """
    Write the list of lines with the trajectory to an output file.
    
    Parameters
    ----------
    trajec_file_name
        Name of the XYZ trajectory file written as output
    traj_lin_list : list
        List of lines to write the XYZ trajectory
    """
    
    trjf = open(trajec_file_name, "w")
    
    for l in traj_lin_list:
        
        trjf.writelines(l + "\n")
    
    trjf.close()
    
    return
#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#

#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#
def main(argv):
    """
    This is code is designed to extract molecular geometries out of a Psi4 
    geometry optimization output file to create and XYZ trajectory file.
    """
        
    # Check arguments given at execution:
    print("No arguments passed. Check documentation: `python Psi4OptimizeToXYZ.py --help`") if (not argv) else ""
    
    # Call parser creator to process options introduced by the user at execution.
    parser = createParser()
    args = parser.parse_args()
    
    # Process parsed arguments into usable objects.
    keywords = processArgs(args)
    
    # Extract the number of atoms and the first and last line of geometry printouts.
    number_of_atoms, cart_geom_lines = getLineMarksNumAtoms(keywords['Psi4GeomOptOutFile'])
    
    # Scrape the geometries out of the Psi4 output.
    traj_lin_list = formatTraj(keywords['Psi4GeomOptOutFile'], number_of_atoms, cart_geom_lines)
    
    # Write the trajectory out to an XYZ file.
    writeTrajFromList(keywords['XYZTrajectoryFile'], traj_lin_list)
    
    # Determine execution time and print exit message.
    finish = time.time() - start
    
    print("Execution terminated.")
    print("Time elapsed = {:.2f} s.".format(finish))
    
    return
#--------+---------+---------+---------+---------+---------+---------+-=---=---+---------+---------+---------+---------#

# Main code execution
if __name__ == "__main__":
    main(sys.argv[1:])