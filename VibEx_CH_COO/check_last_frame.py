#!/bin/env python
import MDAnalysis as mda
import numpy as np
import os
import argparse
import sys

#parse command line arguments
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser._action_groups.pop()
required = parser.add_argument_group("required arguments")
required.add_argument("-i", "--input",   type=str,   help="the name of a dcd file",  required=True)
required = required.add_argument_group("optional arguments")
required.add_argument("-o", "--output",   type=str,   help="output file name",  default="spe_molpro.npz")
args = parser.parse_args()
dcd_file = args.input



# Define a function to check if the "result-test.dat" file exists
def check_results_exist(file_path):
    return os.path.exists(file_path)

# Define a function to calculate the distance between two centers
def calculate_distance(center1, center2):
    return np.linalg.norm(center1 - center2)

# Define a function to calculate the bond angle between three atoms
def calculate_bond_angle(atom1, atom2, atom3):
    vector1 = atom1.positions[0] - atom2.positions[0]
    vector2 = atom3.positions[0] - atom2.positions[0]
    cosine_theta = np.dot(vector1, vector2) / (np.linalg.norm(vector1) * np.linalg.norm(vector2))
    bond_angle = np.arccos(cosine_theta) * (180.0 / np.pi)  # Convert to degrees
    return bond_angle

# Function to process the last frame of a single DCD file
#def process_last_frame_dcd_files(start_num,end_num,num_hcooh,num_dio,num_explosion,num_nonreactive,num_unexpected):
def process_last_frame_dcd_files(dcd_file):
    psf_file = f"ch2oo-2.psf"
    flag_hcooh = False
    flag_dio = False
    flag_nr = False
    flag_explosion = False
    flag_unexpected = False
    flag_termination = False

    try:
        u = mda.Universe(psf_file, dcd_file)
    except IOError:
        print(f"Error: Unable to open DCD file '{dcd_file}'. Please check if the file exists.")

    if len(u.trajectory) < 1:
        print(f"Warning: The DCD file '{dcd_file}' does not contain any frames.")
        exit()

    last_frame = u.trajectory[-1]  # Access the last frame directly

    a = u.select_atoms("name C")
    b = u.select_atoms("name H")
    c = u.select_atoms("name H1")
    d = u.select_atoms("name O")
    e = u.select_atoms("name O1")

    coma = a.center_of_mass()
    comb = b.center_of_mass()
    comc = c.center_of_mass()
    comd = d.center_of_mass()
    come = e.center_of_mass()

    distance1 = calculate_distance(coma, comb) #CH
    distance2 = calculate_distance(coma, comc) #CH1
    distance3 = calculate_distance(coma, comd) #CO
    distance4 = calculate_distance(coma, come) #CO1
    distance5 = calculate_distance(comb, come) #HO1
    distance6 = calculate_distance(comd, come) #OO1
    distance7 = calculate_distance(comc, come) #O1H1
    distance8 = calculate_distance(comc, comd) #OH1
    distance9 = calculate_distance(comb, comd) #OH
    distance10 = calculate_distance(comb, comc) #HH1

    # Calculate bond angles
    angle_coo1 = calculate_bond_angle(a, d, e)
    angle_coo2 = calculate_bond_angle(a, e, d)

    if (
        (distance2 < 1.5 and # h1coo1h 
        distance1 > 1.8 and
        distance3 < 1.8 and
        distance6 < 2.3 and #oo1
        distance5 < 1.4 and #o1h
        angle_coo1 > 90.0 and
        distance4 >= 2.0) or

        (distance1 < 1.5 and # hcoo1h1 
        distance2 > 1.8 and
        distance3 < 1.8 and
        distance6 < 2.3 and
        distance7 < 1.4 and
        angle_coo1 > 90.0 and
        distance4 >= 2.0) or

        (distance2 < 1.5 and # h1co1oh 
        distance1 > 1.8 and
        distance4 < 1.8 and
        distance6 < 2.3 and #oo1
        distance9 < 1.4 and #o1h
        angle_coo2 > 90.0 and
        distance3 >= 2.0) or

        (distance1 < 1.5 and # hco1oh1
        distance2 > 1.8 and
        distance4 < 1.8 and
        distance6 < 2.3 and
        distance8 < 1.4 and
        angle_coo2 > 90.0 and
        distance3 >= 2.0) or

        (distance2 < 1.8 and # h1co+oh 
        distance3 < 1.6 and
        distance6 >= 2.3 and
        distance7 > 1.8 and
        distance5 < 1.4 and
        distance1 > 1.8 and
        distance4 > 1.8) or

        (distance1 < 1.6 and # hco+o1h1 
        distance3 < 1.8 and
        distance6 >= 2.3 and
        distance7 < 1.4 and
        distance5 > 1.8 and
        distance2 > 1.8 and
        distance4 > 1.8) or 

        (distance2 < 1.6 and # h1co + o1h 
        distance3 < 1.6 and
        distance8 < 1.5 and
        distance6 >= 1.9 and
        distance7 > 1.8 and
        distance5 < 1.8 and
        distance4 > 1.8) or 

        (distance1 < 1.6 and # hco1 + oh1 
        distance4 < 1.6 and
        distance8 < 1.6 and
        distance6 >= 1.9 and
        distance7 > 1.8 and
        distance3 > 1.8 and
        distance2 > 1.8) or 

        (distance2 < 1.6 and # h1co1 + oh 
        distance4 < 1.6 and
        distance9 < 1.5 and
        distance6 >= 1.9 and
        distance7 > 1.8 and
        distance1 > 1.8 and
        distance3 > 1.8) or 

        (distance1 > 1.8 and # co + ho1h1 
        distance2 > 1.8 and
        distance3 < 1.6 and
        distance4 > 2.5 and 
        distance5 < 1.4 and
        distance6 >= 1.9 and
        distance7 < 1.4 and
        distance9 > 1.8 ) or

        (distance1 > 1.8 and # co1 + hoh1 
        distance2 > 1.8 and
        distance4 < 1.6 and
        distance3 > 2.5 and 
        distance9 < 1.4 and
        distance6 >= 1.9 and
        distance8 < 1.4 and
        distance5 > 1.8 ) or

        (distance3 < 1.6 and  #CO2 + H2
        distance4 < 1.6 and
        distance1 > 1.9 and
        distance2 > 1.9 and
        distance5 > 1.9 and
        distance7 > 1.9 and
        distance8 > 1.9 and
        distance9 > 1.9 )
    ):
        flag_hcooh = True

    elif (
        distance1 < 1.5 and
        distance2 < 1.5 and
        distance3 < 1.8 and
        distance4 < 1.7 and
        distance6 < 2.0 and
        (angle_coo1 <= 90.0 or angle_coo2 <= 90.0)
    ):
        flag_dio = True

    elif (
        (min(distance1,distance5,distance9,distance10) > 2.8) or #H
        (min(distance2,distance7,distance8,distance10) > 2.8) or #H1
        (min(distance3,distance6,distance8,distance9) > 2.8) or #O
        (min(distance4,distance5,distance6,distance7) > 2.8) or #O1
        (min(distance1,distance2,distance3,distance4) > 2.8)  #C
    ):
        flag_explosion = True

    elif (
        distance1 < 1.8 and 
        distance2 < 1.8 and
        distance3 < 1.8 and
        distance6 < 2.3 and
        distance4 > 1.7 and
        angle_coo1 >= 90.0
    ):
        flag_nr = True
    else:
        flag_unexpected = True

    if (flag_hcooh or flag_dio or flag_explosion):
        flag_termination = True
    return flag_termination
# Loop through the given DCD file
flag_reactive=process_last_frame_dcd_files(dcd_file)
if flag_reactive :
    sys.exit(0)
else:
    sys.exit(-1)

