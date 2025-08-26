#!/bin/env python
import MDAnalysis as mda
import numpy as np
import os

# Define a function to check if the results-test.dat file exists
def check_results_exist(file_path):
    """
    Check if a file exists.

    Args:
        file_path (str): The path to the file to be checked.

    Returns:
        bool: True if the file exists, False otherwise.
    """
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

# Function to process multipole DCD files
def process_dcd_files(start_num, end_num, frame_interval,num_intervals):

    num_hcooh =[0] * num_intervals
    num_dio = [0]*num_intervals
    num_explosion = [0]*num_intervals
    num_unexpected = [0]*num_intervals
    num_nonreactive = 0
    with open("results.dat",'w') as outfile, open("reactiontime.dat","w") as ff:
        for num in range(start_num, end_num + 1):
            print(num) 

            psf_file = f"ch2oo-2.psf"
            # find the last dcd file which is the reacted or clapsed trajectory
            for kk in range(10,0,-1):
                dcd_file = f"nve_seed{num}-{kk}.dcd"
                if check_results_exist(dcd_file) :
                    break

            #check the last frame of the last dcd file for single traj, e.g. nve_seed1-10.dcd
            try:
                u = mda.Universe(psf_file, dcd_file)
            except IOError:
                print(f"Error: Unable to open DCD file '{dcd_file}'. Please check if the file exists.")
                continue

            last_frame = u.trajectory[-1]  # Access the last frame directly

            a = u.select_atoms("name C")
            b = u.select_atoms("name H")
            c = u.select_atoms("name H1")
            d = u.select_atoms("name O")
            e = u.select_atoms("name O1")
            #center of mass
            coma = a.center_of_mass()
            comb = b.center_of_mass()
            comc = c.center_of_mass()
            comd = d.center_of_mass()
            come = e.center_of_mass()
            # intermolecular distance
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
            angle_coo1 = calculate_bond_angle(a, d, e) #COO1
            angle_coo2 = calculate_bond_angle(a, e, d) #CO1O
            flag_hcooh = False
            flag_dio = False
            flag_nr = False
            flag_explosion = False
            flag_unexpected = False

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
                distance1 < 1.5 and #dioxirane
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

            # go to next traj if not reacted for the last frame of the last dcd file for one tr
            if flag_nr:
                outfile.write(f"nonreactive={num} {kk} {distance3} {distance5} {distance4}\n")
                num_nonreactive += 1
                continue # go to next traj

            #if the last frame of the last dcd file is not non-reactive
            for jj in range(1,11):#loop through dcd files in a reverse order
                dcd_file = f"nve_seed{num}-{jj}.dcd"
                #
                if not check_results_exist(dcd_file):
                   break 

                #read the every dcd file for every single trajectory
                try:
                    u = mda.Universe(psf_file, dcd_file)
                except IOError:
                    print(f"Error: Unable to open DCD file '{dcd_file}'. Please check if the file exists.")
                    continue

                if jj == 1:
                    t1 = int(len(u.trajectory) / 100) #simulation time for the first dcd file

                last_frame = u.trajectory[-1]  # Access the last frame directly

                a = u.select_atoms("name C")
                b = u.select_atoms("name H")
                c = u.select_atoms("name H1")
                d = u.select_atoms("name O")
                e = u.select_atoms("name O1")
                #center of mass
                coma = a.center_of_mass()
                comb = b.center_of_mass()
                comc = c.center_of_mass()
                comd = d.center_of_mass()
                come = e.center_of_mass()
                # intermolecular distance
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
                angle_coo1 = calculate_bond_angle(a, d, e) #COO1
                angle_coo2 = calculate_bond_angle(a, e, d) #CO1O
                flag_hcooh = False
                flag_dio = False
                flag_nr = False
                flag_explosion = False
                flag_unexpected = False

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
                    distance1 < 1.5 and #dioxirane
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
                    continue
                else:
                    flag_unexpected = True

                # go to next traj if not reacted for the last frame of the last dcd file for one tr
                #if flag_nr:
                #    if((jj == 1 and len(u.trajectory) == 10000000) or (len(u.trajectory) == 1000000 and ( jj == 7 or jj == 10 ))) :  
                #        outfile.write(f"nonreactive={num} {jj} {distance3} {distance5} {distance4}\n")
                #        num_nonreactive += 1
                #        break # go to next traj
                #    else:
                #        continue # go to next dcd file

                #if reacted, loop through this dcd file
                num_frames = min(len(u.trajectory),(num_intervals * frame_interval))
                ist = t1 // 10 
                flag_hcooh = False
                flag_dio = False
                flag_nr = False
                flag_explosion = False
                flag_unexpected = False
                for i in range(0,num_frames,frame_interval):
                    start_frame = i
                    end_frame = min(i + frame_interval, num_frames)

                    if jj > 1:
                        idx = i // frame_interval + ist + 1 
                    elif jj == 1:
                        idx = i // frame_interval

#                    print("start_frame="+str(i)+" end_frame="+str(end_frame))
                    for ii in range(start_frame,end_frame):
                        ts = u.trajectory[ii]

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

                        if(
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
                            num_hcooh[idx] += 1
                            outfile.write(f"hcooh={num} {jj} {ii} {distance3} {distance5} {distance4}\n")
                            flag_hcooh = True
                            break

                        elif (
                            distance1 < 1.5 and
                            distance2 < 1.5 and
                            distance3 < 1.8 and
                            distance4 < 1.7 and
                            distance6 < 2.0 and
                            (angle_coo1 <= 90.0 or angle_coo2 <= 90.0)
                        ):
                            num_dio[idx] += 1
                            outfile.write(f"dioxirane={num} {jj} {ii} {distance3} {distance5} {distance4}\n")
                            flag_dio = True
                            if jj == 1:
                                rtime = (ii + 1) * 0.01 # ps
                                ff.write(f"{rtime:.3f} \n")
                            else:
                                rtime = (jj - 1) * 100.0 + (ii + 1) * 0.01 #ps
                                ff.write(f"{rtime:.3f} \n")
                            break

                        elif (
                            (min(distance1,distance5,distance9,distance10) > 2.8) or #H
                            (min(distance2,distance7,distance8,distance10) > 2.8) or #H1
                            (min(distance3,distance6,distance8,distance9) > 2.8) or #O
                            (min(distance4,distance5,distance6,distance7) > 2.8) or #O1
                            (min(distance1,distance2,distance3,distance4) > 2.8)  #C
                        ):
                            num_explosion[idx] += 1
                            outfile.write(f"explosion={num} {jj} {ii} {distance3} {distance5} {distance4}\n")
                            flag_explosion = True
                            break

                        elif (
                            distance1 < 1.8 and
                            distance2 < 1.8 and
                            distance3 < 1.8 and
                            distance6 < 2.3 and
                            distance4 > 1.7 
                        ):
                            flag_nr = True 
                            continue
                        else:
                            num_unexpected[idx] += 1
                            outfile.write(f"unexpected={num} {jj} {ii} {distance3} {distance5} {distance4}\n")
                            flag_unexpected = True
                            break

                    if (flag_hcooh == True or flag_dio == True or flag_explosion == True or flag_unexpected == True):
                        break
                if (flag_hcooh == True or flag_dio == True or flag_explosion == True or flag_unexpected == True):
                  break

        for x in range(1,num_intervals):
            num_hcooh[x] += num_hcooh[x-1]
            num_dio[x] += num_dio[x-1]
            num_explosion[x] += num_explosion[x-1]
            num_unexpected[x] += num_unexpected[x-1]

    return num_hcooh, num_dio, num_explosion, num_nonreactive, num_unexpected

# Check if the results file exists and remove it if it does
file_path = "result-test.dat"
if check_results_exist(file_path):
    print(f"The file '{file_path}' exists and will be removed.")
    os.remove(file_path)

# Define the range of file numbers to process
start_num =1 
end_num = 2000  # Adjust the end number as needed
frame_interval = 1000  # Interval to process frames
num_intervals = 100

time_list = [(x+1)*(frame_interval//100) for x in range(num_intervals)]
# Loop through the range of file numbers and process frames at the specified interval

num_hcooh, num_dio, num_explosion, num_nonreactive, num_unexpected = process_dcd_files(start_num, end_num, frame_interval,num_intervals)
num_reactive=[0]*num_intervals
num_dcd = end_num + 1 - start_num
for tmp in range(num_intervals):
    num_reactive[tmp] = num_hcooh[tmp] + num_dio[tmp]

with open('traj-out.dat',"w") as ff:
    ff.write(f"The num of dcd files: {num_dcd}\n")
    ff.write(f"The num of reactive trajs: {num_reactive[-1]}\n")
    ff.write(f"The num of hcooh trajs: {num_hcooh[-1]}\n")
    ff.write(f"The num of dioxirane trajs: {num_dio[-1]}\n")
    ff.write(f"The num of nonreactive trajs: {num_nonreactive}\n")
    ff.write(f"The num of exploded trajs: {num_explosion[-1]}\n")
    ff.write(f"The num of unexpected trajs: {num_unexpected[-1]}\n")

with open("time-tot.dat",'w') as ft:
    for x in range(num_intervals):
        ft.write(f"{time_list[x]}\t{num_hcooh[x]}\t{num_dio[x]}\t{num_reactive[x]} \n")

print(f"The num of dcd files: {num_dcd}\n")
print(f"The num of reactive trajs: {num_reactive[-1]}\n")
print(f"The num of hcooh trajs: {num_hcooh[-1]}\n")
print(f"The num of dioxirane trajs: {num_dio[-1]}\n")
print(f"The num of nonreactive trajs: {num_nonreactive}\n")
print(f"The num of exploded trajs: {num_explosion[-1]}\n")
print(f"The num of unexpected trajs: {num_unexpected[-1]}\n")

