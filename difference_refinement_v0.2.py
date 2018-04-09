#Script for performing difference refinement with a "ground state" model, and "ground state" and "excited state" data sets (MTZ)
#Create a working directory, and put your two MTZ files and one PDB file in the directory
#You must have phenix and ccp4 installed and sourced
#Run with libtbx.python, pass the filenames as arguments in the following order: ground state model, ground state mtz, excited state mtz
#The ground state MTZ file must contain R-free flags
#Example usage: libtbx.python difference_refinement.py ground_coords.pdb ground_data.mtz excited_data.mtz

import os
import sys
import math
import string
from iotbx.reflection_file_reader import any_reflection_file
import subprocess


#set up directories for inputs and outputs

cwd = os.getcwd()

def specify_diff_refinement_input_files( starting_directory ):
  print "Collecting input from command line arguments..."
  ground_coords_fullpath = os.path.join(starting_directory, sys.argv[1])
  ground_data_fullpath = os.path.join(starting_directory, sys.argv[2])
  excited_data_fullpath = os.path.join(starting_directory, sys.argv[3])
  print" Input collected!"



def prepare_directories_for_diff_refinement ( starting_directory ):
  print "Preparing directories for Difference Refinement setup..."
  starting_structure_factor_mtz_path = os.path.join(starting_directory, "starting_structure_factor_mtz")
  os.mkdir(starting_structure_factor_mtz_path, 0755)
  scaled_intensity_mtz_path = os.path.join(starting_directory, "scaled_intensity_mtz")
  os.mkdir(scaled_intensity_mtz_path)
  scaled_structure_factor_mtz_path = os.path.join(starting_directory, "scaled_structure_factor_mtz")
  os.mkdir(scaled_structure_factor_mtz_path)
  difference_refinement_path = os.path.join(starting_directory, "difference_refinement")
  os.mkdir(difference_refinement_path, 0755)
  print "Preparation complete!"


def generate_I_model( coords_fullpath, I_model_output_dir_fullpath ):
  print "Using phenix.fmodel to generate Icalc for input model: %s" % (coords_fullpath)
  get_I_model_command = "phenix.fmodel %s output.type=real output.obs_type=intensities output.label=I high_resolution=1.5 add_sigmas=True scale=0.1 output.file_name=ground_state_f_model_as_intensity.mtz" % (coords_fullpath)
  os.system(get_I_model_command)
  print 

 
def scale_intensities_to_I_model( unscaled_intensity_mtz_fullpath, I_model_mtz_fullpath, scaling_directory_fullpath, starting_directory ):
  unscaled_intensity_mtz_filename = unscaled_intensity_mtz_fullpath.split('/')[-1]
  move_unscaled_intensity_mtz_file_to_scaling_directory = "mv %s %s" % (unscaled_intensity_mtz_fullpath, scaling_directory_fullpath)
  os.system(move_unscaled_intensity_mtz_file_to_scaling_directory)
  I_model_mtz_filename = I_model_mtz_fullpath.split('/')[-1]
  move_I_model_file_to_scaling_directory = "mv %s %s" % (I_model_mtz_fullpath, scaling_directory_fullpath)
  os.system(move_I_model_file_to_scaling_directory)
  os.chdir(scaling_directory_fullpath)
  scaled_output_filename_prefix = unscaled_mtz_filename.split('.')[-2]
  scale_to_I_model_command = "phenix.scale_and_merge reference_data=%s data=%s anomalous=False output.file_name=%s_scaled.mtz" % (starting_structure_factor_mtz_path, starting_structure_factor_mtz_path)
  os.system(scale_to_I_model_command)
  os.chdir(starting_directory)


def convert_scaled_intensities_to_structure_factors( scaled_intensity_mtz_fullpath, scaled_structure_factor_directory_fullpath, starting_directory ):
  scaled_intensity_mtz_filename = scaled_intensity_mtz_fullpath.split('/')[-1]
  scaled_structure_factor_mtz_prefix = scaled_intensity_mtz_filename.split('.')[-2]

  hkl_file = any_reflection_file(scaled_intensity_mtz_fullpath)
  info = hkl_file.file_content()
  xtals =info.crystals()
  if len(xtals) > 1:
    xtal = xtals[0]
  else:
    xtal =xtals
  space_group = info.space_group_number()
  unit_cell = xtal.unit_cell_parameters()
  uc_a, uc_b, uc_c, uc_al, uc_be, uc_ga = unit_cell

  uc_str = str(unit_cell)[1:-1]

  scaled_structure_factor_mtz_filename = scaled_structure_factor_mtz_prefix+"_asF.mtz"
  scaled_structure_factor_mtz_fullpath = os.path.join(scaled_structure_factor_directory_fullpath, scaled_structure_factor_mtz_filename)

  RFE_def_file_fullpath = os.path.join(starting_directory, "reflection_file_editor.def")
  RFE_def_file = open(RFE_def_file_fullpath, 'w')
 
  RFE_def_text = """show_arrays = False
dry_run = False
verbose = True
mtz_file {
  crystal_symmetry {
    unit_cell = %s
    space_group = %s
    output_unit_cell = None
    output_space_group = None
    change_of_basis = None
    eliminate_invalid_indices = False
    expand_to_p1 = False
    disable_unit_cell_check = True
    disable_space_group_check = False
    eliminate_sys_absent = True
  }
  d_max = None
  d_min = None
  wavelength = None
  output_file = %s
  job_title = None
  resolve_label_conflicts = True
  exclude_reflection = None
  miller_array {
    file_name = %s
    labels = Iobs,SIGIobs
    output_labels = Fobs SIGFobs
    column_root_label = None
    d_min = None
    d_max = None
    output_as = auto intensities *amplitudes_fw amplitudes
    anomalous_data = *Auto merged anomalous
    force_type = *auto amplitudes intensities
    scale_max = None
    scale_factor = None
    remove_negatives = False
    massage_intensities = False
    filter_by_signal_to_noise = None
    add_b_iso = None
    add_b_aniso = 0 0 0 0 0 0
    shuffle_values = False
    reset_values_to = None
   } 
  r_free_flags {
    generate = False
    force_generate = False
    new_label = "FreeR_flag"
    fraction = 0.1
    max_free = 2000
    lattice_symmetry_max_delta = 5
    use_lattice_symmetry = True
    use_dataman_shells = False
    n_shells = 20
    random_seed = None
    extend = True
    old_test_flag_value = None
    export_for_ccp4 = False
    preserve_input_values = True
    warn_if_all_same_value = True
    adjust_fraction = False
    d_eps = 0.0001
    relative_to_complete_set = False
    remediate_mismatches = False
  }
}

""" % (uc_str, space_group, output_mtz_file, input_intensity_mtz)
  
  RFE_def_file.write(RFE_def_text)
  RFE_def_file.close()

  print "Wrote RFE def file successfully! - Now converting intensities to structure factor amplitudes with RFE."

  subprocess.call(['iotbx.reflection_file_editor', 'reflection_file_editor.def'])

  print "Conversion to Structure Factors is finished!"


#generate f_model mtz from ground state model using phenix.refine and clean up unnecessary outputs



def convert_scaled_intensities_to_structure_factors_FandW( scaled_intensity_mtz_fullpath, scaled_structure_factor_directory_fullpath, starting_directory ):


#use phenix.mtz.dump to make mtzs readable as text and store F and sigF values in hkl-keyed dictionaries

ground_scaled_F_mtz = os.path.join(cwd, "scaled_starting_structure_factor_mtz/ground_state_SF_amplitudes_scaled.mtz")
excited_scaled_F_mtz = os.path.join(cwd, "scaled_starting_structure_factor_mtz/excited_state_SF_amplitudes_scaled.mtz")
f_model_mtz = os.path.join(cwd, "starting_structure_factor_mtz/f_model.mtz")

for sfa_mtz in [ground_scaled_F_mtz, excited_scaled_F_mtz, f_model_mtz]:
  if sfa_mtz == ground_scaled_F_mtz:
    output_prefix = "ground_scaled_F"
  if sfa_mtz == excited_scaled_F_mtz:
    output_prefix = "excited_scaled_F"
  if sfa_mtz == f_model_mtz:
    output_prefix = "f_model"
  mtz_dump_command = "phenix.mtz.dump -c -f s %s > %s.mtz.txt" % (sfa_mtz, output_prefix)
  os.system(mtz_dump_command)

ground_scaled_F_as_text_filename = os.path.join(cwd, "ground_scaled_F.mtz.txt")
ground_scaled_F_as_text = open(ground_scaled_F_as_text_filename, 'r')

ground_scaled_F_as_dict = {}
list_of_indices_ground = []

for line in ground_scaled_F_as_text.readlines()
  h = line.split(',')[0]
  k = line.split(',')[1]
  l = line.split(',')[2]
  F = line.split(',')[3]
  sigF = line.split(',')[4]
  F_sigF_list = [F, sigF]
  miller_index = h+" "+k+" "+l
  list_of_indices_ground.append(miller_index)
  reflection = {miller_index : F_sigF_list}
  ground_scaled_F_as_dict.update(reflection)


excited_scaled_F_as_text_filename = os.path.join(cwd, "excited_scaled_F.mtz.txt")
excited_scaled_F_as_text = open(excited_scaled_F_as_text_filename, 'r')

excited_scaled_F_as_dict = {}
list_of_indices_excited = []

for line in excited_scaled_F_as_text.readlines()
  h = line.split(',')[0]
  k = line.split(',')[1]
  l = line.split(',')[2]
  F = line.split(',')[3]
  sigF = line.split(',')[4]
  F_sigF_list = [F, sigF]
  miller_index = h+" "+k+" "+l
  list_of_indices_excited.append(miller_index)
  reflection = {miller_index : F_sigF_list}
  excited_scaled_F_as_dict.update(reflection)


f_model_as_text_filename = os.path.join(cwd, "f_model.mtz.txt")
f_model_SF_as_text = open(f_model_as_text_filename, 'r')

f_model_as_dict = {}

for line in excited_SF_as_text.readlines()
  h = line.split(',')[0]
  k = line.split(',')[1]
  l = line.split(',')[2]
  F = line.split(',')[3]
  sigF = line.split(',')[4]
  F_sigF_list = [F, sigF]
  miller_index = h+" "+k+" "+l
  reflection = {miller_index : F_sigF_list}
  f_model_as_dict.update(reflection)


#make list of common reflections in ground and excited states

common_reflections = []

for ground_index in list_of_indices_ground:
  if ground_index in list_of_indices_excited:
    common_reflections.append(ground_index)
#  else:
#    print ground_index

#calculate difference Fs and write to ascii .hkl file

diff_hkl_file_fullpath = os.path.join(difference_refinement_path, "difference_structure_factors.hkl")
diff_hkl_file = open(diff_hkl_file_fullpath, 'w')

for reflection in common_reflections:
  h_diff = reflection.split(' ')[0]
  k_diff = reflection.split(' ')[1]
  l_diff = reflection.split(' ')[2]
  ground_F = float(ground_scaled_F_as_dict[reflection][0])
  ground_sigF = float(ground_scaled_F_as_dict[reflection][1])
  excited_F = float(excited_scaled_F_as_dict[reflection][0])
  excited_sigF = float(excited_scaled_F_as_dict[reflection][1])
  f_model = float(f_model_as_dict[reflection][0])
  diff_F = excited_F - (ground_F - f_model)
  diff_sigF = sqrt((ground_sigF**2)+(excited_sigF**2))
  hkl_line = "%4s%4s%4s%8.2f%8.2f\n" % (h_diff, k_diff, l_diff, diff_F, diff_sigF)
  diff_hkl_file.write(hkl_line)

diff_hkl_file.close()

#convert hkl file to mtz using ccp4 f2mtz
#note that the space group is hard-coded in here. this can be improved in future

f2mtz_script_fullpath = os.path.join(difference_refinement_path, "f2mtz.sh")
f2mtz_script = open(diff_hkl_file_fullpath, 'w')

f2mtz_command = """#!/bin/sh
f2mtz hklin %s/difference_structure_factors.hkl hklout %s/difference_structure_factors.mtz <<EOF
SYMMETRY 96
LABOUT H  K  L  F  sigF
CTYPE  H  H  H  F  Q
EOF
""" % (difference_refinement_path, difference_refinement_path) 

f2mtz_script.write(f2mtz_command)
f2mtz_script.close()

os.chdir(difference_refinement_path)
os.system("./f2mtz.sh")


#add R-free flags back to the difference structure factors, again with reflection file editor

RFE2_def_file_fullpath = os.path.join(cwd, "reflection_file_editor.def")
RFE2_def_file = open(RFE2_def_file_fullpath, 'w')

input_diff_f_mtz = "%s/difference_structure_factors.mtz" % (difference_refinement_path)
hkl_file = any_reflection_file(input_diff_f_mtz)
info = hkl_file.file_content()
xtals =info.crystals()
if len(xtals) > 1:
  xtal = xtals[0]
else:
  xtal =xtals
space_group = info.space_group_number()
unit_cell = xtal.unit_cell_parameters()
uc_a, uc_b, uc_c, uc_al, uc_be, uc_ga = unit_cell

uc_str = str(unit_cell)[1:-1]
 
RFE_def_text = """show_arrays = False
dry_run = False
verbose = True
mtz_file {
  crystal_symmetry {
    unit_cell = %s
    space_group = %s
    output_unit_cell = None
    output_space_group = None
    change_of_basis = None
    eliminate_invalid_indices = False
    expand_to_p1 = False
    disable_unit_cell_check = True
    disable_space_group_check = False
    eliminate_sys_absent = True
  }
  d_max = None
  d_min = None
  wavelength = None
  output_file = %s/difference_structure_factors_flags.mtz
  job_title = None
  resolve_label_conflicts = True
  exclude_reflection = None
  miller_array {
    file_name = %s/difference_structure_factors.mtz
    labels = IOBS,SIGIOBS
    output_labels = IOBS SIGIOBS
    column_root_label = None
    d_min = None
    d_max = None
    output_as = *auto intensities amplitudes_fw amplitudes
    anomalous_data = *Auto merged anomalous
    force_type = *auto amplitudes intensities
    scale_max = None
    scale_factor = None
    remove_negatives = False
    massage_intensities = False
    filter_by_signal_to_noise = None
    add_b_iso = None
    add_b_aniso = 0 0 0 0 0 0
    shuffle_values = False
    reset_values_to = None
  }
  miller_array {
    file_name = %s
    labels = R-free-flags
    output_labels = FreeR_flag
    column_root_label = None
    d_min = None
    d_max = None
    output_as = *auto intensities amplitudes_fw amplitudes
    anomalous_data = *Auto merged anomalous
    force_type = *auto amplitudes intensities
    scale_max = None
    scale_factor = None
    remove_negatives = False
    massage_intensities = False
    filter_by_signal_to_noise = None
    add_b_iso = None
    add_b_aniso = 0 0 0 0 0 0
    shuffle_values = False
    reset_values_to = None
  }
  r_free_flags {
    generate = True
    force_generate = False
    new_label = "FreeR_flag"
    fraction = 0.1
    max_free = 2000
    lattice_symmetry_max_delta = 5
    use_lattice_symmetry = True
    use_dataman_shells = False
    n_shells = 20
    random_seed = None
    extend = True
    old_test_flag_value = None
    export_for_ccp4 = False
    preserve_input_values = True
    warn_if_all_same_value = True
    adjust_fraction = False
    d_eps = 0.0001
    relative_to_complete_set = False
    remediate_mismatches = False
  }
}

""" % (uc_str, space_group, difference_refinement_path, difference_refinement_path, ground_data_fullpath)
  
RFE2_def_file.write(RFE2_def_text)
RFE2_def_file.close()

print "wrote RFE2 def file successfully!"

subprocess.call(['iotbx.reflection_file_editor', 'reflection_file_editor.def'])


#perform difference refinement with phenix.refine (7 macro cycles)
#I wonder what the optimal parameters are for this refinement. We should find out and make this step more rigorous by creating a full .def file
#play with weights

difference_refinement_command = "phenix.refine %s ./difference_structure_factors_flags.mtz refinement.main.number_of_macrocycles=7"
os.system(difference_refinement_command)


#return to original working directory and clean up by moving .txt files to new subdirectory

os.chdir(..)

os.system("mkdir ./mtz_as_text")
os.system("mv ./*.mtz.txt ./mtz_as_text/")