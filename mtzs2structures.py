#Script for copying a set of Rfree flags across multiple isomorphous data sets and
#perform rigid-body refinement to initiate rebuild.
#Usage: python mtzs2structures.py [isomorphous.pdb] [r_free.mtz]
#Run from a directory containing all of the derivative MTZ files

import os
import sys
import math
import string

#Take command line input (isomorphous-ish PDB file and an MTZ file with desired FreeR_flags)
#assess MTZ files
#make a directory for each MTZ

isomorphous_input_pdb = sys.argv[1]
isomorphous_input_pdb_fullpath = os.path.join(cwd, isomorphous_input_pdb)
r_free_mtz = sys.argv[2]
r_free_mtz_fullpath = os.path.join(cwd, r_free_mtz)

root_dir = os.getcwd()

original_mtz_files = []
original_mtz_filenames = []
dir_paths = []

for root, dirs, filenames in os.walk('./'):
  for filename in filenames:
    if filename.split('.')[-1] == 'mtz':
      fullpath = os.path.join(cwd, filename)
      original_mtz_files.append(fullpath)
      original_mtz_filenames.append(filename.split('.')[-2])

for i, mtz_filename in enumerate(original_mtz_filenames):
  dir_path = os.path.join(root_dir, mtz_filename)
  dir_paths.append(dir_path)
  os.mkdir(dir_path, 0755)
  

#Run phenix.xtriage
#retrieve unit cell parameters and space group
#create input for reflection_file_editor and run
#create new pdb with modified CRYST1
#set up phenix.refine for rigid body refinement and run

for n, path in enumerate(dir_paths):
  
  os.chdir(path)

  cwd = os.getcwd()

  xtriage_cmd = 'phenix.xtriage ./'+str(original_mtz_filenames[n])+'.mtz'
  os.system(xtriage_cmd)

  xtriage_outfile = open('./logfile.log', 'r')
  unit_cell_line = xtriage_outfile.readlines()[20]
  unit_cell = unit_cell_line.split('=')[-1]
  uc_a = unit_cell_line.split(' ')[-6]
  uc_b = unit_cell_line.split(' ')[-5]
  uc_c = unit_cell_line.split(' ')[-4]
  uc_al = unit_cell_line.split(' ')[-3]
  uc_be = unit_cell_line.split(' ')[-2]
  uc_ga = unit_cell_line.split(' ')[-1]

  space_group_line = xtriage_outfile.readlines()[21]
  space_group = space_group_line.split('"')[-2]

  input_refl_mtz_file = original_mtz_files[n]
  output_mtz_filename = original_mtz_filenames[n]+"_flags.mtz"
  output_mtz_file = os.path.join(cwd, flags_mtz_filename)

  RFE_def_file_fullpath = os.path.join(cwd, "reflection_file_editor.def")
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
    labels = IMEAN,SIGIMEAN
    output_labels = IMEAN SIGIMEAN
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

""" % (unit_cell, space_group, output_mtz_file, input_refl_mtz_file, r_free_mtz_fullpath)
  
  RFE_def_file.write(RFE_def_text)
  RFE_def_file.close()

  os.system('iotbx.reflection_file_editor reflection_file_editor.def')

  prepare input pdb file
      pdbtools
      cryst1

  cryst1_sg = space_group.rjust(11)
  cryst1_line = "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%s    \n"

  orig_pdb_file = open(isomorphous_input_pdb_fullpath, 'r')
  output_pdb_file_fullpath = os.path.join(cwd, "coords_in_current_cell.pdb")
  output_pdb_file = open(output_pdb_file_fullpath, 'w')

  output_pdb_lines = []

  for line in orig_pdb_file.readlines():
    if line.startswith('ATOM'):
      output_pdb_lines.append(line)

  output_pdb_file.write(cryst1_line)

  for atom in output_pdb_lines:
    output_pdb_file.write(atom)

  output_pdb_file.close()

  phenix_refine_command = "phenix.refine %s %s strategy=rigid_body" % (output_mtz_file, output_pdb_file)
  os.system(phenix_refine_command)

  os.chdir('..')


