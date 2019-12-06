# needs to run with libtbx.python so source your phenix install before running
# Use a recent nightly build of phenix for functional phenix.scale_and_merge
# Usage: libtbx.python weighted_difference_map.py ground_state.mtz excited_state.mtz scaling_reference.pdb

import os
import sys
import math
import string
from iotbx.reflection_file_reader import any_reflection_file
import subprocess
import numpy as np 
import operator


def generate_I_model( coords_fullpath, I_model_output_dir_fullpath ):

  print "Using phenix.fmodel to generate Icalc for input model: %s" % (coords_fullpath)

  coords_filename = coords_fullpath.split('/')[-1]
  I_model_mtz_filename_prefix = coords_filename.split('.')[-2]
  I_model_mtz_filename = I_model_mtz_filename_prefix+"_I_model.mtz"
  generate_I_model_command = """phenix.fmodel %s output.type=real output.obs_type=intensities output.label=I high_resolution=1.5 add_sigmas=True scale=0.1 output.file_name="%s_I_model.mtz" """ % (coords_fullpath, I_model_mtz_filename_prefix)
  os.system(generate_I_model_command)
  relocate_I_model_command = "mv ./%s_I_model.mtz %s" % (I_model_mtz_filename_prefix, I_model_output_dir_fullpath)
  os.system(relocate_I_model_command)

  return I_model_mtz_filename

  print "An MTZ file containing Imodel has been created: %s/%s_I_model.mtz" % (I_model_output_dir_fullpath, I_model_mtz_filename_prefix)

 
def scale_intensities_to_I_model( unscaled_intensity_mtz_filename, I_model_mtz_filename ):

  start_message = """Scaling experimental MTZ file: %s

to MTZ containing Imodel: %s

Using phenix.scale_and_merge

""" % (unscaled_intensity_mtz_filename, I_model_mtz_filename)
  print start_message

  scaled_output_filename_prefix = unscaled_intensity_mtz_filename.split('.')[-2]
  scale_to_I_model_command = """phenix.scale_and_merge reference_data=%s data=%s anomalous=False output.file_name="%s_scaled.mtz" """ % (I_model_mtz_filename, unscaled_intensity_mtz_filename, scaled_output_filename_prefix)
  os.system(scale_to_I_model_command)
  scaled_intensity_mtz_filename = scaled_output_filename_prefix+"_scaled.mtz"

  return scaled_intensity_mtz_filename

  print "Scaling is complete!!"


def convert_scaled_intensities_to_structure_factors( scaled_intensity_mtz_filename ):

  cwd = os.getcwd()
  scaled_intensity_mtz_fullpath = os.path.join(cwd, scaled_intensity_mtz_filename)
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
  output_file = %s/%s_asF.mtz
  job_title = None
  resolve_label_conflicts = True
  exclude_reflection = None
  miller_array {
    file_name = %s
    labels = Iobs,SIGIobs
    output_labels = FOBS SIGFOBS
    column_root_label = None
    d_min = None
    d_max = None
    output_as = auto intensities amplitudes_fw *amplitudes
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

""" % (uc_str, space_group, cwd, scaled_structure_factor_mtz_prefix, scaled_intensity_mtz_fullpath)
  
  RFE_def_file.write(RFE_def_text)
  RFE_def_file.close()

  print "Wrote RFE def file successfully! - Now converting intensities to structure factor amplitudes with RFE."

  subprocess.call(['iotbx.reflection_file_editor', 'reflection_file_editor.def'])

  output_mtz_filename = "%s_asF.mtz" % (scaled_structure_factor_mtz_prefix)

  return output_mtz_filename

  print "Conversion to Structure Factors is finished!"


def convert_scaled_intensities_to_structure_factors_FandW( scaled_intensity_mtz_filename ):

  cwd = os.getcwd()
  scaled_intensity_mtz_fullpath = os.path.join(cwd, scaled_intensity_mtz_filename)
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
  output_file = %s/%s_asF.mtz
  job_title = None
  resolve_label_conflicts = True
  exclude_reflection = None
  miller_array {
    file_name = %s
    labels = Iobs,SIGIobs
    output_labels = FOBS SIGFOBS
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

""" % (uc_str, space_group, cwd, scaled_structure_factor_mtz_prefix, scaled_intensity_mtz_fullpath)
  
  RFE_def_file.write(RFE_def_text)
  RFE_def_file.close()

  print "Wrote RFE def file successfully! - Now converting intensities to structure factor amplitudes with RFE."

  subprocess.call(['iotbx.reflection_file_editor', 'reflection_file_editor.def'])

  output_mtz_filename = "%s_asF.mtz" % (scaled_structure_factor_mtz_prefix)
  
  return output_mtz_filename

  print "Conversion to Structure Factors is finished!"


def convert_mtz_to_text( structure_factor_mtz_filename ):

  print "Writing reflection data in %s as ASCII text..." % (structure_factor_mtz_filename)

  structure_factor_mtz_prefix = structure_factor_mtz_filename.split('.')[-2]

  mtz_dump_command = "phenix.mtz.dump -c -f s ./%s > %s_mtz.txt" % (structure_factor_mtz_filename, structure_factor_mtz_prefix)
  os.system(mtz_dump_command)

  output_filename = "%s_mtz.txt" % (structure_factor_mtz_prefix)
  return output_filename

  print "The data have been converted from MTZ to text."


def store_Fhkl_from_text_as_dict( input_reflections_as_text_fullpath ):

  reflections_as_text = open(input_reflections_as_text_fullpath, 'r')
  reflections_as_dict = {}
  for i, line in enumerate(reflections_as_text.readlines()):
    if i > 0:
      h = line.split(',')[0]
      k = line.split(',')[1]
      l = line.split(',')[2]
      F = line.split(',')[3]
      sigF = line.split(',')[4]
      F_sigF_list = [F, sigF]
      miller_index = h+" "+k+" "+l
      reflection = {miller_index : F_sigF_list}
      if F != '':
        reflections_as_dict.update(reflection)

  return reflections_as_dict

def identify_common_reflections( reflection_dict_1, reflection_dict_2 ):

  # input parameters are dictionaries...
  common_reflections = []
  for reflection in reflection_dict_1.keys():
    if reflection in reflection_dict_2.keys():
      common_reflections.append(reflection)

  n_common_ref = len(common_reflections)

  print "There are %i reflections common to both data sets" % (n_common_ref)

  return common_reflections


def calculate_weighted_delta_F_as_dict( ground_state_F_as_text_fullpath, excited_state_F_as_text_fullpath ):

  ground_state_F_as_dict = store_Fhkl_from_text_as_dict(ground_state_F_as_text_fullpath)
  excited_state_F_as_dict = store_Fhkl_from_text_as_dict(excited_state_F_as_text_fullpath)
  common_reflections = identify_common_reflections(ground_state_F_as_dict, excited_state_F_as_dict)
  Fo_minus_Fo_and_sigma_dict = {}
  sum_Fo_minus_Fo = 0.0
  sum_sig_Fo_minus_Fo = 0.0
  n_Fo_minus_Fo = 0.0

  for miller_index in common_reflections:
    fobs_1 = float(ground_state_F_as_dict.get(miller_index)[0])
    sigfobs_1 = float(ground_state_F_as_dict.get(miller_index)[1])
    fobs_2 = float(excited_state_F_as_dict.get(miller_index)[0])
    sigfobs_2 = float(excited_state_F_as_dict.get(miller_index)[1])
    f = []
    Fo_minus_Fo = fobs_2 - fobs_1
    f.append(Fo_minus_Fo)
    sum_Fo_minus_Fo = sum_Fo_minus_Fo + Fo_minus_Fo
    sig_Fo_minus_Fo = sqrt((sigfobs_2**2)+(sigfobs_1**2))
    sum_sig_Fo_minus_Fo = sum_sig_Fo_minus_Fo + sig_Fo_minus_Fo
    f.append(sig_Fo_minus_Fo)
    values = {miller_index : f}
    Fo_minus_Fo_and_sigma_dict.update(values)
    n_Fo_minus_Fo = n_Fo_minus_Fo + 1.0

  avg_delta_F = sum_Fo_minus_Fo / n_Fo_minus_Fo
  avg_sig_delta_F = sum_sig_Fo_minus_Fo / n_Fo_minus_Fo

  weighted_delta_F_as_dict = {}

  for miller_index in Fo_minus_Fo_and_sigma_dict.keys():
    delta_F = float(Fo_minus_Fo_and_sigma_dict.get(miller_index)[0])
    sig_delta_F = float(Fo_minus_Fo_and_sigma_dict.get(miller_index)[1])
    weight = 1.0 / (1+ ((sig_delta_F**2)/(avg_sig_delta_F**2)) + (0.05*((delta_F**2)/(avg_delta_F**2))))
    weighted_delta_F = weight * delta_F
    values = {miller_index : weighted_delta_F}
    weighted_delta_F_as_dict.update(values)


  return weighted_delta_F_as_dict


def create_hkl_file( F_dict, hkl_file_path, output_prefix ):

  hkl_filename = "%s.hkl" % (output_prefix)
  hkl_file_fullpath = os.path.join(hkl_file_path, hkl_filename)
  hkl_file = open(hkl_file_fullpath, 'w')

  for reflection in F_dict.keys():
    h = reflection.split(' ')[0]
    k = reflection.split(' ')[1]
    l = reflection.split(' ')[2]
    F_hkl = F_dict.get(reflection)[0]
    sigF_hkl = F_dict.get(reflection)[1]
    hkl_line = "%4s%4s%4s%8.2f%8.2f\n" % (h, k, l, F_hkl, sigF_hkl)
    hkl_file.write(hkl_line)

  hkl_file.close()

  return hkl_filename


def convert_hkl_to_mtz( hkl_filename, output_prefix ):

  working_directory = os.getcwd()
  f2mtz_script_fullpath = os.path.join(working_directory, "f2mtz.sh")
  f2mtz_script = open(f2mtz_script_fullpath, 'w')

  output_mtz_filename = "%s.mtz" % (output_prefix)

  f2mtz_command = """#!/bin/sh
  f2mtz hklin %s/%s hklout %s/%s <<EOF
  SYMMETRY 96
  CELL 79.523  79.523  38.2332 
  LABOUT H  K  L  FOBS  SIGFOBS
  CTYPE  H  H  H  F  Q
  END
  EOF
  """ % (working_directory, hkl_filename, working_directory, output_mtz_filename) 

# CELL SHOULD NOT BE HARD CODED IN FINAL VERSION

  f2mtz_script.write(f2mtz_command)
  f2mtz_script.close()
  os.chmod(f2mtz_script_fullpath, 0755)
  os.system("./f2mtz.sh")

  return output_mtz_filename


def make_maps( mtz_filename, pdb_filename, output_prefix ):

  working_directory = os.getcwd()

  pdb_prefix = pdb_filename.split('.')[-2]
  mtz_prefix = mtz_filename.split('.')[-2]
  recip_space_arrays_output_filename = pdb_prefix+"_"+mtz_prefix+".mtz"

  recip_space_arrays_command = "phenix.reciprocal_space_arrays %s %s" % (pdb_filename, mtz_filename)
  os.system(recip_space_arrays_command)

  fft_script_fullpath = os.path.join(working_directory, "fft.sh")
  fft_script = open(fft_script_fullpath, 'w')

  output_map_filename = "%s.ccp4" % (output_prefix)

  fft_command = """#!/bin/sh
  fft hklin %s/%s mapout %s/%s <<eof-fft
  LABIN F1=FOBS PHI=PHIMODEL
  END
  eof-fft
  """ % (working_directory, recip_space_arrays_output_filename, working_directory, output_map_filename) 

  fft_script.write(fft_command)
  fft_script.close()
  os.system("chmod ./fft.sh 755")
  os.system("./fft.sh")


def main():

  output_prefix = "weighted_delta_F"

  starting_dir = os.getcwd()
  ground_state_data_input_filename = sys.argv[1]
  ground_state_data_fullpath = os.path.join(starting_dir, ground_state_data_input_filename)
  excited_state_data_input_filename = sys.argv[2]
  excited_state_data_fullpath = os.path.join(starting_dir, excited_state_data_input_filename)
  pdb_scaling_reference_input_filename = sys.argv[3]
  pdb_scaling_reference_fullpath = os.path.join(starting_dir, pdb_scaling_reference_input_filename)
  pdb_scaling_reference_prefix = pdb_scaling_reference_input_filename.split('.')[-2]
  temp_dir_path = os.path.join(starting_dir, "temp")
  os.mkdir(temp_dir_path, 0755)

  I_model_mtz_filename = generate_I_model(pdb_scaling_reference_fullpath, temp_dir_path)

  os.chdir(temp_dir_path)

  reflection_text_files = []

  for filename in [ground_state_data_input_filename, excited_state_data_input_filename]:
    copy_unscaled_mtz_to_subdir = "cp ../%s %s" % (filename, temp_dir_path)
    os.system(copy_unscaled_mtz_to_subdir)
    scaled_I_mtz = scale_intensities_to_I_model(filename, I_model_mtz_filename)
    scaled_F_mtz = convert_scaled_intensities_to_structure_factors_FandW(scaled_I_mtz)
    reflection_text_file = convert_mtz_to_text(scaled_F_mtz)
    reflection_text_fullpath = os.path.join(temp_dir_path, reflection_text_file)
    reflection_text_files.append(reflection_text_fullpath)

  os.chdir(starting_dir)

  weighted_delta_F_dict = calculate_weighted_delta_F_as_dict(reflection_text_files[1], reflection_text_files[0])

  hkl_filename = create_hkl_file(weighted_delta_F_dict, starting_dir, output_prefix) #Need Args
  output_mtz_filename = convert_hkl_to_mtz(hkl_filename, output_prefix)
  make_maps(output_mtz_filename, pdb_scaling_reference_input_filename, output_prefix)

  os.system("rm -r ./temp")
