import os
import sys
import math
import string
from iotbx.reflection_file_reader import any_reflection_file
import subprocess
import numpy as np



def specify_input_files( starting_directory ):

  print "Collecting input from command line arguments..."

  ground_coords_fullpath = os.path.join(starting_directory, sys.argv[1])
  ground_data_fullpath = os.path.join(starting_directory, sys.argv[2])
  excited_data_fullpath = os.path.join(starting_directory, sys.argv[3])

  return ground_coords_fullpath, ground_data_fullpath, excited_data_fullpath

  print "Input collected!"


def prepare_directories( starting_directory ):

  print "Preparing directories..."

  scaling_directory_path = os.path.join(starting_directory, "scaling")
  os.mkdir(scaling_directory_path, 0755)
  extrapolated_SF_path = os.path.join(starting_directory, "extrapolated_SF")
  os.mkdir(extrapolated_SF_path, 0755)

  return scaling_directory_path, extrapolated_SF_path

  print "Preparation complete!"


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
  scale_to_I_model_command = """phenix.scale_and_merge reference_data=%s data=%s data_labels=Iobs,SIGIobs anomalous=False """ % (I_model_mtz_filename, unscaled_intensity_mtz_filename)
  os.system(scale_to_I_model_command)
  scaled_intensity_mtz_filename = scaled_output_filename_prefix+"_scaled.mtz"
  rename_command = "mv ./scaled_data.mtz ./%s" % (scaled_intensity_mtz_filename)
  os.system(rename_command)

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
    labels = I,SIGI
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
    labels = I,SIGI
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


def calculate_extrapolated_F_with_weights_as_dict( ground_state_Fobs_as_text_fullpath, excited_state_Fobs_as_text_fullpath, extrapolation_factor ):
  
  ground_state_F_as_dict = store_Fhkl_from_text_as_dict(ground_state_F_as_text_fullpath)
  excited_state_F_as_dict = store_Fhkl_from_text_as_dict(excited_state_F_as_text_fullpath)
  common_reflections = identify_common_reflections(ground_state_F_as_dict, excited_state_F_as_dict)
  delta_F_dict = {}
  sum_delta_F = 0.0
  sum_sig_delta_F = 0.0
  n_delta_F = 0.0
  extrapolated_F_dict = {}

  for miller_index in common_reflections:
    fobs_gr = float(ground_state_F_as_dict.get(miller_index)[0])
    sigfobs_gr = float(ground_state_F_as_dict.get(miller_index)[1])
    fobs_ex = float(excited_state_F_as_dict.get(miller_index)[0])
    sigfobs_ex = float(excited_state_F_as_dict.get(miller_index)[1])
    f = []
    delta_F = fobs_ex - fobs_gr
    f.append(delta_F)
    abs_delta_F = abs(delta_F)
    sum_abs_delta_F = sum_abs_delta_F + abs_delta_F
    sig_delta_F = math.sqrt((sigfobs_ex**2)+(sigfobs_gr**2))
    f.append(sig_delta_F)
    sum_sig_delta_F = sum_sig_delta_F + sig_delta_F
    values = {miller_index : f}
    delta_F_dict.update(values)
    n_delta_F = n_delta_F + 1.0

  avg_abs_delta_F = sum_abs_delta_F / n_delta_F
  avg_sig_delta_F = sum_sig_delta_F / n_delta_F

  for miller_index in delta_F_dict.keys():
    fobs_gr = float(ground_state_F_as_dict.get(miller_index)[0])
    sigfobs_gr = float(ground_state_F_as_dict.get(miller_index)[1])
    delta_F = float(delta_F_dict.get(miller_index)[0])
    abs_delta_F = abs(delta_F)
    sig_delta_F = float(delta_F_dict.get(miller_index)[1])
    weight = 1.0 / (1+ ((sig_delta_F**2)/(avg_sig_delta_F**2)) + (0.05*((abs_delta_F**2)/(avg_abs_delta_F**2))))
    weighted_delta_F = weight * delta_F
    N = float(extrapolation_factor)
    e = []
    extrap_F = fobs_gr + (N * weighted_delta_F)
    e.append(extrap_F)
    sig_extrap_F = math.sqrt(((N**2)*(weight**2)*(sig_delta_F**2))+(sigfobs_gr**2))
    e.append(sig_extrap_F)
    values = {miller_index : e}
    extrapolated_F_dict.update(values)

  return extrapolated_F_dict


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


def add_R_free_flags_to_mtz( data_mtz_filename, R_free_reference_mtz_fullpath ):

  working_directory = os.getcwd()
  RFE_def_file_fullpath = os.path.join(working_directory, "reflection_file_editor_flags.def")
  RFE_def_file = open(RFE_def_file_fullpath, 'w')

  output_mtz_prefix = data_mtz_filename.split('.')[-3]
  output_mtz_filename = "%s_flags.mtz" % (output_mtz_prefix)

  input_diff_f_mtz = "%s/%s" % (working_directory, data_mtz_filename)
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
  output_file = %s/%s
  job_title = None
  resolve_label_conflicts = True
  exclude_reflection = None
  miller_array {
    file_name = %s/%s
    labels = FOBS,SIGFOBS
    output_labels = FOBS SIGFOBS
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
    labels = FreeR_flag
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

""" % (uc_str, space_group, working_directory, output_mtz_filename, working_directory, data_mtz_filename, R_free_reference_mtz_fullpath)
  
  RFE_def_file.write(RFE_def_text)
  RFE_def_file.close()

  print "wrote RFE def file successfully!"

  subprocess.call(['iotbx.reflection_file_editor', 'reflection_file_editor_flags.def'])


def run_simple_phenix_refinement( model, data ):

  refinement_command = "phenix.refine %s %s strategy=individual_sites+individual_adp refinement.main.number_of_macro_cycles=7 refinement.main.nproc=8 refinement.modify_start_model.modify.adp.convert_to_isotropic=True refinement.modify_start_model.modify.adp.set_b_iso=20.00 refinement.modify_start_model.modify.sites.shake=0.5 refinement.target_weights.optimize_xyz_weight=True refinement.target_weights.optimize_adp_weight=True" % (model, data) # add weight opt??
  os.system(refinement_command)


def main():

  ground_state_pdb_filename = sys.argv[1]
  ground_state_mtz_filename = sys.argv[2]
  excited_state_mtz_filename = sys.argv[3]
  R_free_reference_mtz_filename = sys.argv[4]

  extrapolation_factors = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0]

  starting_directory = os.getcwd()
  ground_coords_fullpath, ground_data_fullpath, excited_data_fullpath = specify_input_files(starting_directory)
  scaling_directory_path, extrapolated_SF_path = prepare_directories(starting_directory)

  R_free_reference_mtz_fullpath = os.path.join(starting_directory, R_free_reference_mtz_filename)

  for data_file in [ground_state_mtz_filename, excited_state_mtz_filename]:
    move_command = "cp ./%s %s" % (data_file, scaling_directory_path)
    os.system(move_command)

  coords_fullpath = os.path.join(starting_directory, ground_state_pdb_filename)
  I_model_mtz_filename = generate_I_model(coords_fullpath, scaling_directory_path)

  os.chdir(scaling_directory_path)

  gr_scaled_intensity_mtz_filename = scale_intensities_to_I_model(ground_state_mtz_filename, I_model_mtz_filename)
  convert_scaled_intensities_to_structure_factors_FandW(gr_scaled_intensity_mtz_filename)
  scaled_ground_F_filename = gr_scaled_intensity_mtz_filename.split('.')[-2]+"_asF.mtz"
  scaled_ground_F_as_text_filename = convert_mtz_to_text(scaled_ground_F_filename)
  scaled_ground_F_as_text_fullpath = os.path.join(scaling_directory_path, scaled_ground_F_as_text_filename)


  ex_scaled_intensity_mtz_filename = scale_intensities_to_I_model(excited_state_mtz_filename, I_model_mtz_filename)  
  convert_scaled_intensities_to_structure_factors_FandW(ex_scaled_intensity_mtz_filename)
  scaled_excited_F_filename = ex_scaled_intensity_mtz_filename.split('.')[-2]+"_asF.mtz"
  scaled_excited_F_as_text_filename = convert_mtz_to_text(scaled_excited_F_filename)
  scaled_excited_F_as_text_fullpath = os.path.join(scaling_directory_path, scaled_excited_F_as_text_filename)

  os.chdir(extrapolated_SF_path)
   
  for extrapolation_factor in extrapolation_factors:
    ex_fac_str = str(extrapolation_factor)
    ex_fac_dir_name = "extrapolated_SF"+ex_fac_str
    ex_fac_file_prefix = ex_fac_dir_name
    ex_fac_dir_path = os.path.join(extrapolated_SF_path, ex_fac_dir_name)
    os.mkdir(ex_fac_dir_path, 0755)
    os.chdir(ex_fac_dir_path)

    extrapolated_F_dict = calculate_extrapolated_F_with_weights_as_dict(scaled_ground_F_as_text_fullpath, scaled_excited_F_as_text_fullpath, extrapolation_factor)

    extrap_hkl_filename = create_hkl_file(extrapolated_F_dict, ex_fac_dir_path, ex_fac_file_prefix)
    extrap_mtz_filename = convert_hkl_to_mtz(extrap_hkl_filename, ex_fac_file_prefix)

    add_R_free_flags_to_mtz(extrap_mtz_filename, R_free_reference_mtz_fullpath)

    get_pdb_file = "cp %s/%s ." % (starting_directory, ground_state_pdb_filename)
    os.system(get_pdb_file)

    data_filename_for_refinement = ex_fac_file_prefix.split('.')[-2]+"_flags.mtz"
    run_simple_phenix_refinement(ground_state_pdb_filename, data_filename_for_refinement)

    os.chdir(extrapolated_SF_path)

main()
