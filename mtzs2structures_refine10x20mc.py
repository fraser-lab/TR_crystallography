#Script for copying a set of Rfree flags across multiple isomorphous data sets and
#perform rigid-body refinement to initiate rebuild.
#Usage: libtbx.python mtzs2structures.py [isomorphous.pdb] [r_free.mtz]
#Run from a directory containing all of the derivative MTZ files

import os
import sys
import math
import string
from iotbx.reflection_file_reader import any_reflection_file
import subprocess
from random import randint

os.system("source /netapp/home/mthomp/phenix-installer-1.13-2998-intel-linux-2.6-x86_64-centos5/phenix_env.sh")


#Take command line input (isomorphous-ish PDB file and an MTZ file with desired FreeR_flags)
#assess MTZ files
#make a directory for each MTZ
starting_directory = os.getcwd()
isomorphous_input_pdb = sys.argv[1]
isomorphous_input_pdb_fullpath = os.path.join(starting_directory, isomorphous_input_pdb)
r_free_mtz = sys.argv[2]
r_free_mtz_fullpath = os.path.join(starting_directory, r_free_mtz)

original_mtz_files = []
original_mtz_filenames = []
dir_paths = []

for root, dirs, filenames in os.walk('./'):
  for filename in filenames:
    if filename != r_free_mtz:
      if filename.split('.')[-1] == 'mtz':
        fullpath = os.path.join(starting_directory, filename)
        original_mtz_files.append(fullpath)
        original_mtz_filenames.append(filename.split('.')[-2])
    else:
      pass

for i, mtz_filename in enumerate(original_mtz_filenames):
  dir_path = os.path.join(starting_directory, mtz_filename)
  dir_paths.append(dir_path)
  os.mkdir(dir_path, 0755)
  

#retrieve unit cell parameters and space group
#create input for reflection_file_editor and run
#create new pdb with modified CRYST1
#set up phenix.refine for rigid body refinement and run

for n, path in enumerate(dir_paths):
  print path
  os.chdir(path)

  # xtriage_cmd = 'phenix.xtriage ../'+str(original_mtz_filenames[n])+'.mtz'
  # os.system(xtriage_cmd)
  hkl_file = any_reflection_file("../"+str(original_mtz_filenames[n])+'.mtz')
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

  # xtriage_outfile = open('./logfile.log', 'r')
  # # unit_cell_line = xtriage_outfile.readlines()[20]
  # for line in xtriage_outfile.readlines():
  #   print line
  # unit_cell = unit_cell_line.split('=')[-1]
  # uc_a = unit_cell_line.split(' ')[-6]
  # uc_b = unit_cell_line.split(' ')[-5]
  # uc_c = unit_cell_line.split(' ')[-4]
  # uc_al = unit_cell_line.split(' ')[-3]
  # uc_be = unit_cell_line.split(' ')[-2]
  # uc_ga = unit_cell_line.split(' ')[-1]

  # space_group_line = xtriage_outfile.readlines()[21]
  # space_group = space_group_line.split('"')[-2]


  input_refl_mtz_file = original_mtz_files[n]
  output_mtz_filename = original_mtz_filenames[n]+"_flags.mtz"
  # output_mtz_file = os.path.join(cwd, flags_mtz_filename)
  output_mtz_file = os.path.join(path, output_mtz_filename)

  RFE_def_file_fullpath = os.path.join(path, "reflection_file_editor.def")
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

""" % (uc_str, space_group, output_mtz_file, input_refl_mtz_file, r_free_mtz_fullpath)
  
  RFE_def_file.write(RFE_def_text)
  RFE_def_file.close()

  print "wrote RFE def file successfully!"

  subprocess.call(['iotbx.reflection_file_editor', 'reflection_file_editor.def'])

  # prepare input pdb file
  #     pdbtools
  #     cryst1

  cryst1_sg = str(space_group).rjust(10)
  cryst1_line = "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%s    \n" % (uc_a, uc_b, uc_c, uc_al, uc_be, uc_ga, space_group)
  # cryst1_line = "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%s    \n" % (unit_cell, space_group)
  # cryst1_line = "CRYST1%s%s    \n" % (unit_cell,space_group)

  orig_pdb_file = open(isomorphous_input_pdb_fullpath, 'r')
  output_pdb_file_fullpath = os.path.join(path, "coords_in_current_cell.pdb")
  output_pdb_file = open(output_pdb_file_fullpath, 'w')

  output_pdb_lines = []

  for line in orig_pdb_file.readlines():
    if line.startswith('ATOM'):
      output_pdb_lines.append(line)

  output_pdb_file.write(cryst1_line)

  for atom in output_pdb_lines:
    output_pdb_file.write(atom)

  output_pdb_file.close()

  pdbtools_cmd = "phenix.pdbtools %s convert_to_isotropic=True set_b_iso=15.00 remove_alt_confs=True remove='not protein' output.file_name='coords_in_current_cell.pdb'" % (output_pdb_file_fullpath)
  os.system(pdbtools_cmd)

  phenix_refine_rigid_body_cmd = "phenix.refine %s coords_in_current_cell.pdb strategy=rigid_body --space_group %s --unit_cell %s refinement.main.number_of_macro_cycles=3 output.prefix=%s_rigid_body_refine" % (output_mtz_file, space_group, uc_str.replace(" ",""), original_mtz_filenames[n])
  os.system(phenix_refine_rigid_body_cmd)


  rb_output_pdb_filename = '%s_rigid_body_refine_001.pdb' % (original_mtz_filenames[n])

  run1_path = os.path.join(path, "run1")
  os.mkdir(run1_path, 0755)
  run2_path = os.path.join(path, "run2")
  os.mkdir(run2_path, 0755)
  run3_path = os.path.join(path, "run3")
  os.mkdir(run3_path, 0755)
  run4_path = os.path.join(path, "run4")
  os.mkdir(run4_path, 0755)
  run5_path = os.path.join(path, "run5")
  os.mkdir(run5_path, 0755)
  run6_path = os.path.join(path, "run6")
  os.mkdir(run6_path, 0755)
  run7_path = os.path.join(path, "run7")
  os.mkdir(run7_path, 0755)
  run8_path = os.path.join(path, "run8")
  os.mkdir(run8_path, 0755)
  run9_path = os.path.join(path, "run9")
  os.mkdir(run9_path, 0755)
  run10_path = os.path.join(path, "run10")
  os.mkdir(run10_path, 0755)

  
  copy_pdb_1 = "cp ./%s ./run1/" % (rb_output_pdb_filename)
  os.system(copy_pdb_1)
  copy_pdb_2 = "cp ./%s ./run2/" % (rb_output_pdb_filename)
  os.system(copy_pdb_2)
  copy_pdb_3 = "cp ./%s ./run3/" % (rb_output_pdb_filename)
  os.system(copy_pdb_3)
  copy_pdb_4 = "cp ./%s ./run4/" % (rb_output_pdb_filename)
  os.system(copy_pdb_4)
  copy_pdb_5 = "cp ./%s ./run5/" % (rb_output_pdb_filename)
  os.system(copy_pdb_5)
  copy_pdb_6 = "cp ./%s ./run6/" % (rb_output_pdb_filename)
  os.system(copy_pdb_6)
  copy_pdb_7 = "cp ./%s ./run7/" % (rb_output_pdb_filename)
  os.system(copy_pdb_7)
  copy_pdb_8 = "cp ./%s ./run8/" % (rb_output_pdb_filename)
  os.system(copy_pdb_8)
  copy_pdb_9 = "cp ./%s ./run9/" % (rb_output_pdb_filename)
  os.system(copy_pdb_9)
  copy_pdb_10 = "cp ./%s ./run10/" % (rb_output_pdb_filename)
  os.system(copy_pdb_10)

  copy_mtz_1 = "cp ./%s ./run1/" % (output_mtz_filename)
  os.system(copy_mtz_1)
  copy_mtz_2 = "cp ./%s ./run2/" % (output_mtz_filename)
  os.system(copy_mtz_2)
  copy_mtz_3 = "cp ./%s ./run3/" % (output_mtz_filename)
  os.system(copy_mtz_3)
  copy_mtz_4 = "cp ./%s ./run4/" % (output_mtz_filename)
  os.system(copy_mtz_4)
  copy_mtz_5 = "cp ./%s ./run5/" % (output_mtz_filename)
  os.system(copy_mtz_5)
  copy_mtz_6 = "cp ./%s ./run6/" % (output_mtz_filename)
  os.system(copy_mtz_6)
  copy_mtz_7 = "cp ./%s ./run7/" % (output_mtz_filename)
  os.system(copy_mtz_7)
  copy_mtz_8 = "cp ./%s ./run8/" % (output_mtz_filename)
  os.system(copy_mtz_8)
  copy_mtz_9 = "cp ./%s ./run9/" % (output_mtz_filename)
  os.system(copy_mtz_9)
  copy_mtz_10 = "cp ./%s ./run10/" % (output_mtz_filename)
  os.system(copy_mtz_10)

  os.chdir(run1_path)
  qsub_script_fullpath = os.path.join(run1_path, "refinement.sub")
  qsub_script = open(qsub_script_fullpath, 'w')
  rand_seed_1 = randint(1000000, 9999999)
  phenix_refine_qsub_text_1 = """#$ -S /bin/bash
#$ -o %s/out
#$ -e %s/err
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=8G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=1G
#$ -l h_rt=8:00:00

source /netapp/home/mthomp/phenix-installer-1.13-2998-intel-linux-2.6-x86_64-centos5/phenix_env.sh

phenix.refine %s/%s %s refinement.modify_start_model.modify.sites.shake=1.0 refinement.main.number_of_macro_cycles=20 refinement.main.ordered_solvent=True refinement.ordered_solvent.mode=second_half refinement.target_weights.optimize_xyz_weight=True refinement.target_weights.optimize_adp_weight=True refinement.main.random_seed=%s output.prefix=%s_run1_refine

date
hostname
""" % (run1_path, run1_path, run1_path, rb_output_pdb_filename, output_mtz_file, rand_seed_1, original_mtz_filenames[n])
  qsub_script.write(phenix_refine_qsub_text_1)
  qsub_script.close()
  os.system("qsub refinement.sub")
  os.chdir(path)

  os.chdir(run2_path)
  qsub_script_fullpath = os.path.join(run2_path, "refinement.sub")
  qsub_script = open(qsub_script_fullpath, 'w')
  rand_seed_2 = randint(1000000, 9999999)
  phenix_refine_qsub_text_2 = """#$ -S /bin/bash
#$ -o %s/out
#$ -e %s/err
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=8G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=1G
#$ -l h_rt=8:00:00

source /netapp/home/mthomp/phenix-installer-1.13-2998-intel-linux-2.6-x86_64-centos5/phenix_env.sh

phenix.refine %s/%s %s refinement.modify_start_model.modify.sites.shake=1.0 refinement.main.number_of_macro_cycles=20 refinement.main.ordered_solvent=True refinement.ordered_solvent.mode=second_half refinement.target_weights.optimize_xyz_weight=True refinement.target_weights.optimize_adp_weight=True refinement.main.random_seed=%s output.prefix=%s_run2_refine

date
hostname
""" % (run2_path, run2_path, run2_path, rb_output_pdb_filename, output_mtz_file, rand_seed_2, original_mtz_filenames[n])
  qsub_script.write(phenix_refine_qsub_text_2)
  qsub_script.close()
  os.system("qsub refinement.sub")
  os.chdir(path)

  os.chdir(run3_path)
  qsub_script_fullpath = os.path.join(run3_path, "refinement.sub")
  qsub_script = open(qsub_script_fullpath, 'w')
  rand_seed_3 = randint(1000000, 9999999)
  phenix_refine_qsub_text_3 = """#$ -S /bin/bash
#$ -o %s/out
#$ -e %s/err
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=8G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=1G
#$ -l h_rt=8:00:00

source /netapp/home/mthomp/phenix-installer-1.13-2998-intel-linux-2.6-x86_64-centos5/phenix_env.sh

phenix.refine %s/%s %s refinement.modify_start_model.modify.sites.shake=1.0 refinement.main.number_of_macro_cycles=20 refinement.main.ordered_solvent=True refinement.ordered_solvent.mode=second_half refinement.target_weights.optimize_xyz_weight=True refinement.target_weights.optimize_adp_weight=True refinement.main.random_seed=%s output.prefix=%s_run3_refine

date
hostname
""" % (run3_path, run3_path, run3_path, rb_output_pdb_filename, output_mtz_file, rand_seed_3, original_mtz_filenames[n])
  qsub_script.write(phenix_refine_qsub_text_3)
  qsub_script.close()
  os.system("qsub refinement.sub")
  os.chdir(path)

  os.chdir(run4_path)
  qsub_script_fullpath = os.path.join(run4_path, "refinement.sub")
  qsub_script = open(qsub_script_fullpath, 'w')
  rand_seed_4 = randint(1000000, 9999999)
  phenix_refine_qsub_text_4 = """#$ -S /bin/bash
#$ -o %s/out
#$ -e %s/err
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=8G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=1G
#$ -l h_rt=8:00:00

source /netapp/home/mthomp/phenix-installer-1.13-2998-intel-linux-2.6-x86_64-centos5/phenix_env.sh

phenix.refine %s/%s %s refinement.modify_start_model.modify.sites.shake=1.0 refinement.main.number_of_macro_cycles=20 refinement.main.ordered_solvent=True refinement.ordered_solvent.mode=second_half refinement.target_weights.optimize_xyz_weight=True refinement.target_weights.optimize_adp_weight=True refinement.main.random_seed=%s output.prefix=%s_run4_refine

date
hostname
""" % (run4_path, run4_path, run4_path, rb_output_pdb_filename, output_mtz_file, rand_seed_4, original_mtz_filenames[n])
  qsub_script.write(phenix_refine_qsub_text_4)
  qsub_script.close()
  os.system("qsub refinement.sub")
  os.chdir(path)

  os.chdir(run5_path)
  qsub_script_fullpath = os.path.join(run5_path, "refinement.sub")
  qsub_script = open(qsub_script_fullpath, 'w')
  rand_seed_5 = randint(1000000, 9999999)
  phenix_refine_qsub_text_5 = """#$ -S /bin/bash
#$ -o %s/out
#$ -e %s/err
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=8G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=1G
#$ -l h_rt=8:00:00

source /netapp/home/mthomp/phenix-installer-1.13-2998-intel-linux-2.6-x86_64-centos5/phenix_env.sh

phenix.refine %s/%s %s refinement.modify_start_model.modify.sites.shake=1.0 refinement.main.number_of_macro_cycles=20 refinement.main.ordered_solvent=True refinement.ordered_solvent.mode=second_half refinement.target_weights.optimize_xyz_weight=True refinement.target_weights.optimize_adp_weight=True refinement.main.random_seed=%s output.prefix=%s_run5_refine

date
hostname
""" % (run5_path, run5_path, run5_path, rb_output_pdb_filename, output_mtz_file, rand_seed_5, original_mtz_filenames[n])
  qsub_script.write(phenix_refine_qsub_text_5)
  qsub_script.close()
  os.system("qsub refinement.sub")
  os.chdir(path)

  os.chdir(run6_path)
  qsub_script_fullpath = os.path.join(run6_path, "refinement.sub")
  qsub_script = open(qsub_script_fullpath, 'w')
  rand_seed_6 = randint(1000000, 9999999)
  phenix_refine_qsub_text_6 = """#$ -S /bin/bash
#$ -o %s/out
#$ -e %s/err
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=8G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=1G
#$ -l h_rt=8:00:00

source /netapp/home/mthomp/phenix-installer-1.13-2998-intel-linux-2.6-x86_64-centos5/phenix_env.sh

phenix.refine %s/%s %s refinement.modify_start_model.modify.sites.shake=1.0 refinement.main.number_of_macro_cycles=20 refinement.main.ordered_solvent=True refinement.ordered_solvent.mode=second_half refinement.target_weights.optimize_xyz_weight=True refinement.target_weights.optimize_adp_weight=True refinement.main.random_seed=%s output.prefix=%s_run6_refine

date
hostname
""" % (run6_path, run6_path, run6_path, rb_output_pdb_filename, output_mtz_file, rand_seed_6, original_mtz_filenames[n])
  qsub_script.write(phenix_refine_qsub_text_6)
  qsub_script.close()
  os.system("qsub refinement.sub")
  os.chdir(path)

  os.chdir(run7_path)
  qsub_script_fullpath = os.path.join(run7_path, "refinement.sub")
  qsub_script = open(qsub_script_fullpath, 'w')
  rand_seed_7 = randint(1000000, 9999999)
  phenix_refine_qsub_text_7 = """#$ -S /bin/bash
#$ -o %s/out
#$ -e %s/err
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=8G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=1G
#$ -l h_rt=8:00:00

source /netapp/home/mthomp/phenix-installer-1.13-2998-intel-linux-2.6-x86_64-centos5/phenix_env.sh

phenix.refine %s/%s %s refinement.modify_start_model.modify.sites.shake=1.0 refinement.main.number_of_macro_cycles=20 refinement.main.ordered_solvent=True refinement.ordered_solvent.mode=second_half refinement.target_weights.optimize_xyz_weight=True refinement.target_weights.optimize_adp_weight=True refinement.main.random_seed=%s output.prefix=%s_run7_refine

date
hostname
""" % (run7_path, run7_path, run7_path, rb_output_pdb_filename, output_mtz_file, rand_seed_7, original_mtz_filenames[n])
  qsub_script.write(phenix_refine_qsub_text_7)
  qsub_script.close()
  os.system("qsub refinement.sub")
  os.chdir(path)

  os.chdir(run8_path)
  qsub_script_fullpath = os.path.join(run8_path, "refinement.sub")
  qsub_script = open(qsub_script_fullpath, 'w')
  rand_seed_8 = randint(1000000, 9999999)
  phenix_refine_qsub_text_8 = """#$ -S /bin/bash
#$ -o %s/out
#$ -e %s/err
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=8G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=1G
#$ -l h_rt=8:00:00

source /netapp/home/mthomp/phenix-installer-1.13-2998-intel-linux-2.6-x86_64-centos5/phenix_env.sh

phenix.refine %s/%s %s refinement.modify_start_model.modify.sites.shake=1.0 refinement.main.number_of_macro_cycles=20 refinement.main.ordered_solvent=True refinement.ordered_solvent.mode=second_half refinement.target_weights.optimize_xyz_weight=True refinement.target_weights.optimize_adp_weight=True refinement.main.random_seed=%s output.prefix=%s_run8_refine

date
hostname
""" % (run8_path, run8_path, run8_path, rb_output_pdb_filename, output_mtz_file, rand_seed_8, original_mtz_filenames[n])
  qsub_script.write(phenix_refine_qsub_text_8)
  qsub_script.close()
  os.system("qsub refinement.sub")
  os.chdir(path)

  os.chdir(run9_path)
  qsub_script_fullpath = os.path.join(run9_path, "refinement.sub")
  qsub_script = open(qsub_script_fullpath, 'w')
  rand_seed_9 = randint(1000000, 9999999)
  phenix_refine_qsub_text_9 = """#$ -S /bin/bash
#$ -o %s/out
#$ -e %s/err
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=8G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=1G
#$ -l h_rt=8:00:00

source /netapp/home/mthomp/phenix-installer-1.13-2998-intel-linux-2.6-x86_64-centos5/phenix_env.sh

phenix.refine %s/%s %s refinement.modify_start_model.modify.sites.shake=1.0 refinement.main.number_of_macro_cycles=20 refinement.main.ordered_solvent=True refinement.ordered_solvent.mode=second_half refinement.target_weights.optimize_xyz_weight=True refinement.target_weights.optimize_adp_weight=True refinement.main.random_seed=%s output.prefix=%s_run9_refine

date
hostname
""" % (run9_path, run9_path, run9_path, rb_output_pdb_filename, output_mtz_file, rand_seed_9, original_mtz_filenames[n])
  qsub_script.write(phenix_refine_qsub_text_9)
  qsub_script.close()
  os.system("qsub refinement.sub")
  os.chdir(path)

  os.chdir(run10_path)
  qsub_script_fullpath = os.path.join(run10_path, "refinement.sub")
  qsub_script = open(qsub_script_fullpath, 'w')
  rand_seed_10 = randint(1000000, 9999999)
  phenix_refine_qsub_text_10 = """#$ -S /bin/bash
#$ -o %s/out
#$ -e %s/err
#$ -cwd
#$ -r y
#$ -j y
#$ -l mem_free=8G
#$ -l arch=linux-x64
#$ -l netapp=1G,scratch=1G
#$ -l h_rt=8:00:00

source /netapp/home/mthomp/phenix-installer-1.13-2998-intel-linux-2.6-x86_64-centos5/phenix_env.sh

phenix.refine %s/%s %s refinement.modify_start_model.modify.sites.shake=1.0 refinement.main.number_of_macro_cycles=20 refinement.main.ordered_solvent=True refinement.ordered_solvent.mode=second_half refinement.target_weights.optimize_xyz_weight=True refinement.target_weights.optimize_adp_weight=True refinement.main.random_seed=%s output.prefix=%s_run10_refine

date
hostname
""" % (run10_path, run10_path, run10_path, rb_output_pdb_filename, output_mtz_file, rand_seed_10, original_mtz_filenames[n])
  qsub_script.write(phenix_refine_qsub_text_10)
  qsub_script.close()
  os.system("qsub refinement.sub")
  os.chdir(path)


  os.chdir(starting_directory)


