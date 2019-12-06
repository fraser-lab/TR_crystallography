import os
import sys
import math
import string
import subprocess
import numpy as np 
import matplotlib.pyplot as plt
import operator


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


def calculate_Fo_minus_Fo_as_dict( ground_state_F_as_text_fullpath, excited_state_F_as_text_fullpath ):

  ground_state_F_as_dict = store_Fhkl_from_text_as_dict(ground_state_F_as_text_fullpath)
  excited_state_F_as_dict = store_Fhkl_from_text_as_dict(excited_state_F_as_text_fullpath)
  common_reflections = identify_common_reflections(ground_state_F_as_dict, excited_state_F_as_dict)
  Fo_minus_Fo_dict = {}
  for miller_index in common_reflections:
    fobs_1 = ground_state_F_as_dict.get(miller_index)[0]
    fobs_2 = excited_state_F_as_dict.get(miller_index)[0]
    Fo_minus_Fo = float(fobs_2) - float(fobs_1)
    values = {miller_index : Fo_minus_Fo}
    Fo_minus_Fo_dict.update(values)

  return Fo_minus_Fo_dict

def create_Fo_minus_Fo_rank_order_list( cwd, Fo_minus_Fo_dict ):

  output_file_path = os.path.join(cwd, "Fo_minus_Fo_rank_order.txt")
  output_file = open(output_file_path, 'w')
  ranked_Fo_minus_Fo = sorted(Fo_minus_Fo_dict.items(), key=operator.itemgetter(1))
  for difference in ranked_Fo_minus_Fo:
    output_line = "%s %s \n" % (difference[0], difference[1])
    output_file.write(output_line)
  output_file.close()


def generate_Fo_minus_Fo_histogram( Fo_minus_Fo_dict ):

  Fo_minus_Fo_list = []
  for Fo_minus_Fo in Fo_minus_Fo_dict.itervalues():
    Fo_minus_Fo_list.append(Fo_minus_Fo)
  Fo_minus_Fo_array = np.array(Fo_minus_Fo_list)
  plt.hist(Fo_minus_Fo_array, 500)
  plt.savefig("Fo_minus_Fo_histogram.png")
  plt.show()



def main():

  cwd = os.getcwd()

  input_reflections1_as_text_fullpath = os.path.join(cwd, "lys_200us_540uJ_flags_mtz.txt")
  input_reflections2_as_text_fullpath = os.path.join(cwd, "lys_20ns_540uJ_flags_mtz.txt")
  
  out = calculate_Fo_minus_Fo_as_dict( input_reflections1_as_text_fullpath, input_reflections2_as_text_fullpath )

  create_Fo_minus_Fo_rank_order_list( cwd, out )
  generate_Fo_minus_Fo_histogram( out )

  
main()
 