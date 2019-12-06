#!/bin/sh
  f2mtz hklin /mnt/c/Users/mctuc/Documents/GitHub/TR_crystallography/test_diff_ref/difference_refinement/test_diffSF.hkl hklout /mnt/c/Users/mctuc/Documents/GitHub/TR_crystallography/test_diff_ref/difference_refinement/test_diffSF.mtz <<EOF
  SYMMETRY 96
  LABOUT H  K  L  FOBS  SIGFOBS
  CTYPE  H  H  H  F  Q
  EOF
  