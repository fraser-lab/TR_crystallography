import string
import math
import os


names = ["260", "265", "270", "275", "280", "260-2", "260-3", "270-2", "270-3", "280-2", "280-3"]
temps = [260, 265, 270, 275, 280, 260, 260, 270, 270, 280, 280]
path_to_images = "/Volumes/DatumsDepot/2018/Mike/......."
common_prefix = ""
base_dir = os.getcwd()
ref_data_set = ""

for name in names:
  dirname = os.path.join(base_dir, name)
  os.mkdir(dirname, 0755)
  w1_dirname = os.path.join(dirname, "wedge1")
  os.mkdir(w1_dirname, 0755)
  w2_dirname = os.path.join(dirname, "wedge2")
  os.mkdir(w2_dirname, 0755)
  merge_dirname = os.path.join(dirname, "merge")
  os.mkdir(merge_dirname, 0755)

  image_template_w1 = path_to_images+common_prefix+name+"_1_00001.cbf"
  image_template_w2 = path_to_images+common_prefix+name+"_2_00001.cbf"
  w1_path_to_hkl = w1_dirname+"/XDS_ASCII.HKL"
  w2_path_to_hkl = w2_dirname+"/XDS_ASCII.HKL"


  os.chdir(w1_dirname)
  inp1_file_fullpath = os.path.join(w1_dirname, "XDS.INP")
  inp1_file = open(inp1_file_fullpath, 'w')
 
  inp1_text = """!*****************************************************************************
 DETECTOR=PILATUS         MINIMUM_VALID_PIXEL_VALUE=0  OVERLOAD=1048500
 SENSOR_THICKNESS=0.45        !SILICON=-1.0
 GAIN=1.0
!AIR=0.001 !Air absorption coefficient of x-rays is computed by XDS by default
!NX=number of fast pixels (along X); QX=length of an X-pixel (mm)
!NY=number of slow pixels (along Y); QY=length of a  Y-pixel (mm)
 NX=2463 NY=2527 QX=0.172  QY=0.172  !PILATUS 6M
  UNTRUSTED_ELLIPSE=1184 1289    1218 1322  ! ellipse enclosed by X1 X2 Y1 Y2
 UNTRUSTED_RECTANGLE= 487  495     0 2528
 UNTRUSTED_RECTANGLE= 981  989     0 2528
 UNTRUSTED_RECTANGLE=1475 1483     0 2528
 UNTRUSTED_RECTANGLE=1969 1977     0 2528
 UNTRUSTED_RECTANGLE=   0 2464   195  213
 UNTRUSTED_RECTANGLE=   0 2464   407  425
 UNTRUSTED_RECTANGLE=   0 2464   619  637
 UNTRUSTED_RECTANGLE=   0 2464   831  849
 UNTRUSTED_RECTANGLE=   0 2464  1043 1061
 UNTRUSTED_RECTANGLE=   0 2464  1255 1273
 UNTRUSTED_RECTANGLE=   0 2464  1467 1485
 UNTRUSTED_RECTANGLE=   0 2464  1679 1697
 UNTRUSTED_RECTANGLE=   0 2464  1891 1909
 UNTRUSTED_RECTANGLE=   0 2464  2103 2121
 UNTRUSTED_RECTANGLE=   0 2464  2315 2333
!UNTRUSTED_QUADRILATERAL=565 574  1519 1552  1508 1533  566 1536
!MINIMUM_FRACTION_OF_BACKGROUND_REGION=0.01
!X-GEO_CORR=GD_6M_X06SA_SLS_27022007_X.pck CCP4
!Y-GEO_CORR=GD_6M_X06SA_SLS_27022007_Y.pck CCP4

 DIRECTION_OF_DETECTOR_X-AXIS= 1.0 0.0 0.0
 DIRECTION_OF_DETECTOR_Y-AXIS= 0.0 1.0 0.0 ! 0.0 cos(2theta) sin(2theta)
 TRUSTED_REGION=0.0 1.5 !Relative radii limiting trusted detector region

!====================== JOB CONTROL PARAMETERS ===============================
 JOB= DEFPIX INTEGRATE CORRECT

!MAXIMUM_NUMBER_OF_JOBS=4  !Speeds-up COLSPOT & INTEGRATE on a Linux-cluster
!MAXIMUM_NUMBER_OF_PROCESSORS=8!<100;ignored by single cpu version of xds
!CLUSTER_NODES=bragg01 bragg02 bragg03 bragg04 bragg05 bragg06 bragg07 bragg08
!SECONDS=0  !Maximum number of seconds to wait until data image must appear
!TEST=1     !Test flag. 1,2 additional diagnostics and images

!====================== GEOMETRICAL PARAMETERS ===============================
!ORGX and ORGY are often close to the image center, i.e. ORGX=NX/2, ORGY=NY/2
 ORGX= 1233.318115 ORGY= 1253.498535   !Detector origin (pixels).  ORGX=NX/2; ORGY=NY/2
 DETECTOR_DISTANCE= 250.0   !(mm)

 ROTATION_AXIS= 1.0 0.0 0.0

! Optimal choice is 0.5*mosaicity (REFLECTING_RANGE_E.S.D.= mosaicity)
 OSCILLATION_RANGE=0.5            !degrees (>0)

 X-RAY_WAVELENGTH=1.11583          !Angstroem
 INCIDENT_BEAM_DIRECTION=0.0 0.0 1.0
 FRACTION_OF_POLARIZATION=0.98 !default=0.5 for unpolarized beam
 POLARIZATION_PLANE_NORMAL= 0.0 1.0 0.0

!======================= CRYSTAL PARAMETERS =================================

  SPACE_GROUP_NUMBER=0   !0 for unknown crystals; cell constants are ignored.
! UNIT_CELL_CONSTANTS=157.5  162.8  597.2     90 90 90

! You may specify here the x,y,z components for the unit cell vectors if
! known from a previous run using the same crystal in the same orientation
!UNIT_CELL_A-AXIS= 
!UNIT_CELL_B-AXIS=
!UNIT_CELL_C-AXIS=

!Optional reindexing transformation to apply on reflection indices
!REIDX=   0  0 -1  0  0 -1  0  0 -1  0  0  0


!==================== SELECTION OF DATA IMAGES ==============================
!Generic file name and format (optional) of data images
 NAME_TEMPLATE_OF_DATA_FRAMES= %s

 DATA_RANGE=1 5      !Numbers of first and last data image collected

 BACKGROUND_RANGE=1 5 !Numbers of first and last data image for background

 SPOT_RANGE=1 5       !First and last data image number for finding spots

!==================== DATA COLLECTION STRATEGY (XPLAN) ======================
!                       !!! Warning !!!
! If you processed your data for a crystal with unknown cell constants and
! space group symmetry, XPLAN will report the results for space group P1.

!STARTING_ANGLE=  0.0      STARTING_FRAME=1
!used to define the angular origin about the rotation axis.
!Default:  STARTING_ANGLE=  0 at STARTING_FRAME=first data image

!STARTING_ANGLES_OF_SPINDLE_ROTATION= 0 180 10
!TOTAL_SPINDLE_ROTATION_RANGES=30.0 120 15
!RESOLUTION_SHELLS=10 6 5 4 3 2 1.5 1.3 1.2

!====================== INDEXING PARAMETERS =================================
!Additional parameters for fine tuning that rarely need to be changed
!RGRID=0.0005     ! grid size for difference vector histogram (1/Angstrom)
 SEPMIN=4 CLUSTER_RADIUS=2
!INDEX_ERROR=0.05 INDEX_MAGNITUDE=8 INDEX_QUALITY=0.8
!MERGE_TREE=0.1   ! for unknown crystals
!INDEX_ORIGIN= 0 0 0 ! used by "IDXREF" to add a final index offset !!!
!MAXIMUM_ERROR_OF_SPOT_POSITION=3.0
!MAXIMUM_ERROR_OF_SPINDLE_POSITION=2.0
!MINIMUM_FRACTION_OF_INDEXED_SPOTS=0.5

!============== DECISION CONSTANTS FOR FINDING CRYSTAL SYMMETRY =============
!Decision constants for detection of lattice symmetry (IDXREF, CORRECT)
!MAX_CELL_AXIS_ERROR=0.03 ! Maximum relative error in cell axes tolerated 
!MAX_CELL_ANGLE_ERROR=2.0 ! Maximum cell angle error tolerated

!Decision constants for detection of space group symmetry (CORRECT). 
!Resolution range for accepting reflections for space group determination in
!the CORRECT step. It should cover a sufficient number of strong reflections.
 TEST_RESOLUTION_RANGE=8.0 1.5 
!MIN_RFL_Rmeas= 50 ! Minimum #reflections needed for calculation of Rmeas
!MAX_FAC_Rmeas=2.0 ! Sets an upper limit for acceptable Rmeas
 
!================= PARAMETERS CONTROLLING REFINEMENTS =======================
!REFINE(IDXREF)=BEAM AXIS ORIENTATION CELL !POSITION
!REFINE(INTEGRATE)=!POSITION BEAM ORIENTATION ! CELL !AXIS
!REFINE(CORRECT)=POSITION BEAM ORIENTATION CELL AXIS

!================== CRITERIA FOR ACCEPTING REFLECTIONS ======================
 VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS= 6000 30000 !Used by DEFPIX
           !for excluding shaded parts of the detector.

 INCLUDE_RESOLUTION_RANGE=50.0  1.4  !Angstroem; used by DEFPIX,INTEGRATE,CORRECT

!used by CORRECT to exclude ice-reflections
!EXCLUDE_RESOLUTION_RANGE= 3.93 3.87 !ice-ring at 3.897 Angstrom
!EXCLUDE_RESOLUTION_RANGE= 3.70 3.64 !ice-ring at 3.669 Angstrom
!EXCLUDE_RESOLUTION_RANGE= 3.47 3.41 !ice-ring at 3.441 Angstrom
!EXCLUDE_RESOLUTION_RANGE= 2.70 2.64 !ice-ring at 2.671 Angstrom
!EXCLUDE_RESOLUTION_RANGE= 2.28 2.22 !ice-ring at 2.249 Angstrom
!EXCLUDE_RESOLUTION_RANGE= 2.102 2.042 !ice-ring at 2.072 Angstrom - strong
!EXCLUDE_RESOLUTION_RANGE= 1.978 1.918 !ice-ring at 1.948 Angstrom - weak
!EXCLUDE_RESOLUTION_RANGE= 1.948 1.888 !ice-ring at 1.918 Angstrom - strong
!EXCLUDE_RESOLUTION_RANGE= 1.913 1.853 !ice-ring at 1.883 Angstrom - weak
!EXCLUDE_RESOLUTION_RANGE= 1.751 1.691 !ice-ring at 1.721 Angstrom - weak

!MINIMUM_ZETA=0.05 !Defines width of 'blind region' (XPLAN,INTEGRATE,CORRECT)

!WFAC1=1.0  !This controls the number of rejected MISFITS in CORRECT;
      !a larger value leads to fewer rejections.
!REJECT_ALIEN=20.0 ! Automatic rejection of very strong reflections

!============== INTEGRATION AND PEAK PROFILE PARAMETERS =====================
!Specification of the peak profile parameters below overrides the automatic
!determination from the images
!Suggested values are listed near the end of INTEGRATE.LP
!BEAM_DIVERGENCE=   0.80         !arctan(spot diameter/DETECTOR_DISTANCE)
!BEAM_DIVERGENCE_E.S.D.=   0.080 !half-width (Sigma) of BEAM_DIVERGENCE
!RELRAD=5.0 ! automatic :: BEAM_DIVERGENCE=2*RELRAD*BEAM_DIVERGENCE_E.S.D.
!REFLECTING_RANGE=  0.780 !for crossing the Ewald sphere on shortest route
!REFLECTING_RANGE_E.S.D.=  0.113 !half-width (mosaicity) of REFLECTING_RANGE
!RELBET=3.5 ! automatic :: REFLECTING_RANGE=2*RELBET*REFLECTING_RANGE_E.S.D.
 NUMBER_OF_PROFILE_GRID_POINTS_ALONG_ALPHA/BETA=13!used by: INTEGRATE
 NUMBER_OF_PROFILE_GRID_POINTS_ALONG_GAMMA= 13     !used by: INTEGRATE
!PROFILE_FITTING=TRUE 
!CUT=2.0    !defines the integration region for profile fitting
!MINPK=75.0 !minimum required percentage of observed reflection intensity
!DELPHI= 6.0!controls the number of reference profiles and scaling factors

!======= PARAMETERS CONTROLLING CORRECTION FACTORS (used by: CORRECT) =======
!FRIEDEL'S_LAW=FALSE !Default is TRUE.
!MINIMUM_I/SIGMA=3.0 !minimum intensity/sigma required for scaling reflections
!NBATCH=-1  !controls the number of correction factors along image numbers
!REFLECTIONS/CORRECTION_FACTOR=50   !minimum #reflections/correction needed
!PATCH_SHUTTER_PROBLEM=TRUE         !FALSE is default
!STRICT_ABSORPTION_CORRECTION=TRUE  !FALSE is default
!CORRECTIONS= DECAY MODULATION ABSORPTION
 REFERENCE_DATA_SET= %s   !Name of a reference data set (optional)
!FIT_B-FACTOR_TO_REFERENCE_DATA_SET=TRUE ! default is FALSE

!=========== PARAMETERS DEFINING BACKGROUND AND PEAK PIXELS =================
!STRONG_PIXEL=3.0                              !used by: COLSPOT
!A 'strong' pixel to be included in a spot must exceed the background
!by more than the given multiple of standard deviations.

!MAXIMUM_NUMBER_OF_STRONG_PIXELS=1500000       !used by: COLSPOT

!SPOT_MAXIMUM-CENTROID=3.0                     !used by: COLSPOT

 MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT=3          !used by: COLSPOT
!This allows to suppress spurious isolated pixels from entering the
!spot list generated by "COLSPOT".

!NBX=3  NBY=3  !Define a rectangle of size (2*NBX+1)*(2*NBY+1)
!The variation of counts within the rectangle centered at each image pixel
!is used for distinguishing between background and spot pixels.

!BACKGROUND_PIXEL=6.0                          !used by: COLSPOT,INTEGRATE
!An image pixel does not belong to the background region if the local
!pixel variation exceeds the expected variation by the given number of
!standard deviations.

!SIGNAL_PIXEL=3.0                              !used by: INTEGRATE
!A pixel above the threshold contributes to the spot centroid

 DATA_RANGE_FIXED_SCALE_FACTOR= 1 360 1.0 ! used by : INIT,INTEGRATE

""" % (image_template_w1, ref_data_set)
  
  inp1_file.write(inp1_text)
  inp1_file.close()
  os.system('xds')

  inp1_recycle_file_fullpath = os.path.join(w1_dirname, "XDS.INP")
  inp1_recycle_file = open(inp1_recycle_file_fullpath, 'w')
 
  inp1_recycle_text = """!*****************************************************************************
 DETECTOR=PILATUS         MINIMUM_VALID_PIXEL_VALUE=0  OVERLOAD=1048500
 SENSOR_THICKNESS=0.45        !SILICON=-1.0
 GAIN=1.0
!AIR=0.001 !Air absorption coefficient of x-rays is computed by XDS by default
!NX=number of fast pixels (along X); QX=length of an X-pixel (mm)
!NY=number of slow pixels (along Y); QY=length of a  Y-pixel (mm)
 NX=2463 NY=2527 QX=0.172  QY=0.172  !PILATUS 6M
  UNTRUSTED_ELLIPSE=1184 1289    1218 1322  ! ellipse enclosed by X1 X2 Y1 Y2
 UNTRUSTED_RECTANGLE= 487  495     0 2528
 UNTRUSTED_RECTANGLE= 981  989     0 2528
 UNTRUSTED_RECTANGLE=1475 1483     0 2528
 UNTRUSTED_RECTANGLE=1969 1977     0 2528
 UNTRUSTED_RECTANGLE=   0 2464   195  213
 UNTRUSTED_RECTANGLE=   0 2464   407  425
 UNTRUSTED_RECTANGLE=   0 2464   619  637
 UNTRUSTED_RECTANGLE=   0 2464   831  849
 UNTRUSTED_RECTANGLE=   0 2464  1043 1061
 UNTRUSTED_RECTANGLE=   0 2464  1255 1273
 UNTRUSTED_RECTANGLE=   0 2464  1467 1485
 UNTRUSTED_RECTANGLE=   0 2464  1679 1697
 UNTRUSTED_RECTANGLE=   0 2464  1891 1909
 UNTRUSTED_RECTANGLE=   0 2464  2103 2121
 UNTRUSTED_RECTANGLE=   0 2464  2315 2333
!UNTRUSTED_QUADRILATERAL=565 574  1519 1552  1508 1533  566 1536
!MINIMUM_FRACTION_OF_BACKGROUND_REGION=0.01
!X-GEO_CORR=GD_6M_X06SA_SLS_27022007_X.pck CCP4
!Y-GEO_CORR=GD_6M_X06SA_SLS_27022007_Y.pck CCP4

 DIRECTION_OF_DETECTOR_X-AXIS= 1.0 0.0 0.0
 DIRECTION_OF_DETECTOR_Y-AXIS= 0.0 1.0 0.0 ! 0.0 cos(2theta) sin(2theta)
 TRUSTED_REGION=0.0 1.5 !Relative radii limiting trusted detector region

!====================== JOB CONTROL PARAMETERS ===============================
!JOB= XYCORR INIT COLSPOT IDXREF DEFPIX XPLAN INTEGRATE CORRECT

!MAXIMUM_NUMBER_OF_JOBS=4  !Speeds-up COLSPOT & INTEGRATE on a Linux-cluster
!MAXIMUM_NUMBER_OF_PROCESSORS=8!<100;ignored by single cpu version of xds
!CLUSTER_NODES=bragg01 bragg02 bragg03 bragg04 bragg05 bragg06 bragg07 bragg08
!SECONDS=0  !Maximum number of seconds to wait until data image must appear
!TEST=1     !Test flag. 1,2 additional diagnostics and images

!====================== GEOMETRICAL PARAMETERS ===============================
!ORGX and ORGY are often close to the image center, i.e. ORGX=NX/2, ORGY=NY/2
 ORGX= 1233.318115 ORGY= 1253.498535   !Detector origin (pixels).  ORGX=NX/2; ORGY=NY/2
 DETECTOR_DISTANCE= 250.0   !(mm)

 ROTATION_AXIS= 1.0 0.0 0.0

! Optimal choice is 0.5*mosaicity (REFLECTING_RANGE_E.S.D.= mosaicity)
 OSCILLATION_RANGE=0.5            !degrees (>0)

 X-RAY_WAVELENGTH=1.11583          !Angstroem
 INCIDENT_BEAM_DIRECTION=0.0 0.0 1.0
 FRACTION_OF_POLARIZATION=0.98 !default=0.5 for unpolarized beam
 POLARIZATION_PLANE_NORMAL= 0.0 1.0 0.0

!======================= CRYSTAL PARAMETERS =================================

  SPACE_GROUP_NUMBER=0   !0 for unknown crystals; cell constants are ignored.
! UNIT_CELL_CONSTANTS=157.5  162.8  597.2     90 90 90

! You may specify here the x,y,z components for the unit cell vectors if
! known from a previous run using the same crystal in the same orientation
!UNIT_CELL_A-AXIS= 
!UNIT_CELL_B-AXIS=
!UNIT_CELL_C-AXIS=

!Optional reindexing transformation to apply on reflection indices
!REIDX=   0  0 -1  0  0 -1  0  0 -1  0  0  0


!==================== SELECTION OF DATA IMAGES ==============================
!Generic file name and format (optional) of data images
 NAME_TEMPLATE_OF_DATA_FRAMES= %s

 DATA_RANGE=1 5      !Numbers of first and last data image collected

 BACKGROUND_RANGE=1 5 !Numbers of first and last data image for background

 SPOT_RANGE=1 5       !First and last data image number for finding spots

!==================== DATA COLLECTION STRATEGY (XPLAN) ======================
!                       !!! Warning !!!
! If you processed your data for a crystal with unknown cell constants and
! space group symmetry, XPLAN will report the results for space group P1.

!STARTING_ANGLE=  0.0      STARTING_FRAME=1
!used to define the angular origin about the rotation axis.
!Default:  STARTING_ANGLE=  0 at STARTING_FRAME=first data image

!STARTING_ANGLES_OF_SPINDLE_ROTATION= 0 180 10
!TOTAL_SPINDLE_ROTATION_RANGES=30.0 120 15
!RESOLUTION_SHELLS=10 6 5 4 3 2 1.5 1.3 1.2

!====================== INDEXING PARAMETERS =================================
!Additional parameters for fine tuning that rarely need to be changed
!RGRID=0.0005     ! grid size for difference vector histogram (1/Angstrom)
 SEPMIN=4 CLUSTER_RADIUS=2
!INDEX_ERROR=0.05 INDEX_MAGNITUDE=8 INDEX_QUALITY=0.8
!MERGE_TREE=0.1   ! for unknown crystals
!INDEX_ORIGIN= 0 0 0 ! used by "IDXREF" to add a final index offset !!!
!MAXIMUM_ERROR_OF_SPOT_POSITION=3.0
!MAXIMUM_ERROR_OF_SPINDLE_POSITION=2.0
!MINIMUM_FRACTION_OF_INDEXED_SPOTS=0.5

!============== DECISION CONSTANTS FOR FINDING CRYSTAL SYMMETRY =============
!Decision constants for detection of lattice symmetry (IDXREF, CORRECT)
!MAX_CELL_AXIS_ERROR=0.03 ! Maximum relative error in cell axes tolerated 
!MAX_CELL_ANGLE_ERROR=2.0 ! Maximum cell angle error tolerated

!Decision constants for detection of space group symmetry (CORRECT). 
!Resolution range for accepting reflections for space group determination in
!the CORRECT step. It should cover a sufficient number of strong reflections.
 TEST_RESOLUTION_RANGE=8.0 1.5 
!MIN_RFL_Rmeas= 50 ! Minimum #reflections needed for calculation of Rmeas
!MAX_FAC_Rmeas=2.0 ! Sets an upper limit for acceptable Rmeas
 
!================= PARAMETERS CONTROLLING REFINEMENTS =======================
!REFINE(IDXREF)=BEAM AXIS ORIENTATION CELL !POSITION
!REFINE(INTEGRATE)=!POSITION BEAM ORIENTATION ! CELL !AXIS
!REFINE(CORRECT)=POSITION BEAM ORIENTATION CELL AXIS

!================== CRITERIA FOR ACCEPTING REFLECTIONS ======================
 VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS= 6000 30000 !Used by DEFPIX
           !for excluding shaded parts of the detector.

 INCLUDE_RESOLUTION_RANGE=50.0  1.4  !Angstroem; used by DEFPIX,INTEGRATE,CORRECT

!used by CORRECT to exclude ice-reflections
!EXCLUDE_RESOLUTION_RANGE= 3.93 3.87 !ice-ring at 3.897 Angstrom
!EXCLUDE_RESOLUTION_RANGE= 3.70 3.64 !ice-ring at 3.669 Angstrom
!EXCLUDE_RESOLUTION_RANGE= 3.47 3.41 !ice-ring at 3.441 Angstrom
!EXCLUDE_RESOLUTION_RANGE= 2.70 2.64 !ice-ring at 2.671 Angstrom
!EXCLUDE_RESOLUTION_RANGE= 2.28 2.22 !ice-ring at 2.249 Angstrom
!EXCLUDE_RESOLUTION_RANGE= 2.102 2.042 !ice-ring at 2.072 Angstrom - strong
!EXCLUDE_RESOLUTION_RANGE= 1.978 1.918 !ice-ring at 1.948 Angstrom - weak
!EXCLUDE_RESOLUTION_RANGE= 1.948 1.888 !ice-ring at 1.918 Angstrom - strong
!EXCLUDE_RESOLUTION_RANGE= 1.913 1.853 !ice-ring at 1.883 Angstrom - weak
!EXCLUDE_RESOLUTION_RANGE= 1.751 1.691 !ice-ring at 1.721 Angstrom - weak

!MINIMUM_ZETA=0.05 !Defines width of 'blind region' (XPLAN,INTEGRATE,CORRECT)

!WFAC1=1.0  !This controls the number of rejected MISFITS in CORRECT;
      !a larger value leads to fewer rejections.
!REJECT_ALIEN=20.0 ! Automatic rejection of very strong reflections

!============== INTEGRATION AND PEAK PROFILE PARAMETERS =====================
!Specification of the peak profile parameters below overrides the automatic
!determination from the images
!Suggested values are listed near the end of INTEGRATE.LP
!BEAM_DIVERGENCE=   0.80         !arctan(spot diameter/DETECTOR_DISTANCE)
!BEAM_DIVERGENCE_E.S.D.=   0.080 !half-width (Sigma) of BEAM_DIVERGENCE
!RELRAD=5.0 ! automatic :: BEAM_DIVERGENCE=2*RELRAD*BEAM_DIVERGENCE_E.S.D.
!REFLECTING_RANGE=  0.780 !for crossing the Ewald sphere on shortest route
!REFLECTING_RANGE_E.S.D.=  0.113 !half-width (mosaicity) of REFLECTING_RANGE
!RELBET=3.5 ! automatic :: REFLECTING_RANGE=2*RELBET*REFLECTING_RANGE_E.S.D.
 NUMBER_OF_PROFILE_GRID_POINTS_ALONG_ALPHA/BETA=13!used by: INTEGRATE
 NUMBER_OF_PROFILE_GRID_POINTS_ALONG_GAMMA= 13     !used by: INTEGRATE
!PROFILE_FITTING=TRUE 
!CUT=2.0    !defines the integration region for profile fitting
!MINPK=75.0 !minimum required percentage of observed reflection intensity
!DELPHI= 6.0!controls the number of reference profiles and scaling factors

!======= PARAMETERS CONTROLLING CORRECTION FACTORS (used by: CORRECT) =======
!FRIEDEL'S_LAW=FALSE !Default is TRUE.
!MINIMUM_I/SIGMA=3.0 !minimum intensity/sigma required for scaling reflections
!NBATCH=-1  !controls the number of correction factors along image numbers
!REFLECTIONS/CORRECTION_FACTOR=50   !minimum #reflections/correction needed
!PATCH_SHUTTER_PROBLEM=TRUE         !FALSE is default
!STRICT_ABSORPTION_CORRECTION=TRUE  !FALSE is default
!CORRECTIONS= DECAY MODULATION ABSORPTION
 REFERENCE_DATA_SET= %s   !Name of a reference data set (optional)
!FIT_B-FACTOR_TO_REFERENCE_DATA_SET=TRUE ! default is FALSE

!=========== PARAMETERS DEFINING BACKGROUND AND PEAK PIXELS =================
!STRONG_PIXEL=3.0                              !used by: COLSPOT
!A 'strong' pixel to be included in a spot must exceed the background
!by more than the given multiple of standard deviations.

!MAXIMUM_NUMBER_OF_STRONG_PIXELS=1500000       !used by: COLSPOT

!SPOT_MAXIMUM-CENTROID=3.0                     !used by: COLSPOT

 MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT=3          !used by: COLSPOT
!This allows to suppress spurious isolated pixels from entering the
!spot list generated by "COLSPOT".

!NBX=3  NBY=3  !Define a rectangle of size (2*NBX+1)*(2*NBY+1)
!The variation of counts within the rectangle centered at each image pixel
!is used for distinguishing between background and spot pixels.

!BACKGROUND_PIXEL=6.0                          !used by: COLSPOT,INTEGRATE
!An image pixel does not belong to the background region if the local
!pixel variation exceeds the expected variation by the given number of
!standard deviations.

!SIGNAL_PIXEL=3.0                              !used by: INTEGRATE
!A pixel above the threshold contributes to the spot centroid

 DATA_RANGE_FIXED_SCALE_FACTOR= 1 360 1.0 ! used by : INIT,INTEGRATE

""" % (image_template_w1, ref_data_set)
  
  inp1_recycle_file.write(inp1_recycle_text)
  inp1_recycle_file.close()

  os.system("mv ./GXPARM.XDS ./XPARM.XDS")
  os.system('xds')
  os.system("mv ./GXPARM.XDS ./XPARM.XDS")
  os.system('xds')
  os.system("mv ./GXPARM.XDS ./XPARM.XDS")
  os.system('xds')

  os.chdir(w2_dirname)
  inp2_file_fullpath = os.path.join(w2_dirname, "XDS.INP")
  inp2_file = open(inp2_file_fullpath, 'w')
 
  inp2_text = """!*****************************************************************************
 DETECTOR=PILATUS         MINIMUM_VALID_PIXEL_VALUE=0  OVERLOAD=1048500
 SENSOR_THICKNESS=0.45        !SILICON=-1.0
 GAIN=1.0
!AIR=0.001 !Air absorption coefficient of x-rays is computed by XDS by default
!NX=number of fast pixels (along X); QX=length of an X-pixel (mm)
!NY=number of slow pixels (along Y); QY=length of a  Y-pixel (mm)
 NX=2463 NY=2527 QX=0.172  QY=0.172  !PILATUS 6M
  UNTRUSTED_ELLIPSE=1184 1289    1218 1322  ! ellipse enclosed by X1 X2 Y1 Y2
 UNTRUSTED_RECTANGLE= 487  495     0 2528
 UNTRUSTED_RECTANGLE= 981  989     0 2528
 UNTRUSTED_RECTANGLE=1475 1483     0 2528
 UNTRUSTED_RECTANGLE=1969 1977     0 2528
 UNTRUSTED_RECTANGLE=   0 2464   195  213
 UNTRUSTED_RECTANGLE=   0 2464   407  425
 UNTRUSTED_RECTANGLE=   0 2464   619  637
 UNTRUSTED_RECTANGLE=   0 2464   831  849
 UNTRUSTED_RECTANGLE=   0 2464  1043 1061
 UNTRUSTED_RECTANGLE=   0 2464  1255 1273
 UNTRUSTED_RECTANGLE=   0 2464  1467 1485
 UNTRUSTED_RECTANGLE=   0 2464  1679 1697
 UNTRUSTED_RECTANGLE=   0 2464  1891 1909
 UNTRUSTED_RECTANGLE=   0 2464  2103 2121
 UNTRUSTED_RECTANGLE=   0 2464  2315 2333
!UNTRUSTED_QUADRILATERAL=565 574  1519 1552  1508 1533  566 1536
!MINIMUM_FRACTION_OF_BACKGROUND_REGION=0.01
!X-GEO_CORR=GD_6M_X06SA_SLS_27022007_X.pck CCP4
!Y-GEO_CORR=GD_6M_X06SA_SLS_27022007_Y.pck CCP4

 DIRECTION_OF_DETECTOR_X-AXIS= 1.0 0.0 0.0
 DIRECTION_OF_DETECTOR_Y-AXIS= 0.0 1.0 0.0 ! 0.0 cos(2theta) sin(2theta)
 TRUSTED_REGION=0.0 1.5 !Relative radii limiting trusted detector region

!====================== JOB CONTROL PARAMETERS ===============================
 JOB= DEFPIX INTEGRATE CORRECT

!MAXIMUM_NUMBER_OF_JOBS=4  !Speeds-up COLSPOT & INTEGRATE on a Linux-cluster
!MAXIMUM_NUMBER_OF_PROCESSORS=8!<100;ignored by single cpu version of xds
!CLUSTER_NODES=bragg01 bragg02 bragg03 bragg04 bragg05 bragg06 bragg07 bragg08
!SECONDS=0  !Maximum number of seconds to wait until data image must appear
!TEST=1     !Test flag. 1,2 additional diagnostics and images

!====================== GEOMETRICAL PARAMETERS ===============================
!ORGX and ORGY are often close to the image center, i.e. ORGX=NX/2, ORGY=NY/2
 ORGX= 1233.318115 ORGY= 1253.498535   !Detector origin (pixels).  ORGX=NX/2; ORGY=NY/2
 DETECTOR_DISTANCE= 250.0   !(mm)

 ROTATION_AXIS= 1.0 0.0 0.0

! Optimal choice is 0.5*mosaicity (REFLECTING_RANGE_E.S.D.= mosaicity)
 OSCILLATION_RANGE=0.5            !degrees (>0)

 X-RAY_WAVELENGTH=1.11583          !Angstroem
 INCIDENT_BEAM_DIRECTION=0.0 0.0 1.0
 FRACTION_OF_POLARIZATION=0.98 !default=0.5 for unpolarized beam
 POLARIZATION_PLANE_NORMAL= 0.0 1.0 0.0

!======================= CRYSTAL PARAMETERS =================================

  SPACE_GROUP_NUMBER=0   !0 for unknown crystals; cell constants are ignored.
! UNIT_CELL_CONSTANTS=157.5  162.8  597.2     90 90 90

! You may specify here the x,y,z components for the unit cell vectors if
! known from a previous run using the same crystal in the same orientation
!UNIT_CELL_A-AXIS= 
!UNIT_CELL_B-AXIS=
!UNIT_CELL_C-AXIS=

!Optional reindexing transformation to apply on reflection indices
!REIDX=   0  0 -1  0  0 -1  0  0 -1  0  0  0


!==================== SELECTION OF DATA IMAGES ==============================
!Generic file name and format (optional) of data images
 NAME_TEMPLATE_OF_DATA_FRAMES= %s

 DATA_RANGE=1 5      !Numbers of first and last data image collected

 BACKGROUND_RANGE=1 5 !Numbers of first and last data image for background

 SPOT_RANGE=1 5       !First and last data image number for finding spots

!==================== DATA COLLECTION STRATEGY (XPLAN) ======================
!                       !!! Warning !!!
! If you processed your data for a crystal with unknown cell constants and
! space group symmetry, XPLAN will report the results for space group P1.

!STARTING_ANGLE=  0.0      STARTING_FRAME=1
!used to define the angular origin about the rotation axis.
!Default:  STARTING_ANGLE=  0 at STARTING_FRAME=first data image

!STARTING_ANGLES_OF_SPINDLE_ROTATION= 0 180 10
!TOTAL_SPINDLE_ROTATION_RANGES=30.0 120 15
!RESOLUTION_SHELLS=10 6 5 4 3 2 1.5 1.3 1.2

!====================== INDEXING PARAMETERS =================================
!Additional parameters for fine tuning that rarely need to be changed
!RGRID=0.0005     ! grid size for difference vector histogram (1/Angstrom)
 SEPMIN=4 CLUSTER_RADIUS=2
!INDEX_ERROR=0.05 INDEX_MAGNITUDE=8 INDEX_QUALITY=0.8
!MERGE_TREE=0.1   ! for unknown crystals
!INDEX_ORIGIN= 0 0 0 ! used by "IDXREF" to add a final index offset !!!
!MAXIMUM_ERROR_OF_SPOT_POSITION=3.0
!MAXIMUM_ERROR_OF_SPINDLE_POSITION=2.0
!MINIMUM_FRACTION_OF_INDEXED_SPOTS=0.5

!============== DECISION CONSTANTS FOR FINDING CRYSTAL SYMMETRY =============
!Decision constants for detection of lattice symmetry (IDXREF, CORRECT)
!MAX_CELL_AXIS_ERROR=0.03 ! Maximum relative error in cell axes tolerated 
!MAX_CELL_ANGLE_ERROR=2.0 ! Maximum cell angle error tolerated

!Decision constants for detection of space group symmetry (CORRECT). 
!Resolution range for accepting reflections for space group determination in
!the CORRECT step. It should cover a sufficient number of strong reflections.
 TEST_RESOLUTION_RANGE=8.0 1.5 
!MIN_RFL_Rmeas= 50 ! Minimum #reflections needed for calculation of Rmeas
!MAX_FAC_Rmeas=2.0 ! Sets an upper limit for acceptable Rmeas
 
!================= PARAMETERS CONTROLLING REFINEMENTS =======================
!REFINE(IDXREF)=BEAM AXIS ORIENTATION CELL !POSITION
!REFINE(INTEGRATE)=!POSITION BEAM ORIENTATION ! CELL !AXIS
!REFINE(CORRECT)=POSITION BEAM ORIENTATION CELL AXIS

!================== CRITERIA FOR ACCEPTING REFLECTIONS ======================
 VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS= 6000 30000 !Used by DEFPIX
           !for excluding shaded parts of the detector.

 INCLUDE_RESOLUTION_RANGE=50.0  1.4  !Angstroem; used by DEFPIX,INTEGRATE,CORRECT

!used by CORRECT to exclude ice-reflections
!EXCLUDE_RESOLUTION_RANGE= 3.93 3.87 !ice-ring at 3.897 Angstrom
!EXCLUDE_RESOLUTION_RANGE= 3.70 3.64 !ice-ring at 3.669 Angstrom
!EXCLUDE_RESOLUTION_RANGE= 3.47 3.41 !ice-ring at 3.441 Angstrom
!EXCLUDE_RESOLUTION_RANGE= 2.70 2.64 !ice-ring at 2.671 Angstrom
!EXCLUDE_RESOLUTION_RANGE= 2.28 2.22 !ice-ring at 2.249 Angstrom
!EXCLUDE_RESOLUTION_RANGE= 2.102 2.042 !ice-ring at 2.072 Angstrom - strong
!EXCLUDE_RESOLUTION_RANGE= 1.978 1.918 !ice-ring at 1.948 Angstrom - weak
!EXCLUDE_RESOLUTION_RANGE= 1.948 1.888 !ice-ring at 1.918 Angstrom - strong
!EXCLUDE_RESOLUTION_RANGE= 1.913 1.853 !ice-ring at 1.883 Angstrom - weak
!EXCLUDE_RESOLUTION_RANGE= 1.751 1.691 !ice-ring at 1.721 Angstrom - weak

!MINIMUM_ZETA=0.05 !Defines width of 'blind region' (XPLAN,INTEGRATE,CORRECT)

!WFAC1=1.0  !This controls the number of rejected MISFITS in CORRECT;
      !a larger value leads to fewer rejections.
!REJECT_ALIEN=20.0 ! Automatic rejection of very strong reflections

!============== INTEGRATION AND PEAK PROFILE PARAMETERS =====================
!Specification of the peak profile parameters below overrides the automatic
!determination from the images
!Suggested values are listed near the end of INTEGRATE.LP
!BEAM_DIVERGENCE=   0.80         !arctan(spot diameter/DETECTOR_DISTANCE)
!BEAM_DIVERGENCE_E.S.D.=   0.080 !half-width (Sigma) of BEAM_DIVERGENCE
!RELRAD=5.0 ! automatic :: BEAM_DIVERGENCE=2*RELRAD*BEAM_DIVERGENCE_E.S.D.
!REFLECTING_RANGE=  0.780 !for crossing the Ewald sphere on shortest route
!REFLECTING_RANGE_E.S.D.=  0.113 !half-width (mosaicity) of REFLECTING_RANGE
!RELBET=3.5 ! automatic :: REFLECTING_RANGE=2*RELBET*REFLECTING_RANGE_E.S.D.
 NUMBER_OF_PROFILE_GRID_POINTS_ALONG_ALPHA/BETA=13!used by: INTEGRATE
 NUMBER_OF_PROFILE_GRID_POINTS_ALONG_GAMMA= 13     !used by: INTEGRATE
!PROFILE_FITTING=TRUE 
!CUT=2.0    !defines the integration region for profile fitting
!MINPK=75.0 !minimum required percentage of observed reflection intensity
!DELPHI= 6.0!controls the number of reference profiles and scaling factors

!======= PARAMETERS CONTROLLING CORRECTION FACTORS (used by: CORRECT) =======
!FRIEDEL'S_LAW=FALSE !Default is TRUE.
!MINIMUM_I/SIGMA=3.0 !minimum intensity/sigma required for scaling reflections
!NBATCH=-1  !controls the number of correction factors along image numbers
!REFLECTIONS/CORRECTION_FACTOR=50   !minimum #reflections/correction needed
!PATCH_SHUTTER_PROBLEM=TRUE         !FALSE is default
!STRICT_ABSORPTION_CORRECTION=TRUE  !FALSE is default
!CORRECTIONS= DECAY MODULATION ABSORPTION
 REFERENCE_DATA_SET= %s   !Name of a reference data set (optional)
!FIT_B-FACTOR_TO_REFERENCE_DATA_SET=TRUE ! default is FALSE

!=========== PARAMETERS DEFINING BACKGROUND AND PEAK PIXELS =================
!STRONG_PIXEL=3.0                              !used by: COLSPOT
!A 'strong' pixel to be included in a spot must exceed the background
!by more than the given multiple of standard deviations.

!MAXIMUM_NUMBER_OF_STRONG_PIXELS=1500000       !used by: COLSPOT

!SPOT_MAXIMUM-CENTROID=3.0                     !used by: COLSPOT

 MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT=3          !used by: COLSPOT
!This allows to suppress spurious isolated pixels from entering the
!spot list generated by "COLSPOT".

!NBX=3  NBY=3  !Define a rectangle of size (2*NBX+1)*(2*NBY+1)
!The variation of counts within the rectangle centered at each image pixel
!is used for distinguishing between background and spot pixels.

!BACKGROUND_PIXEL=6.0                          !used by: COLSPOT,INTEGRATE
!An image pixel does not belong to the background region if the local
!pixel variation exceeds the expected variation by the given number of
!standard deviations.

!SIGNAL_PIXEL=3.0                              !used by: INTEGRATE
!A pixel above the threshold contributes to the spot centroid

 DATA_RANGE_FIXED_SCALE_FACTOR= 1 360 1.0 ! used by : INIT,INTEGRATE

""" % (image_template_w2, ref_data_set)
  
  inp2_file.write(inp2_text)
  inp2_file.close()
  os.system('xds')

  inp2_recycle_file_fullpath = os.path.join(w1_dirname, "XDS.INP")
  inp2_recycle_file = open(inp2_recycle_file_fullpath, 'w')
 
  inp2_recycle_text = """!*****************************************************************************
 DETECTOR=PILATUS         MINIMUM_VALID_PIXEL_VALUE=0  OVERLOAD=1048500
 SENSOR_THICKNESS=0.45        !SILICON=-1.0
 GAIN=1.0
!AIR=0.001 !Air absorption coefficient of x-rays is computed by XDS by default
!NX=number of fast pixels (along X); QX=length of an X-pixel (mm)
!NY=number of slow pixels (along Y); QY=length of a  Y-pixel (mm)
 NX=2463 NY=2527 QX=0.172  QY=0.172  !PILATUS 6M
  UNTRUSTED_ELLIPSE=1184 1289    1218 1322  ! ellipse enclosed by X1 X2 Y1 Y2
 UNTRUSTED_RECTANGLE= 487  495     0 2528
 UNTRUSTED_RECTANGLE= 981  989     0 2528
 UNTRUSTED_RECTANGLE=1475 1483     0 2528
 UNTRUSTED_RECTANGLE=1969 1977     0 2528
 UNTRUSTED_RECTANGLE=   0 2464   195  213
 UNTRUSTED_RECTANGLE=   0 2464   407  425
 UNTRUSTED_RECTANGLE=   0 2464   619  637
 UNTRUSTED_RECTANGLE=   0 2464   831  849
 UNTRUSTED_RECTANGLE=   0 2464  1043 1061
 UNTRUSTED_RECTANGLE=   0 2464  1255 1273
 UNTRUSTED_RECTANGLE=   0 2464  1467 1485
 UNTRUSTED_RECTANGLE=   0 2464  1679 1697
 UNTRUSTED_RECTANGLE=   0 2464  1891 1909
 UNTRUSTED_RECTANGLE=   0 2464  2103 2121
 UNTRUSTED_RECTANGLE=   0 2464  2315 2333
!UNTRUSTED_QUADRILATERAL=565 574  1519 1552  1508 1533  566 1536
!MINIMUM_FRACTION_OF_BACKGROUND_REGION=0.01
!X-GEO_CORR=GD_6M_X06SA_SLS_27022007_X.pck CCP4
!Y-GEO_CORR=GD_6M_X06SA_SLS_27022007_Y.pck CCP4

 DIRECTION_OF_DETECTOR_X-AXIS= 1.0 0.0 0.0
 DIRECTION_OF_DETECTOR_Y-AXIS= 0.0 1.0 0.0 ! 0.0 cos(2theta) sin(2theta)
 TRUSTED_REGION=0.0 1.5 !Relative radii limiting trusted detector region

!====================== JOB CONTROL PARAMETERS ===============================
!JOB= XYCORR INIT COLSPOT IDXREF DEFPIX XPLAN INTEGRATE CORRECT

!MAXIMUM_NUMBER_OF_JOBS=4  !Speeds-up COLSPOT & INTEGRATE on a Linux-cluster
!MAXIMUM_NUMBER_OF_PROCESSORS=8!<100;ignored by single cpu version of xds
!CLUSTER_NODES=bragg01 bragg02 bragg03 bragg04 bragg05 bragg06 bragg07 bragg08
!SECONDS=0  !Maximum number of seconds to wait until data image must appear
!TEST=1     !Test flag. 1,2 additional diagnostics and images

!====================== GEOMETRICAL PARAMETERS ===============================
!ORGX and ORGY are often close to the image center, i.e. ORGX=NX/2, ORGY=NY/2
 ORGX= 1233.318115 ORGY= 1253.498535   !Detector origin (pixels).  ORGX=NX/2; ORGY=NY/2
 DETECTOR_DISTANCE= 250.0   !(mm)

 ROTATION_AXIS= 1.0 0.0 0.0

! Optimal choice is 0.5*mosaicity (REFLECTING_RANGE_E.S.D.= mosaicity)
 OSCILLATION_RANGE=0.5            !degrees (>0)

 X-RAY_WAVELENGTH=1.11583          !Angstroem
 INCIDENT_BEAM_DIRECTION=0.0 0.0 1.0
 FRACTION_OF_POLARIZATION=0.98 !default=0.5 for unpolarized beam
 POLARIZATION_PLANE_NORMAL= 0.0 1.0 0.0

!======================= CRYSTAL PARAMETERS =================================

  SPACE_GROUP_NUMBER=0   !0 for unknown crystals; cell constants are ignored.
! UNIT_CELL_CONSTANTS=157.5  162.8  597.2     90 90 90

! You may specify here the x,y,z components for the unit cell vectors if
! known from a previous run using the same crystal in the same orientation
!UNIT_CELL_A-AXIS= 
!UNIT_CELL_B-AXIS=
!UNIT_CELL_C-AXIS=

!Optional reindexing transformation to apply on reflection indices
!REIDX=   0  0 -1  0  0 -1  0  0 -1  0  0  0


!==================== SELECTION OF DATA IMAGES ==============================
!Generic file name and format (optional) of data images
 NAME_TEMPLATE_OF_DATA_FRAMES= %s

 DATA_RANGE=1 5      !Numbers of first and last data image collected

 BACKGROUND_RANGE=1 5 !Numbers of first and last data image for background

 SPOT_RANGE=1 5       !First and last data image number for finding spots

!==================== DATA COLLECTION STRATEGY (XPLAN) ======================
!                       !!! Warning !!!
! If you processed your data for a crystal with unknown cell constants and
! space group symmetry, XPLAN will report the results for space group P1.

!STARTING_ANGLE=  0.0      STARTING_FRAME=1
!used to define the angular origin about the rotation axis.
!Default:  STARTING_ANGLE=  0 at STARTING_FRAME=first data image

!STARTING_ANGLES_OF_SPINDLE_ROTATION= 0 180 10
!TOTAL_SPINDLE_ROTATION_RANGES=30.0 120 15
!RESOLUTION_SHELLS=10 6 5 4 3 2 1.5 1.3 1.2

!====================== INDEXING PARAMETERS =================================
!Additional parameters for fine tuning that rarely need to be changed
!RGRID=0.0005     ! grid size for difference vector histogram (1/Angstrom)
 SEPMIN=4 CLUSTER_RADIUS=2
!INDEX_ERROR=0.05 INDEX_MAGNITUDE=8 INDEX_QUALITY=0.8
!MERGE_TREE=0.1   ! for unknown crystals
!INDEX_ORIGIN= 0 0 0 ! used by "IDXREF" to add a final index offset !!!
!MAXIMUM_ERROR_OF_SPOT_POSITION=3.0
!MAXIMUM_ERROR_OF_SPINDLE_POSITION=2.0
!MINIMUM_FRACTION_OF_INDEXED_SPOTS=0.5

!============== DECISION CONSTANTS FOR FINDING CRYSTAL SYMMETRY =============
!Decision constants for detection of lattice symmetry (IDXREF, CORRECT)
!MAX_CELL_AXIS_ERROR=0.03 ! Maximum relative error in cell axes tolerated 
!MAX_CELL_ANGLE_ERROR=2.0 ! Maximum cell angle error tolerated

!Decision constants for detection of space group symmetry (CORRECT). 
!Resolution range for accepting reflections for space group determination in
!the CORRECT step. It should cover a sufficient number of strong reflections.
 TEST_RESOLUTION_RANGE=8.0 1.5 
!MIN_RFL_Rmeas= 50 ! Minimum #reflections needed for calculation of Rmeas
!MAX_FAC_Rmeas=2.0 ! Sets an upper limit for acceptable Rmeas
 
!================= PARAMETERS CONTROLLING REFINEMENTS =======================
!REFINE(IDXREF)=BEAM AXIS ORIENTATION CELL !POSITION
!REFINE(INTEGRATE)=!POSITION BEAM ORIENTATION ! CELL !AXIS
!REFINE(CORRECT)=POSITION BEAM ORIENTATION CELL AXIS

!================== CRITERIA FOR ACCEPTING REFLECTIONS ======================
 VALUE_RANGE_FOR_TRUSTED_DETECTOR_PIXELS= 6000 30000 !Used by DEFPIX
           !for excluding shaded parts of the detector.

 INCLUDE_RESOLUTION_RANGE=50.0  1.4  !Angstroem; used by DEFPIX,INTEGRATE,CORRECT

!used by CORRECT to exclude ice-reflections
!EXCLUDE_RESOLUTION_RANGE= 3.93 3.87 !ice-ring at 3.897 Angstrom
!EXCLUDE_RESOLUTION_RANGE= 3.70 3.64 !ice-ring at 3.669 Angstrom
!EXCLUDE_RESOLUTION_RANGE= 3.47 3.41 !ice-ring at 3.441 Angstrom
!EXCLUDE_RESOLUTION_RANGE= 2.70 2.64 !ice-ring at 2.671 Angstrom
!EXCLUDE_RESOLUTION_RANGE= 2.28 2.22 !ice-ring at 2.249 Angstrom
!EXCLUDE_RESOLUTION_RANGE= 2.102 2.042 !ice-ring at 2.072 Angstrom - strong
!EXCLUDE_RESOLUTION_RANGE= 1.978 1.918 !ice-ring at 1.948 Angstrom - weak
!EXCLUDE_RESOLUTION_RANGE= 1.948 1.888 !ice-ring at 1.918 Angstrom - strong
!EXCLUDE_RESOLUTION_RANGE= 1.913 1.853 !ice-ring at 1.883 Angstrom - weak
!EXCLUDE_RESOLUTION_RANGE= 1.751 1.691 !ice-ring at 1.721 Angstrom - weak

!MINIMUM_ZETA=0.05 !Defines width of 'blind region' (XPLAN,INTEGRATE,CORRECT)

!WFAC1=1.0  !This controls the number of rejected MISFITS in CORRECT;
      !a larger value leads to fewer rejections.
!REJECT_ALIEN=20.0 ! Automatic rejection of very strong reflections

!============== INTEGRATION AND PEAK PROFILE PARAMETERS =====================
!Specification of the peak profile parameters below overrides the automatic
!determination from the images
!Suggested values are listed near the end of INTEGRATE.LP
!BEAM_DIVERGENCE=   0.80         !arctan(spot diameter/DETECTOR_DISTANCE)
!BEAM_DIVERGENCE_E.S.D.=   0.080 !half-width (Sigma) of BEAM_DIVERGENCE
!RELRAD=5.0 ! automatic :: BEAM_DIVERGENCE=2*RELRAD*BEAM_DIVERGENCE_E.S.D.
!REFLECTING_RANGE=  0.780 !for crossing the Ewald sphere on shortest route
!REFLECTING_RANGE_E.S.D.=  0.113 !half-width (mosaicity) of REFLECTING_RANGE
!RELBET=3.5 ! automatic :: REFLECTING_RANGE=2*RELBET*REFLECTING_RANGE_E.S.D.
 NUMBER_OF_PROFILE_GRID_POINTS_ALONG_ALPHA/BETA=13!used by: INTEGRATE
 NUMBER_OF_PROFILE_GRID_POINTS_ALONG_GAMMA= 13     !used by: INTEGRATE
!PROFILE_FITTING=TRUE 
!CUT=2.0    !defines the integration region for profile fitting
!MINPK=75.0 !minimum required percentage of observed reflection intensity
!DELPHI= 6.0!controls the number of reference profiles and scaling factors

!======= PARAMETERS CONTROLLING CORRECTION FACTORS (used by: CORRECT) =======
!FRIEDEL'S_LAW=FALSE !Default is TRUE.
!MINIMUM_I/SIGMA=3.0 !minimum intensity/sigma required for scaling reflections
!NBATCH=-1  !controls the number of correction factors along image numbers
!REFLECTIONS/CORRECTION_FACTOR=50   !minimum #reflections/correction needed
!PATCH_SHUTTER_PROBLEM=TRUE         !FALSE is default
!STRICT_ABSORPTION_CORRECTION=TRUE  !FALSE is default
!CORRECTIONS= DECAY MODULATION ABSORPTION
 REFERENCE_DATA_SET= %s   !Name of a reference data set (optional)
!FIT_B-FACTOR_TO_REFERENCE_DATA_SET=TRUE ! default is FALSE

!=========== PARAMETERS DEFINING BACKGROUND AND PEAK PIXELS =================
!STRONG_PIXEL=3.0                              !used by: COLSPOT
!A 'strong' pixel to be included in a spot must exceed the background
!by more than the given multiple of standard deviations.

!MAXIMUM_NUMBER_OF_STRONG_PIXELS=1500000       !used by: COLSPOT

!SPOT_MAXIMUM-CENTROID=3.0                     !used by: COLSPOT

 MINIMUM_NUMBER_OF_PIXELS_IN_A_SPOT=3          !used by: COLSPOT
!This allows to suppress spurious isolated pixels from entering the
!spot list generated by "COLSPOT".

!NBX=3  NBY=3  !Define a rectangle of size (2*NBX+1)*(2*NBY+1)
!The variation of counts within the rectangle centered at each image pixel
!is used for distinguishing between background and spot pixels.

!BACKGROUND_PIXEL=6.0                          !used by: COLSPOT,INTEGRATE
!An image pixel does not belong to the background region if the local
!pixel variation exceeds the expected variation by the given number of
!standard deviations.

!SIGNAL_PIXEL=3.0                              !used by: INTEGRATE
!A pixel above the threshold contributes to the spot centroid

 DATA_RANGE_FIXED_SCALE_FACTOR= 1 360 1.0 ! used by : INIT,INTEGRATE

""" % (image_template_w2, ref_data_set)
  
  inp2_recycle_file.write(inp2_recycle_text)
  inp2_recycle_file.close()

  os.system("mv ./GXPARM.XDS ./XPARM.XDS")
  os.system('xds')
  os.system("mv ./GXPARM.XDS ./XPARM.XDS")
  os.system('xds')
  os.system("mv ./GXPARM.XDS ./XPARM.XDS")
  os.system('xds')

  os.chdir(merge_dirname)
  merge_file_fullpath = os.path.join(merge_dirname, "XSCALE.INP")
  merge_file = open(merge_file_fullpath, 'w')
 
  merge_text = """


""" & ()

  os.chdir(base_dir)