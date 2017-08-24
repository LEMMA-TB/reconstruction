#ifndef _PARAMTERS_H_
#define _PARAMETERS_H_

#define YEAR 2016
#define MONTH 7

#define MAX_ENTRIES 100000000000

#define N_RUNS 1

#define RUN_NUMB 112551
#define RUN_NUMB2 112597
#define RUN_NUMB3 112600
#define RUN_NUMB4 112410
#define RUN_NUMB5 112413
#define RUN_NUMB6 110103
#define RUN_NUMB7 110104
#define RUN_NUMB8 110105
#define RUN_NUMB9 110106
#define RUN_NUMB10 110107
#define RUN_NUMB11 110108
#define RUN_NUMB12 100497
#define RUN_NUMB13 100500
#define RUN_NUMB14 100501
#define RUN_NUMB15 97735
#define RUN_NUMB16 97736
#define RUN_NUMB17 97736

#define N_FILES -1 //# files for single run to be analyzed; -1 == all files
#define N_FILES_IN 0 //starting number of files to be analyzed

//// OPTIONS
#define SLICESON 0 // slices histograms 0=no 1=yes
#define LOWRES 0 //low resolution histograms 3d 0=off 1=on
#define CRY 0 //on instead of position on telescope use position on crystal
#define SAVE 0 //
#define ONLYTREE 0//

#define CORRECT_TORSION 0 //move the xpos and xpos zeros at the center of the crystal
#define CENTER 0 //move the xpos and xpos zeros at the center of the crystal
#define GONIO_ZERO 0//set the zero of the gonio plot in the first step
#define OFFSET_GONIO 0//
#define PLANE_ROTANG 0

///////////////////////////////////////////////////////////////////////
//////////////// Silicon Detector Position in meters///////////////////
///////////////////////////////////////////////////////////////////////
#define Z110_1   0.          // SD1
#define Z210_1   10.075      // SD2
#define Zg10_1   10.574      // Goniometer-Crystal position
#define Z310_1   20.678      // SD3
#define Z410_1   21.051      // SD4

#define Z112_1   00.000      // SD1
#define Z212_1   10.008      // SD2
#define Zg12_1   10.227      // Goniometer-Crystal position
#define Z312_1   21.178      // SD3

#define Z114_1   0.000      // SD1
#define Z214_1   5.37      // SD2
#define Zg14_1   6.28      // Goniometer-Crystal position
#define Z314_1   11.50      // SD3

#define Z115_1   0.000      // SD1
#define Z215_1   5.21      // SD2
#define Zg15_1   5.57      // Goniometer-Crystal position
#define Z315_1   12.03      // SD3

#define Z116_1   0.000      // SD1
#define Z216_1   5.6      // SD2
#define Zg16_1   6.09      // Goniometer-Crystal position
#define Z316_1   11.4      // SD3

#endif
