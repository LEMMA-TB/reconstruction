#ifndef _FILES_H_
#define _FILES_H_

#define YEAR 2014
#define MONTH 11

#define MAX_ENTRIES 100000000

#define N_RUNS 6

#define RUN_NUMB 100450
#define RUN_NUMB2 100455
#define RUN_NUMB3 100457
#define RUN_NUMB4 100459
#define RUN_NUMB5 100461
#define RUN_NUMB6 100463
#define RUN_NUMB7 100484
#define RUN_NUMB8 100486
#define RUN_NUMB9 100487
#define RUN_NUMB10 100489
#define RUN_NUMB11 100492
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
#define CRY 1 //on instead of position on telescope use position on crystal
#define SAVE 0 //
#define ONLYTREE 0//

#define CORRECT_TORSION 1 //move the xpos and xpos zeros at the center of the crystal
#define CENTER 1 //move the xpos and xpos zeros at the center of the crystal
#define GONIO_ZERO 1//set the zero of the gonio plot in the first step
#define OFFSET_GONIO 1//

#endif
