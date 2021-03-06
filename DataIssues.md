# Data Issues

This is the file DataIssues.md in the AHA-data-processing repository.

We can add notes here whenever we encouter any data issues that need to be looked at.

When an issue is completely resolved, just delete it!

Use level 2 headers (double hash sign) for each folder.  See below for example.

In some folders, trials 10 and higher must be renamed, so the trial number has exactly 4 digits.  For example:
Mocap00010.txt and Treadmill00010.txt must be renamed Mocap0010.txt and Treadmill0010.txt
We will be using T1-T7 trials only for further analysis.

Advantages:
* Ensures that trials are processed in the correct order when folder contents are sorted by name
* Makes it easier for Matlab to generate file names automatically
* If we always have three zeros, we might as well have no zeros.

Important notes: 
* We will be using T1-T7 trials only for further analysis.
* In case we had false start in T1, we will be using T2-T8 for further analysis.


## Par1_PRE
No C3D files yet on Onedrive
No Cortex files were recorded

## Par1_POST
No C3D files yet on Onedrive
No Cortex files were recorded

## Par2_PRE
No C3D files yet on Onedrive
No Cortex files were recorded

## Par2_POST
No C3D files yet on Onedrive
No Cortex files were recorded

The treadmill file names have a spelling error (change Treaddmill to Treadmill).
This issue has been fixed-files are renamed correctly now!

## Par3_PRE
No C3D files yet on Onedrive
No Cortex files were recorded

## Par3_POST
No C3D files yet on Onedrive
No Cortex files were recorded

## Par4_PRE
No C3D files yet on Onedrive

## Par4_POST
There were 15 C3D files, and they were converted without issues, but it seems that no data editing was done in Cortex.  The C3D files appear to have the same amount of data missing as the original TXT files.  
Files were edited and added recently on OneDrive 

## Par5_PRE
No C3D files yet on Onedrive
No Cortex files were recorded

## Par5_POST
No C3D files yet on Onedrive
No Cortex files were recorded

## Par6_PRE
Empty folder

## Par6_POST
Empty folder
No edit was made, we will use the original files on OneDrive 

## Par7_PRE
There were C3D files for trials 1-8 only.

In two files, some markers had less data in the C3D file than in the original TXT file.  The edited TXT file was created, but kept the data from the original TXT file (for those markers).

Par7_PRE\Mocap0001.c3d
* WARNING: C3D has more missing data than TXT for LSHO; C3D data not used.

Par7_PRE\Mocap0008.c3d
* WARNING: C3D has more missing data than TXT for C7; C3D data not used.
* WARNING: C3D has more missing data than TXT for LMEE; C3D data not used.
* WARNING: C3D has more missing data than TXT for LHEE; C3D data not used.
* WARNING: C3D has more missing data than TXT for LMT5; C3D data not used.


## Par7_POST

The folder contains C3D files for trials 1-8

Based on the force plate data, Mocap0002.c3d and Mocap0002.txt appear to be different trials.  Please check if something went wrong with the file naming, or possibly the C3D files came from a different participant or different session?

Also in Mocap0002.c3d, there are 10 markers for which the C3D file has less data than the TXT file.  That's not surprising if they are actually different trials.

The same two issues are present in all other trials with C3D files in this folder.

## Par8_PRE
Empty folder
Files were edited and added recently on OneDrive 


## Par8_POST
In trial 7 the TXT file is only 49.88 seconds long.  The C3D data is 90 seconds, and shows that the participant is not walking in the last 40 seconds.  Treadmill may have been stopped.

In the same trial, some marker data disagrees between the original TXT and the C3D file. Marker labeling seems to have changed.  Please verify that the C3D is good before using it.

It would be good to also check the other trials in this folder.
Par8_T7 was false start- will not be used in the analysis 
Files were edited and added recently on OneDrive 


## Par9_PRE
Empty folder
Files were edited and added recently on OneDrive 


## Par9_POST
Empty folder
Files were edited and added recently on OneDrive 


## Par10_PRE
Empty folder
No edit was made, we will use the original files on OneDrive 


## Par10_POST
Empty folder
No edit was made, we will use the original files on OneDrive 


## Par11_PRE
Empty folder
No edit was made, we will use the original files on OneDrive 


## Par11_POST
Empty folder
No edit was made, we will use the original files on OneDrive 


## Par12_PRE
Empty folder
No edit was made, we will use the original files on OneDrive 


## Par12_POST
No edit was made, we will use the original files on OneDrive 
Folder now has all the c3d files. I used the c3dbatch.m (seems to work fine and FAST!) 







