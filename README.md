# AHA data processing
 Matlab code for data processing in the AHA project

## What is needed to run the software
1. Matlab.  Any version should work. Tests were done with Matlab 2020a.
2. Opensim.  This is only needed for converting the C3D files to TXT files. Tests were done with Opensim 4.0, but later versions will probably work also.
3. Opensim-Matlab interface. To set this up, follow the instructions on https://simtk-confluence.stanford.edu/display/OpenSim/Scripting+with+Matlab.
4. Not strictly needed, but we recommend installing Notepad++ to read and edit text files, such as Mocap files and Log files.

## C3D to TXT file conversion
This is needed when marker data has been edited offline using Cortex. Cortex will save the data as a C3D file.  Marker data from the C3D file must be converted into the Mocap TXT file format.

### c3dtotxt.m 
c3dtotxt.m is a function to convert one file. You need to give it the name (full path) of the C3D file, and the name of the TXT Mocap file that was originally recorded on D-Flow. The original file has the force plate data that is used for the synchronization. This function will extract the marker data from the C3D file, synchronize it with the TXT file, and generate a new TXT file, which contains the force plate data from the original TXT file and the marker data from the C3D file.  The new TXT file has the same name as the original TXT file, but the .txt extension is replaced by _edited.txt.

For testing this function, there is a script c3dtest.m that converts trial 1 from partipant 4 (POST test), and compares the converted file to the original file. This shows that the data is the same and the synchronization worked correctly. In this trial, the C3D data started recording 60 ms after the D-Flow recording, and ended 1050 ms earlier.  So the new TXT file is a little shorter than the original TXT file. However, the time stamps and synchronization are correct, so there is no problem combining this with the data on the treadmill file (which also has time stamps).

### Batch conversion
c3dbatch.m is a script that is set up to converts all c3d trials in a folder, or in a series of folders. To run this, the data needs to be organized as follows:
* The C3D files must be named the same as the original TXT files, but with the .c3d extension e.g. Mocap0001.c3d.
* The C3d Files must be in the same folder as the corresponding TXT files.
* All folders to be be converted must be subfolders of one "data path" folder.  In c3dbatch.m, you should set up this data path for your computer.

During the conversion, a file c3dbatch.log is written with a record of what was done.  This file also reports, for each file, how much marker data was missing in the original TXT file, and how much is missing in the new TXT file. One marker missing in one frame is counted as one piece of missing data.  Inspect this file to make sure that everything went correctly.

If a _edited.txt file already exists for a trial, the conversion is skipped for that trial, because it was already done. If you want to redo the conversion for that trial, simply delete the _edited.txt file before running c3dbatch.m.

### Possible issues
* The C3D conversion has not been tested yet on a trial where the C3D recording started *before* the Mocap file.
* The C3D conversion will not work yet on trials where the TXT file does not have labeled data for all 47 markers.  This can easily be fixed.

## Analysis of normal gait
(describe the code and the results that are produced)

## Analysis of perturbation response
(describe the code and the results that are produced)

## Script for automated data processing of all trials
(describe the code and how the results are stored/presented)
