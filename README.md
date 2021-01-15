# AHA data processing
 Matlab code for data processing in the AHA project

## What is needed to run the software
1. Matlab.  Any version should work. Tests were done with Matlab 2020a.
2. Opensim.  This is only needed for converting the C3D files to TXT files. Tests were done with Opensim 4.0, but later versions will probably work also.
3. Opensim-Matlab interface. To set this up, follow the instructions on https://simtk-confluence.stanford.edu/display/OpenSim/Scripting+with+Matlab.

## C3D to TXT file conversion
There is a function c3dtotxt.m to convert one file. You need to give it the name (full path) of the C3D file, and the name of the TXT Mocap file that was originally recorded on D-Flow. The original file has the force plate data that is used for the synchronization. This function will extract the marker data from the C3D file, synchronize it with the TXT file, and generate a new TXT file, which contains the force plate data from the original TXT file and the marker data from the C3D file.  The new TXT file has the same name as the C3D file, but the .c3d extension is replaced by .txt.

To test that the conversion works correctly, there is a script c3dtest.m that converts trial 1 from partipant 4 (POST test), and compares the converted file to the original file. This shows that the data is the same and the synchronization worked correctly. In this file, the C3D data started recording 60 ms after the D-Flow recording, and ended 1050 ms earlier.  So the new TXT file starts and ends with some frames where the marker data is all zero (missing).

To run the conversion automatically, there is a script c3dconvert.m that converts all trials in a folder, or in a series of folders. (this is still in development).

## Analysis of normal gait
(describe the code and the results that are produced)

## Analysis of perturbation response
(describe the code and the results that are produced)

## Script for automated data processing of all trials
(describe the code and how the results are stored/presented)
