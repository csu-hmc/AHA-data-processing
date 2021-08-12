# AHA data processing
 Matlab code for data processing in the AHA project

## What is needed to run the software
1. Matlab.  Any version should work. Tests were done with Matlab 2020a.
2. Opensim.  This is only needed for converting the C3D files to TXT files. Tests were done with Opensim 4.0, but later versions will probably work also.
3. Opensim-Matlab interface. To set this up, follow the instructions on https://simtk-confluence.stanford.edu/display/OpenSim/Scripting+with+Matlab.
4. Not strictly needed, but we recommend installing Notepad++ to read and edit text files, such as Mocap files and Log files.
5. Software to manage the repository on your computer. For Windows, we recommend Github Desktop

If you have not done so yet, clone this Github repository.  On Windows, the repository will usually show up in your Documents/Github folder. Use the main repository folder, named AHA-data-processing, as the working directory in Matlab when running the code.

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

During the conversion, a file c3dbatch.log is written with a record of what was done.  This file also reports, for each file, how much marker data was missing in the original TXT file, and how much is missing in the new TXT file. One marker missing in one frame is counted as one piece of missing data.  Inspect this file to make sure that everything went correctly.  If something does not look good, report it.

If a _edited.txt file already exists for a trial, the conversion is skipped for that trial, because it was already done. If you want to redo the conversion for that trial, simply delete the _edited.txt file before running c3dbatch.m.

If the conversion fails for one file, try to fix the problem (it may be a bug) or ask Ton for help.  After the problem is fixed, simply restart c3dbatch.m. It will automatically skip those files that were already converted, and continue from where it stopped.

### Possible issues
This has only been tested on trial 1 in the Participant 4, Post-test folder.  It is possible that the conversion does not work for some other files.  If so, please report the problem, which trial it is, and what information is written on the screen.  Some possible issues that may occur:
* If the C3D recording started *before* the Mocap file, there is code to handle that, but this has not been tested because it did not occur in the test file.
* If the original TXT file does not have labeled data for all 47 markers, the conversion should still work but this has not been tested yet either.

## Methods of analysis

There are three tools: response.m, step_analysis.m, and mos.m, for respectively, the response to perturbations, the step time and step length, and the margin of stability.

Each of these three codes is written so they can be easily tested on one trial, by executing them without inputs. You can do this by clicking the "Play" (right triangle) button in the Editor tab. This setup is not yet completed for step_analysis.m and mos.m.

For processing a series of trials, see the "Workflow" section below.

### Response to perturbation

The code in response.m determines the magnitude (T2max) and the duration (t2) of the perturbation response, using the Hotelling T-squared statistic.  (more details and references needed, we may write a paper about this too)

### Step analysis

The code in step_analysis.m determines the right and left step times and step lengths during seconds 10-30 of the trial. Average and standard deviation are calculated of all four variables.

### Margin of stability

The code in mos.m calculates the margin of stability according to the paper by Patricia Young (REFERENCE NEEDED).  (more details needed)  The code is not working yet.

## Workflow for analysis of the trials

This describes the workflow to get all data processed.  It is fully automated, with an option to pause and show details for each trial. This option should be turned on until we know that all trials can be processed automatically.

### Data folders
The trials should be organized into folders, one folder for each test session with up to 15 trials. The folder should contain the mocap files and treadmill files. The only naming requirement is that the filenames should contain a unique trial number that is used to match a mocap file to a treadmill file. 

To set up the processing, create a text file named 0FileList.txt.  In this file, type a list of the mocap files that must be processed.  I have created an example in the Par7_PRE folder.  This allows you to customize the set of files, choose less than 15 trials, choose edited mocap files, etc.

When problems occur during processing, you can still edit the C3D files and make new mocap files using the C3D conversion tools described in the previous section.

### Main program
The program main.m can be set up to process all data, or a subset of the data. Edit the program to select the subject(s) and condition(s) to process. 

This allows you to skip the processing of sessions that were already done.

In this program you can also set the "detail" option.

### Processdata function
The code in processdata.m is called from the main program. It processes all trials from one test session. 

Here you can find settings specific to the analysis, such as the markerset for the response analysis.  The code should be self-explanatory.

The code uses the three analysis tools to generate a number of values from each trial.  These variables are initially stored in a Matlab table, which is then written to an Excel file.  The Excel file is placed in the same session folder where the trials are.  The name of the Excel fils is the same as the folder name.  For instance, after running the processdata function on the Par7_PRE folder, the Par7_PRE folder will contain a file Par7_PRE.xlsx.

### Graphics
Various graphs are generated when "detail" is turned on. The user will be prompted to inspect them for problems, and then hit ENTER to continue.  Some of these graphs may be useful as illustrations in manuscripts or dissertation. You should be able to do Edit->Copy options and Edit-?Copy Figure to export them.  If not, please report the issue.

For each folder, the step analysis results are plotted with trial number on the horizontal axis.  This may show interesting trends, and could be useful to detect fatigue.  This file is also exported as step_analysis.png.

If other graphs, or graphs exported to files, are useful, please ask.
