If running on a lab machine, you should not have to do any additional setup
as the needed libraries are installed.  If you are running on other machines,
you'll want to create a conda environment to run the program in and install 
the necessary libraries.  Once the conda environment is created, activate it 
and run the following command in the directory where you unpacked these files:

make conda-setup-run

This will run a conda command to install the necessary python libraries for
the tools to run.

To invoke the program, you currently have to run the script in the directory
where you have unpacked the files.  The data can be anywhere.  The program
is invoked with the following command:

python PSFfitting.py <tinytim file> <object image file>

