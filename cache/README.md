This cache folder stores contains mean and covariance estimates for datasets that have been processed before.
The file names are hash codes and uniquely related to the data.
This will speed up the processing when a file is processed more than once.

The files in this folders should not be committed to the repository.  A .gitignore file is
placed in this folder to prevent this.

The contents of the .gitignore file is:

*     				# all files will be ignored
!.gitignore  		# except .gitignore itself
!README.md			# and the README.md file

