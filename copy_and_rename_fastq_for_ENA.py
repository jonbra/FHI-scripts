import sys
import os

suffix = sys.argv[1]
# Trenger kanskje ikke bruke denne.
# Evt. benytte meg av sys.argv til andre ting. Som hvilken mappe eller noe.

fullSuffix = "." + suffix

print("Adding " + fullSuffix + " to all filenames!")

# Get the current directory
cwd = os.getcwd()

# Get list of all files in current directory
files = os.listdir(cwd)

# Loop through all files and change filenames
for fileName in files:
	newFilename = filename + fullSuffix
	os.rename(fileName, newFilename)

# Jeg kan legge inn ting her i newFilename. 
