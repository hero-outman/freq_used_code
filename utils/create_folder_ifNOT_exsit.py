import os


MYDIR = ''
CHECK_FOLDER = os.path.isdir(MYDIR)

if not CHECK_FOLDER:
    os.makedirs(MYDIR)
    print("created folder : ", MYDIR)

else:
    print(MYDIR, "folder already exists.")