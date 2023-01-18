#!/usr/bin/env python
""" sweep.py

    Calling Convention: python3 sweep.py [source file] -key= arg1 arg2 ...

    - source file: Should be runs.jl or the master file from which words will be replaced

    - '-key=': OPTIONAL ARGUMENT that defines the keyword the script will look for in the source
               otherwise, default set to "l_scattering" 

    - arg1, arg2, etc: values to newly map to the keyword; should include units!

    - Example: python3 sweep.py runs.jl -key=l_scattering 10*nm 20*nm 30*nm 

    Output:
        N files correspondings to the N arguments given to the script

        Naming convention: runs_KEYWORD_VALUE.jl
"""

import sys

displayInfo = True
displayError = True
customKeyWord = False

def debug():
    if(displayInfo == True):
        print("Number of arguments passed: ", n)
        print("Name of Python Script: ", sys.argv[0])
        print("Name of Julia Script: ", juliafile)
        if customKeyWord:
            print("Keyword = " + key)
        print('Arguments: ' + str(arguments))

def replace(fileName, argList, keyWord='l_scattering'):
    text = []
    files = []

    with open(fileName, 'r') as reader:
        text = reader.readlines()
        # print(text)
        
        # text = [line.replace('l_scattering = 10*nm', 'l_scattering = 10*nm') for line in text]

        for arg in argList:
            ######### OLD IMPLEMENTATION #########
            # files.append([line.replace(keyWord, 'l_scattering = ' + arg + '*nm') for line in text])

            file = []
            for line in text:
                if keyWord in line:
                    size = len(line) - len(line.lstrip())
                    line = line[0:size] + keyWord + ' = ' + arg + ',\n'
                file.append(line)
            files.append(file)


    for i, file in enumerate(files):       
        with open('runs_' + keyWord.upper() + '_' + argList[i] + '.jl', 'w') as writer:
            for line in file:
                writer.write(line)



n = len(sys.argv)

if n == 2 and sys.argv[1] == '--help':
    print(__doc__)
    sys.exit()


if n == 1 and displayError:
    sys.exit("No Julia file or arguments have been passed! Type 'python3 sweep.py --help' for documentation.")


juliafile = sys.argv[1]

if n == 2 and displayError:
    sys.exit("No arguments have been passed! Type 'python3 sweep.py --help' for documentation.")

# Pass in arguments as space seperated strings! EX: ... runs.jl 10 20 30 40
arguments = []
key = sys.argv[2]
if '-key=' in key:
    key = key.split('=')[1]
    arguments = sys.argv[3:n]
    customKeyWord = True
else:
    arguments = sys.argv[2:n]
    key = None

if displayInfo:
    debug()

# python3 -this file -.jl files -params to sweep - start value -end value - Nsteps (subdivisions) 
# OR
# ... -sweep.py ... --string 1 -string 2 -->
if key != None: 
    replace(fileName=juliafile, keyWord=key, argList=arguments)
else:
    replace(fileName=juliafile, argList=arguments)

            