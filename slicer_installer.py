#!/usr/bin/env python3
#pyright: reportUnboundVariable=false


import shutil
import subprocess
import sys
import os


# Function to search the input_file for a line containing old_string, then replace that whole line with the new_string
def replace_path(input_file, old_string, new_string):
    old_file = open(input_file, "r")

    new_file = ""

    for line in old_file.readlines():
        if old_string in line:
            new_line = new_string + "\n"
        else:
            new_line = line
        
        new_file += new_line
        old_file.close()

    writing_file = open(input_file, "w")
    writing_file.write(new_file)
    writing_file.close()


conda_file = shutil.which('conda') # Get the full path to conda executable

if 'anaconda3' in conda_file.split('/'):
    conda_dir = '{}/anaconda3'.format(conda_file.split('/anaconda3/')[0]) # Get the parent directory of anaconda 3

else:
    conda_dir = '/'.join(conda_file.split('/')[0:(conda_file.split('/').index('bin'))])

    if conda_dir + '/bin/conda' != conda_file:   

        while True: # Infinite loop until broken by user action
            user_conda = input('Installer could not identify anaconda parent directory. \
                Your anaconda installation probably has a different parent directory name. \
                To continue, please type in the full path to your anaconda parent directory (e.g., "/home/username/Anaconda"). \
                The more sure way to do this is to call "which conda" in your terminal, and copy and paste the output up to the anaconda parent directory \
                (e.g., if "which conda" output = "/home/username/Anaconda/bin/conda", then paste in "/home/username/Anaconda".')

            if user_conda.split('/') in conda_file.split('/'):
                conda_dir = user_conda
                print('Anaconda directory confirmed. Continuing installation.')
                break

            else:
                print('Anaconda directory is incorrect or not found. Please try again.') # Prompt them to try again

print('Anaconda directory is {}'.format(conda_dir))


while True: # Infinite loop until broken by user action
    user_answer = input('Do you want to install SLICER in a specific environment? (Y/N) ') # Proceed only after the user provided Y or N answer to this question

    if user_answer.lower() == 'y': # If user is installing SLICER in a newly created environment
        user_environment = input('Please type in the name of your SLICER environment: ') # Ask what the user wants the newly created conda environment to be named with
        slicer_env_dir = '{}/envs/{}'.format(conda_dir, user_environment)
        prefix_string = 'prefix: {}'.format(slicer_env_dir) # Path to the directory in which the newly created environment components are installed
        name_string = 'name: {}'.format(user_environment) # The name of the newly created environment
        input('SLICER will be installed in {} environment. Press ENTER to continue'.format(user_environment)) # Notify the user
        break # Dialogue finished. Proceed with the installation.

    if user_answer.lower() == 'n': # If user is installing SLICER in the base environment (where Anaconda 3 is installed)
        user_environment = 'base' # The default name of the environment where Anaconda 3 is installed
        slicer_env_dir = '{}'.format(conda_dir)
        prefix_string = 'prefix: {}'.format(slicer_env_dir) # Path to the directory in which the newly created environment components are installed
        name_string = 'name: {}'.format(user_environment) # The name of the newly created environment
        input('SLICER will be installed in base environment. Press ENTER to continue') # Notify the user
        break # Dialogue finished. Proceed with the installation.

    else: # If user typed in anything else other than Y or N (case-insensitive)
        print('Invalid input. Please try again.') # Prompt them to try again


# shutil.copy('{}/slicer_env_linux.yml'.format(sys.path[0]), '{}/slicer_{}.yml'.format(sys.path[0], user_environment)) # Copy paste the provided slicer_env_linux.yml file
# with open('{}/slicer_{}.yml'.format(sys.path[0], user_environment), 'r') as original_yml: yml_contents = original_yml.readlines() # Read the copy pasted file contents as list of lines

if sys.platform == "linux" or sys.platform == "linux2":
    shutil.copy('{}/slicer_env_linux.yml'.format(sys.path[0]), '{}/slicer_{}.yml'.format(sys.path[0], user_environment)) # Copy paste the provided slicer_env_linux.yml file
    with open('{}/slicer_{}.yml'.format(sys.path[0], user_environment), 'r') as original_yml: yml_contents = original_yml.readlines() # Read the copy pasted file contents as list of lines

elif sys.platform == "darwin":
    subprocess.run('sudo xcodebuild -license accept', shell = True)
    subprocess.run('xcode-select --install', shell = True)
    subprocess.run('sudo chown -R $(whoami) /usr/local || sudo chown -R $(whoami) $(brew --prefix)/*', shell = True)
    # subprocess.run('sudo chown -R $(whoami) /usr/local/share/zsh /usr/local/share/zsh/site-functions', shell = True)
    # subprocess.run('chmod u+w /usr/local/share/zsh /usr/local/share/zsh/site-functions', shell = True)
    subprocess.run('git -C /usr/local/Homebrew/Library/Taps/homebrew/homebrew-core fetch --unshallow || git -C /usr/local/Homebrew/Library/Taps/homebrew/homebrew-core fetch --all', shell = True)
    subprocess.run('git -C /usr/local/Homebrew/Library/Taps/homebrew/homebrew-cask fetch --unshallow || git -C /usr/local/Homebrew/Library/Taps/homebrew/homebrew-cask fetch --all', shell = True)
    subprocess.run('/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"', shell = True)
    # subprocess.run('sudo chown -R $(whoami) /usr/local/var/homebrew', shell = True)
    subprocess.run('brew update', shell = True)
    subprocess.run('brew doctor', shell = True)
    subprocess.run('brew install wget', shell = True)
    subprocess.run('brew upgrade', shell = True)
    subprocess.run('brew install gnu-sed', shell = True)
    shutil.copy('{}/slicer_env_macos.yml'.format(sys.path[0]), '{}/slicer_{}.yml'.format(sys.path[0], user_environment))
    with open('{}/slicer_{}.yml'.format(sys.path[0], user_environment), 'r') as original_yml: yml_contents = original_yml.readlines() # Read the copy pasted file contents as list of lines
    subprocess.run('pip install pandas pysam seaborn pyfastx venn scipy Levenshtein', shell = True)


new_yml_contents = []    
for yml_contents_row in yml_contents:
    if 'name: ' not in yml_contents_row and 'prefix: ' not in yml_contents_row:
        new_yml_contents.append(yml_contents_row) # Append all lines into the list except the line that defines the name and the directory of the environment

with open('{}/slicer_{}.yml'.format(sys.path[0], user_environment), 'w') as modified_yml: 
    modified_yml.write('{}\n'.format(name_string) + ''.join(new_yml_contents) + '{}\n'.format(prefix_string)) # Insert the environment name as the first item in the list and the environment directory (prefix) as the last item in the list

if user_answer.lower() == 'y': # Execute this if user is installing SLICER in a newly created environment
    subprocess.run('conda env create -f {}/slicer_{}.yml'.format(sys.path[0], user_environment), shell = True)

if user_answer.lower() == 'n': # Execute this if user is installing SLICER in the base environment (where Anaconda 3 is installed)
    subprocess.run('conda env update -f {}/slicer_{}.yml'.format(sys.path[0], user_environment), shell = True)

subprocess.run('chmod +x {}/slicer.py'.format(sys.path[0]), shell = True) # Mark all scripts in folder "slicer_scripts" as executable

path_to_slicer_scripts = 'PATH=$PATH:{}'.format(sys.path[0]) # Define the path to the directory "slicer_scripts"

path_exist = 0
# with open(os.path.expanduser('~/.bashrc'), "r+") as file:
#     for line in file:
#         if '{}\n'.format(path_to_slicer_scripts) == line: # Look for the path to the directory "slicer_scripts"
#             path_exist = 1 # If it is already there
#             break # Stop looking and proceed with the next process

#     if path_exist == 0: # If it is not found until the end of the lines
#         file.write('{}\n'.format(path_to_slicer_scripts)) # Write the path to the directory "slicer_scripts" so the scripts within can be accessed from any directory
        
# subprocess.run('exec bash', shell = True) # reload bashrc

if sys.platform == "linux" or sys.platform == "linux2":
    # Linux
    with open(os.path.expanduser('~/.bashrc'), "r+") as file:
        for line in file:
            if '{}\n'.format(path_to_slicer_scripts) == line: # Look for the path to the directory "slicer_scripts"
                path_exist = 1 # If it is already there
                break # Stop looking and proceed with the next process

        if path_exist == 0: # If it is not found until the end of the lines
            file.write('{}\n'.format(path_to_slicer_scripts)) # Write the path to the directory "slicer_scripts" so the scripts within can be accessed from any directory

    subprocess.run('exec bash', shell = True) # Reload bashrc

elif sys.platform == "darwin":
    # OSX

    if os.path.isfile('~/.bash_profile'):
        print (".bash_profile exist")
        with open(os.path.expanduser('~/.bash_profile'), "r+") as file:
            for line in file:
                if '{}\n'.format(path_to_slicer_scripts) == line:
                    break
                else: # Not found, we are at the EOF
                    cwd = os.getcwd()
                    file.write('export PATH="{}/slicer:$PATH"'.format(cwd))
    else:
        print (".bash_profile doesnt exist, making it")
        with open(os.path.expanduser('~/.bash_profile'), "w+") as file:
            cwd = os.getcwd()
            file.write('export PATH="{}/slicer:$PATH"'.format(cwd))
    subprocess.run('source ~/.bash_profile', shell = True) # Reload bash_profile

subprocess.run('rm slicer_{}.yml'.format(user_environment), shell = True) # Housekeeping. Remove auto-generated .yml file of the new environment