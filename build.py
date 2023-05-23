import os
import sys
import shutil
import subprocess

PROJ_DIR = os.path.realpath(os.path.dirname(__file__))
OUTPUT_DIR = 'bof'
SRC_DIR = 'src'
INC_DIR = 'common'

NO_COLOR = '\033[0m'
RED_COLOR = '\033[31m'
GREEN_COLOR = '\033[32m'

BOFs = open('bofs.txt', 'r').readlines()

number_bofs_rebuilt = 0
number_of_errors = 0

def build(bof_name) -> bool:
    os.chdir(PROJ_DIR)

    src_path = os.path.join(SRC_DIR, bof_name)
    output_file = os.path.join(OUTPUT_DIR, f'{bof_name}.o')
    if os.path.isdir(src_path):
        src_file = os.path.join(src_path, bof_name + '.c')
    else:
        src_file = src_path + '.c'
    cmd = f'gcc.exe -O2 -I "{INC_DIR}" -I "{src_path}" -c {src_file} -o {output_file}'

    # Up to date
    if os.path.exists(output_file) and os.path.getmtime(src_file) < os.path.getmtime(output_file):
        return False

    print(f'Building {bof} ... ', end='')

    if not os.path.isfile(src_file):
        print(f'{RED_COLOR}ERR: BOF src_file of "{bof_name}" dose not exist{NO_COLOR}')
        return False

    result = subprocess.run(cmd, capture_output=True)
    if result.returncode != 0:
        global number_of_errors
        print(RED_COLOR + 'ERR: compiler returned', result.returncode, NO_COLOR)
        print(result.stdout.decode())
        print(result.stderr.decode())
        number_of_errors += 1
        return False

    print(GREEN_COLOR + 'OK!' + NO_COLOR)
    return True

if __name__ == '__main__':
    if not os.path.exists(OUTPUT_DIR):
        print(f'[*] {OUTPUT_DIR} dose not exist, create it.')
        os.mkdir(OUTPUT_DIR)
    for bof in BOFs:
        bof = bof.strip()
        sys.stdout.flush()
        if build(bof.strip()):
            number_bofs_rebuilt += 1

    if number_of_errors == 0 and number_bofs_rebuilt == 0:
        print('All BOFs are up to date.')
    else:
        print(f'{number_of_errors} BOFs have error.')
        print(f'{number_bofs_rebuilt} BOFs ware built.')

