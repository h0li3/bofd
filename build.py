# Desc: BOF sources building script
# Author: hello

import os
import sys
import shutil
import subprocess

PROJ_DIR = os.path.realpath(os.path.dirname(__file__))
OUTPUT_DIR = 'bof'
SRC_DIR = 'src'
INC_DIRS = [SRC_DIR]

NO_COLOR = '\033[0m'
RED_COLOR = '\033[31m'
GREEN_COLOR = '\033[32m'

BOFs = open('bofs.txt', 'r').readlines()

number_bofs_rebuilt = 0
number_of_errors = 0

clang_root = 'c:/Users/Administrator/scoop/apps/msys2/current/clang64'
flags = f'-fplugin=pass.dll -fpass-plugin=pass.dll -mllvm -bl={clang_root}/lib -mllvm -bren -fno-exceptions -D _UNICODE -D UNICODE'


def check_time_changes(bof_name: str, output_name: str) -> bool:
    src_path = os.path.join(PROJ_DIR, SRC_DIR, bof_name)
    out_path = os.path.join(PROJ_DIR, OUTPUT_DIR, output_name)
    common_path = os.path.join(PROJ_DIR, SRC_DIR, 'common')

    if not os.path.exists(out_path):
        return False

    if not os.path.isdir(src_path):
        return False

    for p in os.scandir(src_path):
        if os.path.getmtime(p.path) > os.path.getmtime(out_path):
            return False

    for p in os.scandir(common_path):
        if os.path.getmtime(p.path) > os.path.getmtime(out_path):
            return False

    return True


def build(bof_name) -> bool:
    os.chdir(PROJ_DIR)
    arch = 'x64'
    out_name = f'{bof_name}.{arch}.o'
    clang_bin = clang_root + '/bin/clang++'

    src_path = os.path.join(SRC_DIR, bof_name)
    output_file = os.path.join(OUTPUT_DIR, out_name)
    if os.path.isdir(src_path):
        src_file = os.path.join(src_path, bof_name)
    else:
        src_file = src_path
    if os.path.isfile(src_file + '.cc'):
        src_file += '.cc'
    elif os.path.isfile(src_file + '.cpp'):
        src_file += '.cpp'
    elif os.path.isfile(src_file + '.c'):
        clang_bin = clang_root + '/bin/clang'
        src_file += '.c'

    inc_flags = ''
    if INC_DIRS:
        inc_flags = '-I ' + ' -I'.join(INC_DIRS)

    cmd = f'{clang_bin} {inc_flags} -I"{src_path}" {flags} -c {src_file} -o {os.path.join(PROJ_DIR, output_file)}'

    print(f'Building {bof} ... ', end='')
    #print(cmd)

    if not os.path.isfile(src_file):
        print(f'{RED_COLOR}ERR: BOF src_file of "{bof_name}" dose not exist{NO_COLOR}')
        return False

    # Up to date
    if check_time_changes(bof_name, out_name):
        print(f'{GREEN_COLOR}Updated{NO_COLOR}')
        return False

    result = subprocess.run(cmd, capture_output=True)
    if result.returncode != 0:
        global number_of_errors
        print(RED_COLOR + 'ERR: compiler returned', result.returncode, NO_COLOR)
        print(result.stdout.decode())
        print(result.stderr.decode())
        number_of_errors += 1
        return False

    print(f'{GREEN_COLOR}OK! {output_file}{NO_COLOR}')
    return True

if __name__ == '__main__':
    if not os.path.exists(OUTPUT_DIR):
        print(f'[*] {OUTPUT_DIR} dose not exist, create it.')
        os.mkdir(OUTPUT_DIR)
    for bof in BOFs:
        bof = bof.strip()
        if not bof or bof.startswith('#'):
            continue
        if build(bof):
            number_bofs_rebuilt += 1
        sys.stdout.flush()

    if number_of_errors == 0 and number_bofs_rebuilt == 0:
        print('All BOFs are up to date.')
    else:
        print(f'{number_of_errors} BOFs have error.')
        print(f'{number_bofs_rebuilt} BOFs ware built.')

