import os
import sys
import shutil
import subprocess

PROJ_DIR = os.path.realpath(os.path.dirname(__file__))
OUTPUT_DIR = os.path.join(PROJ_DIR, 'bof')
SRC_DIR = os.path.join(PROJ_DIR, 'src')

BOFs = open('bofs.txt', 'r').readlines()

def build(bof_name):
    src_path = os.path.join(SRC_DIR, bof_name)
    if os.path.isdir(src_path):
        src_path = os.path.join(src_path, bof_name + '.c')
    else:
        src_path += '.c'

    if not os.path.isfile(src_path):
        print(f'[-] "{src_path}" dose not exist')
        return
    result = subprocess.run('cl.exe /GS- -c ' + src_path, capture_output=True)
    if result.returncode != 0:
        print('[-] could not build bof:', result.returncode)
        return
    output_path = bof_name + '.obj'
    bof_path = os.path.join(OUTPUT_DIR, bof_name + '.o')
    shutil.move(output_path, bof_path)

if __name__ == '__main__':
    if not os.path.exists(OUTPUT_DIR):
        print(f'[*] {OUTPUT_DIR} dose not exist, build it')
        os.mkdir(OUTPUT_DIR)
    for bof in BOFs:
        bof = bof.strip()
        print(f'[*] building bof "{bof}"')
        build(bof.strip())

