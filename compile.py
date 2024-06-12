#!/usr/bin/env python3

import os
from pathlib import Path


NAMES = [
    'bga',
    'cloneIGA',
    'gen2phen',
    'iga',
    'pga',
    'PGAinternal'
]

current = Path('.')

for name in NAMES:
    _dir = f"bios.{name}"
    os.chdir(_dir)
    print(os.getcwd())
    os.system('cmake -G Xcode .')
    os.system('cmake --build .')
    os.chdir("..")
