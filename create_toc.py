#!/usr/bin/env python3
"""Create table of contents _toc.yml out of _toc_template.yml"""

from glob import glob

import yaml

ham_list = glob('Hamiltonians/*')

print(ham_list)

with open('_toc_template.yml', 'r', encoding='UTF-8') as file:
    toc = yaml.safe_load(file)

toc['chapters'] = [{'file': f'{ham_dir}/doc'} for ham_dir in ham_list]

with open('_toc.yml', 'w', encoding='UTF-8') as file:
    yaml.dump(toc, file)
