#!/usr/bin/env python3
import copy
import glob

import yaml

pipeline_config = yaml.load("""
default:
  tags:
    - k8s
  artifacts:
    expire_in: 1 day

.template_noMPI:
  variables:
    MODE: noMPI
  script:
    - cd Hamiltonians/${HAM_NAME}
    - ./clone_alf.sh
    - cd ALF
    - . configure.sh $MACHINE Devel $MODE $CONFIG_ARGS
    - make
    - cd ../Start
    - ../ALF/Prog/ALF.out

.template_MPI:
  variables:
    MODE: MPI
  script:
    - cd Hamiltonians/${HAM_NAME}
    - ./clone_alf.sh
    - cd ALF
    - . configure.sh $MACHINE Devel $MODE $CONFIG_ARGS
    - make
    - cd ../Start
    - mpiexec -n 4 ../ALF/Prog/ALF.out
""", yaml.Loader)


ENVIRONMENTS = yaml.load("""
Bullseye:
    image: git.physik.uni-wuerzburg.de:25812/alf/alf_docker/pyalf-requirements/bullseye
    variables: {MACHINE: GNU}
Bookworm:
    image: git.physik.uni-wuerzburg.de:25812/alf/alf_docker/pyalf-requirements/bookworm
    variables: {MACHINE: GNU}
Intel21:
    image: git.physik.uni-wuerzburg.de:25812/alf/alf_docker/pyalf-requirements/bullseye-intel
    variables: {MACHINE: INTEL}
Intel-2024.2:
    image: git.physik.uni-wuerzburg.de:25812/alf/alf_docker/pyalf-requirements/bookworm-intel-2024.2
    variables: {MACHINE: INTEL}
IntelLLVM-2024.2:
    image: git.physik.uni-wuerzburg.de:25812/alf/alf_docker/pyalf-requirements/bookworm-intel-2024.2
    variables: {MACHINE: INTELLLVM}
PGI-21-03:
    image: git.physik.uni-wuerzburg.de:25812/alf/alf_docker/pyalf-requirements/bullseye-pgi-21-03
    variables: {MACHINE: PGI}
macGNU:
    tags: ['macos']
    variables: {MACHINE: GNU}
""", yaml.Loader)


def prep_runs(hamiltonian_names, env_name, env_spec):
    for mode in ['noMPI', 'MPI']:
        template_name = '.template_' + mode
        jobname = f'test_{mode}_{env_name}'
        pipeline_config[jobname] = {
            **{'extends': template_name},
            **copy.deepcopy(env_spec),
            'parallel': {'matrix':
                        [{'HAM_NAME': hamiltonian_names}]}
        }

if __name__ == "__main__":
    hamiltonian_names = [i.split('/')[1] for i in glob.glob('Hamiltonians/*/clone_alf.sh')]

    for env_name, env_spec in ENVIRONMENTS.items():
        prep_runs(hamiltonian_names, env_name, env_spec)

    with open('generated-config.yml', 'w', encoding='UTF-8') as f:
        f.write(yaml.dump(pipeline_config))
