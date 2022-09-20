#!/usr/bin/env python
import sys


def read_score(dock_pdb_file, vina_v='vina'):
    fp = open(dock_pdb_file)
    lines = fp.readlines()
    fp.close()
    score_list = list()
    for line in lines:
        if vina_v == 'vina':
            if line.startswith('REMARK VINA RESULT'):
                score = float(line[19:29].strip())
                score_list += [score]
        elif vina_v == 'smina':
            if line.startswith('REMARK minimizedAffinity'):
                score = float(line[24:36].strip())
                score_list += [score]

    return score_list


def main():

    if len(sys.argv) < 3:
        print('gather_score.py smi_list_file out_dir vina_v')
        sys.exit()

    list_file = sys.argv[1]
    out_dir = sys.argv[2]
    vina_v = 'vina'
    if len(sys.argv) > 3:
        if sys.argv[3].lower() == 'smina':  # vina or smina
            vina_v = 'smina'

    fp = open(list_file)
    lines = fp.readlines()
    fp.close()
    for line in lines:
        lis = line.strip().split()
        mol_id = lis[0]
        dock_pdb_file = out_dir + '/' + '/dock_' + mol_id + '.pdb'

        score = read_score(dock_pdb_file, vina_v)
        print(mol_id, score)


if __name__ == "__main__":
    main()
