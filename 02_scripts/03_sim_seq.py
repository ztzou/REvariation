#!/usr/bin/env python

import sys

import numpy
import scipy.stats

from CodonZ import codonListY as codon_list

mat_dir = './00_matrices/'
out_dir_all = './01_sim_seq/'
clade_list = ('03_HsaMmu', '05_Droso')
combin_list = ((0, 0), (1, 1), (0, 1), (1, 0))
d = 0.01
codon_count = 5000000
t_div = 0.05
codon_list = [x for x in codon_list if not x in ('TAA','TAG','TGA')]


def main():
    i1, i2 = combin_list[int(sys.argv[1])]
    rep_index = int(sys.argv[2])
    k_clade, w_clade = clade_list[i1], clade_list[i2]
    p_matrix = numpy.loadtxt(f'{mat_dir}{k_clade}.{w_clade}.mat')
    codon_freq = numpy.loadtxt(f'{mat_dir}{k_clade}.{w_clade}.cdfreq')
    seq1, seq2 = sim_seq(p_matrix, codon_freq, t_div / d, codon_count)
    with open(f'{out_dir_all}{k_clade}.{w_clade}'
              f'.rep{rep_index:02d}.fasta', 'w') as f:
        print(f'>taxa1\n{seq1}\n>taxa2\n{seq2}', file=f)


def sim_seq(p_matrix, equi_freq, t, ncd):
    anc_arr = get_seq_by_multinomial(equi_freq, ncd)
    child_arr = anc_arr.copy()
    for codon_index in range(len(codon_list)):
        positions = numpy.arange(child_arr.shape[0])[anc_arr == codon_index]
        position_count = len(positions)
        codon_child_states = get_seq_by_multinomial(
            single_branch_evo(codon_index, t, p_matrix), position_count
        )
        child_arr[positions] = codon_child_states
    seq1 = ''.join([codon_list[x] for x in child_arr])
    child_arr = anc_arr.copy()
    for codon_index in range(len(codon_list)):
        positions = numpy.arange(child_arr.shape[0])[anc_arr == codon_index]
        position_count = len(positions)
        codon_child_states = get_seq_by_multinomial(
            single_branch_evo(codon_index, t, p_matrix), position_count
        )
        child_arr[positions] = codon_child_states
    seq2 = ''.join([codon_list[x] for x in child_arr])
    return seq1, seq2


def get_seq_by_multinomial(equi_freq, site_count):
    equi_freq = equi_freq / numpy.sum(equi_freq)
    cus = scipy.stats.rv_discrete(values=(range(len(codon_list)), equi_freq))
    return numpy.array(cus.rvs(size=site_count), dtype='uint8')


def single_branch_evo(state_index, t, p_matrix):
    ind_vec = (numpy.arange(len(codon_list)) == state_index).astype('float')
    # for i in range(int(t)):
    # ind_vec = numpy.dot(ind_vec, mat)
    ind_vec = numpy.dot(ind_vec, numpy.linalg.matrix_power(p_matrix, int(t)))
    return ind_vec


if __name__ == '__main__':
    main()