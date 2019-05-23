#!/usr/bin/env python

import numpy
import pandas
import scipy.linalg

from CodonZ import std_table
from CodonZ import codonListY as codon_list
from TsTv import isTs

param_dir = '../01_inference/05_re_calculation/filtered/'
out_dir = './00_matrices/'
clade_list = ('03_HsaMmu', '05_Droso')
combin_list = ((0, 0), (1, 1), (0, 1), (1, 0))
d = 0.01
sense_codon_index = [x for x in range(len(codon_list))
                     if not codon_list[x] in ('TAA','TAG','TGA')]
codon_list = [x for x in codon_list if not x in ('TAA','TAG','TGA')]


def main():
    for i1, i2 in combin_list:
        k_clade, w_clade = clade_list[i1], clade_list[i2]
        kappa_df = pandas.read_csv(f'{param_dir}kappa_ws_filtered.tsv',
                                   sep='\t', index_col=0)
        kappa = kappa_df.loc[k_clade]['k']
        cdfreqs_df = pandas.read_csv(f'{param_dir}cdfreqs_filtered.tsv',
                                     sep='\t', index_col=0)
        cdfreqs = cdfreqs_df.loc[k_clade].values[sense_codon_index]
        onew_df = pandas.read_csv(f'{param_dir}onew_filtered.tsv',
                                  sep='\t', index_col=0)
        onew = onew_df.loc[w_clade]['w0']
        re_df = pandas.read_csv(f'{param_dir}RE_mean_filtered.tsv',
                                sep='\t', index_col=0)
        re = re_df.loc[w_clade]
        w = re * onew
        pmat = get_p_matrix(cdfreqs, kappa, w)
        numpy.savetxt(f'{out_dir}/{k_clade}.{w_clade}.mat', pmat, fmt="%.6f")
        numpy.savetxt(f'{out_dir}/{k_clade}.{w_clade}.cdfreq',
                      cdfreqs, fmt="%.6f")
        print(k_clade, w_clade)


def get_p_matrix(pi, kappa, w):
    ncd = len(codon_list)
    qmat = numpy.ones((ncd, ncd))
    pi = pi / numpy.sum(pi)
    for i in range(ncd):
        for j in range(ncd):
            diffL = [
                x for x in range(3) if codon_list[i][x] != codon_list[j][x]
            ]
            if len(diffL) > 1:
                qmat[i, j] = 0
            elif len(diffL) == 1:
                pos = diffL[0]
                if std_table[codon_list[i]] == std_table[codon_list[j]]:
                    if isTs(codon_list[i][pos] + codon_list[j][pos]):
                        qmat[i, j] *= kappa
                else:
                    aai = std_table[codon_list[i]]
                    aaj = std_table[codon_list[j]]
                    if aai + aaj in w.index:
                        omega = w[aai + aaj]
                    else:
                        omega = w[aaj + aai]
                    if isTs(codon_list[i][pos] + codon_list[j][pos]):
                        qmat[i, j] *= kappa * omega
                    else:
                        qmat[i, j] *= omega
                qmat[i, j] *= pi[j]
    for i in range(ncd):
        qmat[i, i] = 0
        qmat[i, i] = - numpy.sum(qmat[i, :])
    scale_f = numpy.sum(pi * numpy.diag(qmat))
    qmat = qmat / numpy.abs(scale_f)
    pmat = scipy.linalg.expm(qmat * d)
    return pmat


if __name__ == '__main__':
    main()