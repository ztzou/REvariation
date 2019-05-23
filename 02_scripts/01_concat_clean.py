#!/usr/bin/env python

import os
import re
import subprocess

import numpy

in_dir = '02_aligned_cds_pair/'
out_dir = './03_concat_clean/'
record_table = open(f'{out_dir}/record.tsv', 'w')
stops = numpy.array([numpy.fromstring('TGA', dtype='uint8'),
                     numpy.fromstring('TAA', dtype='uint8'),
                     numpy.fromstring('TAG', dtype='uint8')])

def main():
    if not os.path.isdir(out_dir):
        subprocess.call(['mkdir', out_dir])
    dir_list = sorted([
        d for d in os.listdir(in_dir)
        if not os.path.isfile(os.path.join(in_dir, d))
    ])
    for clade in dir_list:
        fasta_list = sorted([
            f for f in os.listdir(f'{in_dir}/{clade}') if f.endswith('fasta')
        ])
        flag = True
        seq_dict = {}
        taxon_list = None
        for i, fasta_file in enumerate(fasta_list):
            taxon_list, cog_seq_dict = read_fasta(
                f'{in_dir}/{clade}/{fasta_file}'
            )
            if flag:
                for taxon in taxon_list:
                    seq_dict[taxon] = cog_seq_dict[taxon]
                flag = False
            else:
                for taxon in taxon_list:
                    seq_dict[taxon] += cog_seq_dict[taxon]
            print(f'{clade}\t{i}')
        alignment_arr = numpy.vstack([
            numpy.fromstring(''.join(seq_dict[taxon]), dtype='uint8')
            for taxon in taxon_list
        ])
        clean_arr = clean_alignment(alignment_arr)
        print(f'{clade}\t{len(taxon_list)}\t{len(fasta_list)}\t'
              f'{alignment_arr.shape[1]}\t{clean_arr.shape[1]}',
              file=record_table)
        with open(f'{out_dir}/{clade}.all.fasta', 'w') as out_f:
            for i, taxon in enumerate(taxon_list):
                print(f'>{taxon}\n{alignment_arr[i, :].tostring().decode()}',
                      file=out_f)
        with open(f'{out_dir}/{clade}.clean.fasta', 'w') as out_f:
            for i, taxon in enumerate(taxon_list):
                print(f'>{taxon}\n{clean_arr[i, :].tostring().decode()}',
                      file=out_f)


def clean_alignment(aln_arr):
    gap_int = numpy.fromstring('-?RYSWKMBDHVN', dtype='uint8')
    arr1 = (aln_arr.reshape(aln_arr.shape[0], -1, 3)[:, :, :, None] == gap_int)
    arr2 = (numpy.sum(arr1.astype('uint8'), axis=(0, 2, 3)) == 0)
    indicator_arr = numpy.vstack((arr2, arr2, arr2)).flatten(order='F')
    clean_arr = aln_arr[:, indicator_arr]
    arr1 = clean_arr.reshape(clean_arr.shape[0], -1, 3)[:, :, None, :] - stops
    arr2 = (numpy.product(arr1.astype('bool').sum(axis=3), axis=(0, 2)) != 0)
    indicator_arr = numpy.vstack((arr2, arr2, arr2)).flatten(order='F')
    clean_arr = clean_arr[:, indicator_arr]
    return clean_arr


def read_fasta(filename):
    with open(filename, 'r') as fasta:
        contents = [x.rstrip() for x in fasta.readlines()]
        taxon_list = [x[1:] for x in contents[::2]]
        cog_seq_dict = dict(zip(taxon_list, contents[1::2]))
    return taxon_list, cog_seq_dict


if __name__ == '__main__':
    main()