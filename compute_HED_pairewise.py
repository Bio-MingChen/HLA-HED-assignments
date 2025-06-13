#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from Bio import SeqIO

# ------------------ Grantham Distance Matrix ------------------
AA_LIST = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
           'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

GRANTHAM_MATRIX = [
    [0,112,111,126,195, 91,107, 60, 86, 94, 96,106, 84,113, 27, 99, 58,148,112, 64],
    [112,  0, 86, 96,180, 43, 54,125, 29, 97,102, 26, 91, 97,103,110, 71,101, 77, 96],
    [111, 86,  0, 23,139, 46, 42, 80, 68,149,153, 94,142,160, 91, 46, 65,174,143,133],
    [126, 96, 23,  0,154, 61, 45, 94, 81,168,172,101,160,177,108, 65, 85,181,160,152],
    [195,180,139,154,  0,154,170,159,174,198,198,202,196,205,169,112,149,215,194,192],
    [ 91, 43, 46, 61,154,  0, 29, 87, 24,109,113, 53,101,116, 76, 68, 42,130, 99, 96],
    [107, 54, 42, 45,170, 29,  0, 98, 40,134,138, 56,126,152, 93, 80, 65,181,140,121],
    [ 60,125, 80, 94,159, 87, 98,  0, 98,135,138,127,127,153, 42, 56, 59,184,147,109],
    [ 86, 29, 68, 81,174, 24, 40, 98,  0, 94, 99, 32, 87,100, 77, 89, 47,115, 83, 84],
    [ 94, 97,149,168,198,109,134,135, 94,  0,  5,102, 10, 21, 95,142, 89, 61, 33, 29],
    [ 96,102,153,172,198,113,138,138, 99,  5,  0,107, 21, 22, 98,145, 92, 61, 36, 32],
    [106, 26, 94,101,202, 53, 56,127, 32,102,107,  0, 95,103,103,121, 78,110, 85, 97],
    [ 84, 91,142,160,196,101,126,127, 87, 10, 21, 95,  0, 28, 87,135, 81, 84, 36, 21],
    [113, 97,160,177,205,116,152,153,100, 21, 22,103, 28,  0,114,155,103, 40, 22, 50],
    [ 27,103, 91,108,169, 76, 93, 42, 77, 95, 98,103, 87,114,  0, 74, 38,147,110, 68],
    [ 99,110, 46, 65,112, 68, 80, 56, 89,142,145,121,135,155, 74,  0, 58,177,144,124],
    [ 58, 71, 65, 85,149, 42, 65, 59, 47, 89, 92, 78, 81,103, 38, 58,  0,128, 92, 69],
    [148,101,174,181,215,130,181,184,115, 61, 61,110, 84, 40,147,177,128,  0, 37, 88],
    [112, 77,143,160,194, 99,140,147, 83, 33, 36, 85, 36, 22,110,144, 92, 37,  0, 55],
    [ 64, 96,133,152,192, 96,121,109, 84, 29, 32, 97, 21, 50, 68,124, 69, 88, 55,  0]
]

GRANTHAM = {
    (AA_LIST[i], AA_LIST[j]): GRANTHAM_MATRIX[i][j]
    for i in range(20) for j in range(20)
}

# ------------------ HED计算函数 ------------------
def grantham_distance(a1, a2):
    return 0 if a1 == a2 else GRANTHAM.get((a1, a2), 0)

# def calculate_hed(s1, s2):
#     if len(s1) != len(s2):
#         raise ValueError(f"Length mismatch: {len(s1)} vs {len(s2)}")
#     return sum(grantham_distance(a, b) for a, b in zip(s1, s2)) / 2.2

def calculate_hed(s1, s2):
    if len(s1) != len(s2):
        raise ValueError(f"Length mismatch: {len(s1)} vs {len(s2)}")
    total = sum(grantham_distance(a, b) for a, b in zip(s1, s2))
    return total / len(s1)  # 注意这里是除以长度，不是 2.2

# ------------------ 序列读取函数 ------------------
def load_allele_sequences(fasta_path):
    seqs = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        allele = record.id.split()[0]
        seq = str(record.seq).replace("\n", "").replace(" ", "")
        seqs[allele] = seq
    return seqs

# ------------------ 主计算函数 ------------------
def compute_hed_between_alleles(allele1, allele2, fasta_path):
    allele_seqs = load_allele_sequences(fasta_path)

    if allele1 not in allele_seqs:
        raise ValueError(f"{allele1} not found in fasta file")
    if allele2 not in allele_seqs:
        raise ValueError(f"{allele2} not found in fasta file")

    s1 = allele_seqs[allele1]
    s2 = allele_seqs[allele2]

    if len(s1) != len(s2):
        raise ValueError(f"Length mismatch between {allele1} and {allele2}")

    hed = calculate_hed(s1, s2)
    print(f"HED({allele1}, {allele2}) = {hed:.2f}")
    return hed

# ------------------ 命令行调用支持 ------------------
def main():
    if len(sys.argv) != 3:
        print("用法: python compute_hed_pairwise.py A*01:01 A*02:01")
        sys.exit(1)

    allele1 = sys.argv[1]
    allele2 = sys.argv[2]
    fasta_path = "./data/hla_exon_sequences.fasta"
    compute_hed_between_alleles(allele1, allele2, fasta_path)

if __name__ == "__main__":
    main()
