#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import pandas as pd
import numpy as np
import click
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
GRANTHAM = {(AA_LIST[i], AA_LIST[j]): GRANTHAM_MATRIX[i][j] for i in range(20) for j in range(20)}

# ------------------ 计算函数 ------------------
def grantham_distance(a1, a2):
    return 0 if a1 == a2 else GRANTHAM.get((a1, a2), 0)

def calculate_hed(s1, s2):
    if len(s1) != len(s2):
        raise ValueError(f"Length mismatch: {len(s1)} vs {len(s2)}")
    total = sum(grantham_distance(a, b) for a, b in zip(s1, s2))
    return total / len(s1)

def calculate_hed_trimmed(s1, s2):
    min_len = min(len(s1), len(s2))
    total = sum(grantham_distance(a, b) for a, b in zip(s1[:min_len], s2[:min_len]))
    return total / min_len

def load_allele_sequences(fasta_path):
    seqs = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        allele = record.id.split()[0]
        seqs[allele] = str(record.seq).replace("\n", "").replace(" ", "")
    return seqs

def summarize_hed_per_locus_with_calculation(directory="data", fasta_file="data/hla_exon_sequences.fasta", trim=False):
    target_genes = ["A", "B", "C", "DRB1", "DQB1", "DQA1", "DPB1", "DPA1", "DRB3", "DRB4", "DRB5"]
    allele_seqs = load_allele_sequences(fasta_file)

    print("\n📊 Summary of HED per HLA locus:")
    print(f"{'HLA Locus':<10} {'Trim':<6} {'Median HED':>12} {'IQR':>20} {'Valid Pairs':>15}")

    for gene in target_genes:
        input_path = os.path.join(directory, f"{gene}.txt")
        output_path = os.path.join(directory, f"{gene}_anno.txt")
        output_trim_path = os.path.join(directory, f"{gene}_anno_trim.txt")

        if not os.path.exists(input_path):
            print(f"⚠️ Input file not found: {input_path}")
            continue

        df = pd.read_csv(input_path, sep="\t")
        if df.shape[1] < 3:
            print(f"⚠️ Invalid columns in {input_path}")
            continue

        col1, col2 = df.columns[1], df.columns[2]
        heds, heds_trim = [], []

        for _, row in df.iterrows():
            a1, a2 = row[col1], row[col2]
            seq1, seq2 = allele_seqs.get(a1), allele_seqs.get(a2)

            # 原始 HED
            if seq1 and seq2 and len(seq1) == len(seq2):
                try:
                    heds.append(calculate_hed(seq1, seq2))
                except:
                    heds.append(None)
            else:
                heds.append(None)

            # 裁剪后 HED
            if seq1 and seq2:
                try:
                    heds_trim.append(calculate_hed_trimmed(seq1, seq2))
                except:
                    heds_trim.append(None)
            else:
                heds_trim.append(None)

        # 输出文件（无论是否 trim 参数开启，都写出裁剪列）
        df["HED"] = heds
        df["HED_trim"] = heds_trim
        df.to_csv(output_trim_path if trim else output_path, sep="\t", index=False)

        # 打印统计
        def report(values, label):
            vals = [v for v in values if v is not None]
            if not vals:
                print(f"{gene:<10} {label:<6} {'N/A':>12} {'N/A':>20} {'0':>15}")
            else:
                q1, med, q3 = np.percentile(vals, [25, 50, 75])
                print(f"{gene:<10} {label:<6} {med:12.2f} ({q1:.2f}-{q3:.2f}){len(vals):15d}")

        report(heds, "False")
        report(heds_trim, "True")

# ------------------ CLI 主入口 ------------------
@click.command()
@click.option('--input', '-i', help='输入的 TSV 文件')
@click.option('--fasta', '-f', default="./data/hla_exon_sequences.fasta", help='FASTA文件')
@click.option('--output', '-o', help='输出 TSV 文件')
@click.option('--summary', is_flag=True, help='是否运行summarize')
@click.option('--trim', is_flag=True, help='是否裁剪不同长度序列用于计算HED')
def main(input, fasta, output, summary, trim):
    if summary:
        summarize_hed_per_locus_with_calculation("data", fasta, trim)
    elif input and output:
        print("单文件处理暂不支持 trim 模式，请使用 --summary 模式。")
    else:
        print("请指定 --summary 或 input/output 文件")

if __name__ == "__main__":
    main()
