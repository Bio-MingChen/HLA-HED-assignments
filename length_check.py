from collections import defaultdict

def parse_fasta_by_gene(fasta_path):
    """
    读取FASTA文件，返回每个基因中不同长度的型别列表
    返回:
        dict: {gene: {length: [alleles]}}
    """
    gene_length_dict = defaultdict(lambda: defaultdict(list))

    current_id = None
    current_seq = []

    with open(fasta_path, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id and current_seq:
                    seq_str = "".join(current_seq)
                    gene = current_id.split("*")[0]
                    gene_length_dict[gene][len(seq_str)].append(current_id)
                current_id = line[1:].split()[0]  # A*01:01
                current_seq = []
            else:
                current_seq.append(line)

        # 最后一条记录处理
        if current_id and current_seq:
            seq_str = "".join(current_seq)
            gene = current_id.split("*")[0]
            gene_length_dict[gene][len(seq_str)].append(current_id)

    return gene_length_dict

def print_length_distribution(gene_length_dict, output_path="./data/hla_seq_length_stats_by_gene.txt"):
    """
    打印每个基因内各个长度的型别数量，并输出到文件
    """
    with open(output_path, "w", encoding="utf-8") as f:
        for gene in sorted(gene_length_dict.keys()):
            print(f"\n=== 基因: {gene} ===")
            f.write(f"\n=== 基因: {gene} ===\n")
            for length in sorted(gene_length_dict[gene].keys()):
                count = len(gene_length_dict[gene][length])
                print(f"长度: {length}, 型别数: {count}")
                f.write(f"长度: {length}, 型别数: {count}\n")

def get_alleles_by_gene_and_length(gene_length_dict, gene_name, seq_length):
    """
    获取指定基因和指定长度的所有等位基因列表
    """
    return gene_length_dict.get(gene_name, {}).get(seq_length, [])

# 主执行函数
def main():
    fasta_path = "./data/hla_exon_sequences.fasta"
    gene_lengths = parse_fasta_by_gene(fasta_path)
    print_length_distribution(gene_lengths)

    # 示例：获取 C 基因中长度为 273 的等位基因
    target_gene = "A"
    target_length = 177
    alleles = get_alleles_by_gene_and_length(gene_lengths, target_gene, target_length)
    print(f"\n✅ {target_gene} 基因中长度为 {target_length} 的等位基因有 {len(alleles)} 个:")
    for a in alleles:
        print(f" - {a}")

if __name__ == "__main__":
    main()
