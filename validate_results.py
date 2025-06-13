def validate_exon_matches(output_fasta, positive_fasta):
    from collections import defaultdict

    def read_fasta(fasta_path, strip_class_info=False):
        seqs = defaultdict(str)
        current_id = None
        with open(fasta_path, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    current_id = line[1:].strip()
                    if current_id.endswith("N"):
                        current_id = None  # 跳过 null allele 型别

                    if strip_class_info:
                        current_id = current_id.split()[0]  # 去掉 class_I
                elif current_id:
                    seqs[current_id] += line
        return seqs

    # 读取两个 fasta 文件，output_fasta 需 strip class info
    output_seqs = read_fasta(output_fasta, strip_class_info=True)
    positive_seqs = read_fasta(positive_fasta)

    matched = 0
    missing = []
    mismatched = []

    for allele, pos_seq in positive_seqs.items():
        if allele not in output_seqs:
            missing.append(allele)
        elif output_seqs[allele] != pos_seq:
            mismatched.append(allele)
        else:
            matched += 1
            print(f"✅ 匹配成功: {allele}")

    # 总结输出
    print(f"\n总计阳性型别: {len(positive_seqs)}")
    print(f"  - 匹配成功: {matched}")
    print(f"  - 缺失: {len(missing)}")
    print(f"  - 序列不一致: {len(mismatched)}")

    if missing:
        print("❌ 缺失型别:", ", ".join(missing))
    if mismatched:
        print("⚠️ 序列不一致型别:", ", ".join(mismatched))

if __name__ == "__main__":
    validate_exon_matches("./data/hla_exon_sequences.fasta", "./data/ABC_prot.fa")
