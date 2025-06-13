import re
from Bio.Seq import Seq

def parse_hla_dat_4res(dat_path, top=None):
    """
    解析 HLA .dat 文件，提取 HLA Class I 的 exon2+3 和 Class II 的 exon2 的氨基酸序列。

    参数：
        dat_path (str): .dat 文件路径
        top (int): 限制解析的记录数，用于调试

    返回：
        dict: {allele4: {"class": "I" or "II", "seq": translated_aa_sequence}}
    """
    seen_allele4 = set()

    # 支持的经典HLA基因
    classical_genes = {"A", "B", "C", "DRB1", "DRB3", "DRB4", "DRB5", "DQA1", "DQB1", "DPA1", "DPB1"}
    mapping = {}

    # 读取并按 record 切分
    with open(dat_path, encoding='utf-8') as f:
        text = f.read()
    records = re.split(r"(?m)^\s*//\s*$", text)

    count = 0
    for rec in records:
        lines = rec.strip().splitlines()
        if not lines or not lines[0].startswith("ID"):
            continue

        # DE 行抓 allele name 和 class
        allele_full = cls = None
        for line in lines:
            if line.startswith("DE") and "Class I" in line:
                cls = "I"
            elif line.startswith("DE") and "Class II" in line:
                cls = "II"
            if line.startswith("DE") and "HLA-" in line:
                match = re.search(r"HLA-([\w\*:\-]+)", line)
                if match:
                    allele_full = match.group(1)
                    break
        if not allele_full or cls is None:
            print("跳过: 无法识别 DE 行")
            continue

        # 跳过 null allele
        if allele_full.endswith("N"):
            print(f"跳过 null allele: {allele_full}")
            continue

        # 只保留经典基因
        locus = allele_full.split("*")[0]
        if locus not in classical_genes:
            print(f"跳过非经典基因: {allele_full}")
            continue

        # 构造 4-digit allele
        fields = allele_full.split("*")[1].split(":")
        if len(fields) < 2:
            print(f"跳过格式错误: {allele_full}")
            continue
        allele4 = f"{locus}*{fields[0]}:{fields[1]}"
        if allele4 in seen_allele4:
            print(f"跳过重复 4-digit allele: {allele4}")
            continue
        print(f"处理 allele: {allele4} (class {cls})")

        # 提取 /translation（多行拼接）和 /codon_start
        translation_aa = ""
        in_translation = False
        codon_start = 1  # 默认值

        for line in lines:
            if '/translation="' in line:
                in_translation = True
                part = line.split('="')[1].strip()
                translation_aa += part.rstrip('"')
                if part.endswith('"'):
                    in_translation = False
            elif in_translation:
                # 去除行首的 "FT" 和空格，保留实际氨基酸序列
                aa_part = line.strip()
                if aa_part.startswith("FT"):
                    aa_part = aa_part[2:].strip()
                translation_aa += aa_part.rstrip('"')
                if line.strip().endswith('"'):
                    in_translation = False
            elif line.strip().startswith('/codon_start='):
                codon_start = int(line.strip().split('=')[1])

        # 最终 strip 一下，保险起见
        translation_aa = translation_aa.strip() if translation_aa else None

        # 获取 CDS 区域坐标（拼接所有 CDS 相关的 FT 行）
        cds_txt_lines = []
        in_cds = False
        for line in lines:
            if line.startswith("FT   CDS"):
                in_cds = True
                cds_txt_lines.append(line)
            elif in_cds and line.startswith("FT"):
                cds_txt_lines.append(line)
            elif in_cds:
                break
        cds_txt = "\n".join(cds_txt_lines)

        coords = re.findall(r"(\d+)\.\.(\d+)", cds_txt)
        if not coords:
            print(f"跳过无 CDS: {allele_full}")
            continue
        cds_coords = [(int(s) - 1, int(e)) for s, e in coords]  # 注意减1：python index从0开始
        print(f"CDS 坐标: {cds_coords}")

        # 获取 exon 坐标（用于 later exon 切片）
        exon_coords = []
        for line in lines:
            if line.startswith("FT   exon"):
                ex_match = re.search(r"(\d+)\.\.(\d+)", line)
                if ex_match:
                    start = int(ex_match.group(1)) - 1  # Python 0-based
                    end = int(ex_match.group(2))
                    exon_coords.append((start, end))

        print(f"exon 区段: {exon_coords}")

        # 获取 DNA 序列
        in_seq = False
        dna_lines = []
        for line in lines:
            if line.startswith("SQ"):
                in_seq = True
                continue
            if in_seq:
                dna_lines.append(re.sub(r"[^ACGTacgt]", "", line))
        dna = "".join(dna_lines).upper()
        if not dna:
            print(f"跳过: 无 DNA 序列")
            continue

        # 提取 CDS 段落并处理 codon_start 偏移（只对第一段偏移）
        cds_segs = []
        for idx, (s, e) in enumerate(cds_coords):
            seg = dna[s:e]
            if idx == 0:
                seg = seg[codon_start - 1:]  # 只对第一段偏移
            cds_segs.append(seg)

        cds_full = "".join(cds_segs)
        cds_full = cds_full[:len(cds_full)//3*3]


        aa_full = str(Seq(cds_full).translate(table=1))
        print(f"翻译长度: {len(aa_full)} 氨基酸")

        # 定位 exon 片段
        # 提取 codon_start (前面已有，不变)
        codon_start = 1  # 默认值
        for line in lines:
            if '/codon_start=' in line:
                codon_start = int(line.strip().split('=')[1])

        # 提取 exon number 和坐标
        exon_number_coords = {}
        for i, line in enumerate(lines):
            if line.strip().startswith("FT   exon"):
                exon_coord_match = re.search(r'(\d+)\.\.(\d+)', line)
                exon_num_match = re.search(r'/number="(\d+)"', lines[i+1] if i+1 < len(lines) else "")
                if exon_coord_match and exon_num_match:
                    exon_number = int(exon_num_match.group(1))
                    start, end = int(exon_coord_match.group(1)) - 1, int(exon_coord_match.group(2))
                    exon_number_coords[exon_number] = (start, end)

        # 根据 exon_number 获取 exon1 (用于记录), exon2, exon3
        exon1 = exon_number_coords.get(1, None)
        exon2 = exon_number_coords.get(2, None)
        exon3 = exon_number_coords.get(3, None)

        # 检查 exon2 是否存在（必需）
        if exon2 is None:
            print(f"跳过: exon2 不存在")
            continue

        # 拼接分析用 exon（只用于 exon2 + exon3）
        analysis_exon_segs = []
        analysis_exon_segs.append(dna[exon2[0]:exon2[1]])
        if cls == "I":
            if exon3:
                analysis_exon_segs.append(dna[exon3[0]:exon3[1]])
            else:
                print(f"跳过: class I exon3 不存在")
                continue

        exon_dna = "".join(analysis_exon_segs)

        # codon_start 偏移仅当 exon1 缺失时使用
        if exon1 is None and codon_start > 1:
            exon_dna = exon_dna[codon_start - 1:]

        # 确保长度为3的倍数
        exon_dna = exon_dna[:len(exon_dna)//3*3]

        # 翻译为氨基酸（分析+验证用）
        print(f"exon_dna is {exon_dna}")

        # 尝试3个reading frame，找最匹配translation_aa的那个
        best_match = ""
        best_offset = -1

        for offset in range(3):
            frame_dna = exon_dna[offset:]
            frame_dna = frame_dna[:len(frame_dna) // 3 * 3]
            frame_aa = str(Seq(frame_dna).translate(table=1))
            if frame_aa in translation_aa:
                best_match = frame_aa
                best_offset = offset
                break  # 完全匹配直接选用
            if len(frame_aa) > len(best_match):
                best_match = frame_aa
                best_offset = offset

        exon_aa = best_match.split("*", 1)[0]
        print(f"✅ 使用 reading frame {best_offset} 翻译")
        print(f"exon_aa is {exon_aa}")

        # 与 translation 注释比对 (严格版本)
        MIN_EXON_AA_LENGTH = 50  # 定义最小合理长度

        if translation_aa:
            if len(exon_aa) < MIN_EXON_AA_LENGTH:
                print(f"跳过: exon氨基酸序列过短 ({len(exon_aa)} 个氨基酸)")
                continue
            if exon_aa not in translation_aa:
                print(f"跳过: exon翻译与 /translation 不匹配")
                print(f"[Translation片段开头] {translation_aa[:60]}...")
                print(f"[Exon翻译片段]      {exon_aa[:60]}...")
                print(f"[Translation长度] {len(translation_aa)}")
                print(f"[Exon翻译长度]    {len(exon_aa)}")
                continue
            else:
                print(f"✅ exon翻译在 /translation 中匹配成功 (长度 {len(exon_aa)} 个氨基酸)")

        # 最终存储（直接使用 exon_aa）
        mapping[allele4] = {"class": cls, "seq": exon_aa}


        seen_allele4.add(allele4)
        count += 1
        print(f"当前已成功解析记录数：{count}")
        if top is not None and count >= top:
            print(f"已达到 top={top} 限制，停止解析。")
            break

    return mapping


def write_mapping_to_fasta(mapping, output_path):
    """
    将解析得到的HLA序列mapping字典写入fasta文件。
    每个allele4为一个记录，记录中包含class和对应exon氨基酸序列。

    参数：
        mapping (dict): {allele4: {"class": "I"/"II", "seq": translated_aa_sequence}}
        output_path (str): 输出 fasta 文件路径
    """
    with open(output_path, 'w') as f:
        for allele4, info in mapping.items():
            seq = info['seq']
            cls = info['class']
            f.write(f">{allele4} class_{cls}\n")
            # fasta 推荐每行不超过60个字符，规范输出
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + '\n')

def main():
    dat_file = "./data/hla.dat"
    # dat_file = "./data/test1.dat"
    mapping = parse_hla_dat_4res(dat_file, top=None)
    print("\nParsed sequences:", list(mapping.keys()))

    # 新增的输出调用
    output_fasta = "./data/hla_exon_sequences.fasta"
    print(mapping)
    write_mapping_to_fasta(mapping, output_fasta)
    print(f"成功将结果写入到 {output_fasta}")

if __name__ == "__main__":
    main()
