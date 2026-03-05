import re
import sys

sys.stderr = open(snakemake.log[0], "w")


# * IO
in_sam = snakemake.input['sam']
fai_file = snakemake.input['fai']
mapped_sam = snakemake.output['mapped']
out_bed = snakemake.output['bed']

# 获取比对到的 SAM 记录
with open(in_sam) as fin, open(mapped_sam, 'w') as fout:
    for line in fin:
        if line.startswith('@'):
            continue
        cols = line.rstrip('\n').split('\t')
        if len(cols) < 2:
            continue
        try:
            flag = int(cols[1])
        except:
            continue
        if (flag & 4) == 0:
            fout.write(line)

# 解析比对结果，生成 BED 文件
ref_lengths = {}
with open(fai_file) as f:
    for line in f:
        p = line.strip().split('\t')
        if len(p) >= 2:
            ref_lengths[p[0]] = int(p[1])
mappings = []
counts = {}
with open(mapped_sam) as f:
    for line in f:
        if line.startswith('@'):
            continue
        cols = line.rstrip('\n').split('\t')
        if len(cols) < 11:
            continue
        qname = cols[0]
        try:
            flag = int(cols[1])
        except:
            continue
        rname = cols[2]
        pos = int(cols[3]) if cols[3].isdigit() else 0
        cigar = cols[5]
        # sum M ops
        m_total = 0
        for m in re.findall(r'(\d+)M', cigar):
            try:
                m_total += int(m)
            except:
                pass
        if m_total == 0 and len(cols) >= 10:
            seq = cols[9]
            m_total = len(seq.replace('*', ''))
        start0 = pos - 1
        end = start0 + m_total
        strand = '-' if (flag & 16) else '+'
        mappings.append((qname, rname, start0, end, strand, cigar))
        counts[qname] = counts.get(qname, 0) + 1
with open(out_bed, 'w') as out:
    for q, r, s, e, st, c in mappings:
        if r not in ref_lengths:
            sys.stderr.write(f"警告: 参考中未找到 contig '{r}'\n")
        else:
            if e > ref_lengths[r]:
                sys.stderr.write(f"警告: 引物 {q} 在 {r} 的 end={e} 超过长度 {ref_lengths[r]}\n")
        out.write(f"{r}\t{s}\t{e}\t{q}\t0\t{st}\n")
sys.stderr.write(f"生成 BED: {out_bed} (records={len(mappings)})")

# 报告多位点映射的引物
multi = {k: v for k, v in counts.items() if v > 1}
if multi:
    sys.stderr.write("注意: 以下引物有多位点映射：")
    for k, v in multi.items():
        sys.stderr.write(f"  {k}: {v} 次")
