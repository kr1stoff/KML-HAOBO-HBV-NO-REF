# 使用 tantan 获取低复杂度区域(Low Complexity Region)
from pathlib import Path
from Bio import SeqIO
import sys

sys.stderr = open(snakemake.log[0], "w")

input_fasta = snakemake.input[0]
output_bed = snakemake.output[0]


def extract_tantan_lcr(input_fasta, output_bed, min_len: int = 6):
    """使用 tantan 获取低复杂度区域(Low Complexity Region)"""
    output_bed = Path(output_bed).resolve()
    # 解析 tantan 输出, 小写字母表示低复杂度区域
    record = list(SeqIO.parse(input_fasta, 'fasta'))[0]
    chrom = record.id
    lower_bases = []
    pos = 0
    for base in record.seq:
        if base.islower():
            lower_bases.append(pos)
        pos += 1

    # 转换为 BED 格式
    ranges = to_ranges(lower_bases)
    with open(output_bed, 'w') as f:
        for start, end in ranges:
            # * 限制低复杂区最小长度
            if end - start >= min_len:
                f.write(f"{chrom}\t{start}\t{end}\n")


def to_ranges(nums):
    """把升序整数列表合并成连续区段字符串"""
    if not nums:
        return []
    nums = sorted(set(nums))          # 去重并排序
    ranges, start = [], nums[0]

    for pre, cur in zip(nums, nums[1:] + [None]):
        if cur is None or cur != pre + 1:      # 连续中断
            if start != pre:                   # 区间
                ranges.append([start, pre])
            start = cur
    return ranges


# main
extract_tantan_lcr(input_fasta, output_bed)
