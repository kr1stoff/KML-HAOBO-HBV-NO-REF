import vcfpy
from collections import defaultdict
import sys

sys.stderr = open(snakemake.log[0], "w")


in_vcf = snakemake.input[0]
out_tab = snakemake.output[0]

reader = vcfpy.Reader.from_path(in_vcf)
f = open(out_tab, "w")
f.write('\t'.join(['Chrom', 'Pos', 'Ref', 'Alt', 'TotalDepth', 'AltDepth', 'AltFreq', 'SAF', 'SAR',
        'ForwardStrandRate', 'ReverseStrandRate', 'RPL', 'RPR', 'PlacedLeftRate', 'PlacedRightRate']) + '\n')

outdict = defaultdict(lambda: {
    'dp': 0,
    'ao': [],
    'saf': [],
    'sar': [],
    'rpl': [],
    'rpr': [],
})

for record in reader:
    alt = record.ALT[0].value
    key = (record.CHROM, record.POS, record.REF, alt)
    outdict[key]['dp'] = record.INFO.get('DP')
    # 由于alt可能是多个所以值为列表，当前项目中经过 bcftools norm + vt decompose_blocksub 后 alt只有一个
    outdict[key]['ao'].append(record.INFO.get('AO', [0])[0])
    outdict[key]['saf'].append(record.INFO.get('SAF', [0])[0])
    outdict[key]['sar'].append(record.INFO.get('SAR', [0])[0])
    outdict[key]['rpl'].append(record.INFO.get('RPL', [0])[0])
    outdict[key]['rpr'].append(record.INFO.get('RPR', [0])[0])

for key in outdict:
    ao = sum(outdict[key]['ao'])
    saf = sum(outdict[key]['saf'])
    sar = sum(outdict[key]['sar'])
    rpl = sum(outdict[key]['rpl'])
    rpr = sum(outdict[key]['rpr'])
    dp = outdict[key]['dp']
    # 计算等位基因频率，链偏倚数值，位置偏倚数值
    af = round(ao / dp if dp > 0 else 0, 4)
    sa_sum = saf + sar
    saf_rate = round(saf / sa_sum if sa_sum > 0 else 0, 4)
    sar_rate = round(sar / sa_sum if sa_sum > 0 else 0, 4)
    rpl_rate = round(rpl / (rpl + rpr) if rpl + rpr > 0 else 0, 4)
    rpr_rate = round(rpr / (rpl + rpr) if rpl + rpr > 0 else 0, 4)

    f.write('\t'.join(
        map(str, list(key) + [outdict[key]['dp'], ao, af, saf, sar, saf_rate, sar_rate, rpl, rpr, rpl_rate, rpr_rate])) + '\n')

f.close()
