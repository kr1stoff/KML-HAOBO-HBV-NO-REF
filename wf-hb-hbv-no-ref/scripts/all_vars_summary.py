import pandas as pd
import sys
sys.stderr = open(snakemake.log[0], "w")

#  * IO
rawtab = snakemake.input[0]
known_sites_file = snakemake.input[1]
outtab = snakemake.output[0]
outexcel = snakemake.output[1]
# thresholds
core_depth_cutoff = int(snakemake.params.core_depth_cutoff)
other_depth_cutoff = int(snakemake.params.other_depth_cutoff)

# 参考的位点
known_sites = pd.read_csv(known_sites_file, dtype={'Chrom': str, 'Pos': str, 'Ref': str, 'Alt': str})

# 整合深度和频率总得变异字典, key 为 (Chrom, Pos, Ref, Alt)
# 0 depth,  1 freq
vardicts = [{}, {}]
alldf = pd.ExcelFile(rawtab)

for sn in alldf.sheet_names:
    df = pd.read_excel(rawtab, sheet_name=sn, usecols=['Chrom', 'Pos', 'Ref', 'Alt', 'AltDepth', 'AltFreq'],
                       dtype={'Chrom': str, 'Pos': str, 'Ref': str, 'Alt': str})
    df = pd.merge(df, known_sites, on=['Chrom', 'Pos', 'Ref', 'Alt'], how='left')
    # 目标位点和非目标位点分开过滤
    target = df[(~df['Target'].isna()) & (df['AltDepth'] > core_depth_cutoff)]
    nontarget = df[(df['Target'].isna()) & (df['AltDepth'] > other_depth_cutoff)]
    df = pd.concat([target, nontarget], axis=0)
    for it in df.itertuples():
        vardicts[0].setdefault((it.Chrom, it.Pos, it.Ref, it.Alt), {})[sn] = it.AltDepth
        vardicts[1].setdefault((it.Chrom, it.Pos, it.Ref, it.Alt), {})[sn] = it.AltFreq

# ! [20260114] 如果所有样本都没有变异的情况, 会导致 depth_df,freq_df 为空, 后续合并会报错
if vardicts == [{}, {}]:
    outdf = pd.DataFrame(columns=['Chrom', 'Pos', 'Ref', 'Alt', 'Target'])
else:
    # 转成 dataframe, 添加样本信息
    depth_df = pd.DataFrame(vardicts[0]).T.astype('Int64')
    depth_df.columns = [f'{col}-Depth' for col in depth_df.columns]
    freq_df = pd.DataFrame(vardicts[1]).T
    freq_df.columns = [f'{col}-Freq' for col in freq_df.columns]
    # 合并深度和频率表
    sum_df = pd.concat([depth_df, freq_df], axis=1)
    # 排序列名, 加入变异位点信息
    sum_df = sum_df[sorted(sum_df.columns)].reset_index()
    sum_df.columns = ['Chrom', 'Pos', 'Ref', 'Alt'] + list(sum_df.columns[4:])
    # 添加 known_sites
    outdf = pd.merge(known_sites, sum_df, on=['Chrom', 'Pos', 'Ref', 'Alt'], how='outer')
    outdf.drop_duplicates(inplace=True)
    # Pos 改回 int, 并排序
    outdf['Pos'] = outdf['Pos'].astype(int)
    outdf.sort_values(by=['Chrom', 'Pos'], inplace=True)

# 输出
outdf.to_csv(outtab, index=False, sep='\t')
outdf.to_excel(outexcel, index=False)
