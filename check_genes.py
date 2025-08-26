#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
check_genes.py – chuyển ENSG→Symbol và kiểm tra gene có trong OncoKB (DrugTarget = Yes)
– chấp nhận cả cột Gene Aliases.
"""

import pandas as pd, mygene, argparse, re

# ---------- 1. Tham số dòng lệnh ----------
parser = argparse.ArgumentParser(description="Check OncoKB (Hugo Symbol + Aliases)")
parser.add_argument('--input',  required=True, help='File .txt kết quả BRW')
parser.add_argument('--oncokb', required=True, help='File Excel OncoKB (có cột Gene Aliases)')
parser.add_argument('--output', default='top100_checked.tsv', help='File TSV đầu ra')
args = parser.parse_args()

# ---------- 2. Đọc Top-100 ENSG ----------
df_res = (pd.read_csv(args.input, sep='\t', nrows=100)
            .rename(columns={'GeneNames': 'Ensembl_ID'}))

# ---------- 3. Tra cứu ENSG → Symbol ----------
mg = mygene.MyGeneInfo()
map_df = (pd.DataFrame(
            mg.querymany(df_res['Ensembl_ID'].tolist(),
                         scopes='ensembl.gene',
                         fields='symbol',
                         species='human'))
          .rename(columns={'query':  'Ensembl_ID',
                           'symbol': 'Symbol'})
          .dropna(subset=['Symbol']))

# ---------- 4. Đọc OncoKB (chỉ Drug target = Yes) ----------
okb = (pd.read_excel(args.oncokb,
                     usecols=['Hugo Symbol', 'Cancer Drug target gene', 'Gene Aliases'])
         .rename(columns={'Hugo Symbol': 'Symbol',
                          'Cancer Drug target gene': 'DrugTarget',
                          'Gene Aliases': 'Aliases'}))

okb_yes = okb[okb['DrugTarget'].astype(str).str.strip().str.lower() == 'yes']

# ---------- 5. Tạo tập hợp Symbol + Alias ----------
def explode_names(row):
    names = {str(row['Symbol']).strip()}
    if pd.notna(row['Aliases']):
        # tách theo ; , hoặc khoảng trắng
        parts = re.split(r'[;,]\s*|\s+', str(row['Aliases']))
        names.update([p.strip() for p in parts if p.strip()])
    return names

symbol_set = set()
okb_yes.apply(lambda r: symbol_set.update(explode_names(r)), axis=1)

# ---------- 6. Gắn cờ & lưu ----------
out = (df_res.merge(map_df, on='Ensembl_ID', how='left')
              .assign(In_OncoKB=lambda d:
                      d['Symbol'].apply(lambda s:
                          'Yes' if pd.notna(s) and s.strip() in symbol_set else 'No'))
              .loc[:, ['Ensembl_ID', 'Symbol', 'Score', 'In_OncoKB']])

out.to_csv(args.output, sep='\t', index=False)
print(out.head(10))

# ---------- 7. Thống kê ----------
hits = (out['In_OncoKB'] == 'Yes').sum()
print(f'✅ Có {hits} / 100 gene khớp OncoKB')
print(f'📄 Kết quả lưu tại: {args.output}')
