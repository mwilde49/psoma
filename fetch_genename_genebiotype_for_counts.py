import os, sys, pandas as pd

if len(sys.argv) != 3:
    print("Usage: python fetch_genename_genebiotype_for_counts.py <CONFIG_DIRECTORY> <REF_FILE>")
    sys.exit(1)

config_directory = sys.argv[1]
ref_file = sys.argv[2]

raw_counts_dir = os.path.join(config_directory, '6_raw_counts_output')
counts_file = os.path.join(config_directory, '6_raw_counts_output', 'raw_htseq_counts.csv')
feature_counts_file = os.path.join(config_directory, '6_raw_counts_output', 'raw_feature_counts.csv')

if not os.path.isdir(raw_counts_dir):
    raise FileNotFoundError(f"Directory does not exist: {raw_counts_dir}")

if not os.path.isfile(counts_file):
    raise FileNotFoundError(f"File does not exist: {counts_file}")

if not os.path.isfile(feature_counts_file):
    raise FileNotFoundError(f"File does not exist: {feature_counts_file}")

if not os.path.isfile(ref_file):
    raise FileNotFoundError(f"File does not exist: {ref_file}")

ref_df = pd.read_csv(ref_file, delimiter=' ', names=['Gene_ID', 'Gene_Name', 'Gene_Biotype'])
counts_df = pd.read_csv(counts_file)
feature_counts_df = pd.read_csv(feature_counts_file)
feature_counts_df.columns = ['Gene_ID'] + list(feature_counts_df.columns[1:])

counts_df = counts_df.merge(ref_df, on=['Gene_ID'], how='left')
feature_counts_df = feature_counts_df.merge(ref_df, on=['Gene_ID'], how='left')

counts_df.to_csv(os.path.join(config_directory, '6_raw_counts_output', 'raw_htseq_counts.csv'), index=False)
feature_counts_df.to_csv(os.path.join(config_directory, '6_raw_counts_output', 'raw_feature_counts.csv'), index=False)
