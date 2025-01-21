import argparse
import pickle
import os
import pandas as pd
import numpy as np
import pysam
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import MinMaxScaler
from xgboost import XGBRegressor

# Utility functions
def calculate_evaluation_metrics(y_true, y_pred):
    mse = mean_squared_error(y_true, y_pred)
    rmse = np.sqrt(mse)
    mae = mean_absolute_error(y_true, y_pred)
    r2 = r2_score(y_true, y_pred)
    return mse, rmse, mae, r2

def calculate_roq(rawP, k=10, maxMAPQ=60):
    return maxMAPQ * (1 / (1 + np.exp(-k * (rawP - 0.5))))

# Argument parser
parser = argparse.ArgumentParser(description='Add ROQ tags to BAM file using a trained model.')
parser.add_argument('-i', '--input', help='Feature summary file (tab-delimited)', required=True)
parser.add_argument('-ib', '--inputBam', help='Input original BAM file', required=True)
parser.add_argument('-o', '--outdir', help='Output directory for results', required=True)
parser.add_argument('-m', '--model', help='Trained model to use', required=False)
args = parser.parse_args()

# Paths and directories
input_bam = args.inputBam
feature_file = args.input
outdir = args.outdir

if not os.path.exists(outdir):
    os.makedirs(outdir)

output_bam = os.path.join(outdir, 'output.bam')

# Load or initialize the model
if args.model:
    with open(args.model, 'rb') as file:
        xgb_model = pickle.load(file)
else:
    xgb_model = XGBRegressor(
        colsample_bytree=0.6,
        learning_rate=0.1,
        max_depth=7,
        min_child_weight=5,
        n_estimators=200,
        subsample=1.0
    )

# Load the feature file
df = pd.read_csv(feature_file, sep='\t', header=None)
df.columns = ['READ_NAME', 'MAPQ', 'ALIGN_SCORE', 'SECONDARY_ALIGN_SCORE', 'MISMATCHES','GAP_OPENS', 'GAP_EXT', 'EDIT_DIST', 'MATE_ALIGNMENT_SCORE', 'ALIGN_SCORE_DIFF','INSERT_SIZE', 'READ_GC_CONT', 'N_LOW_QUALITY_BASE', 'AVG_QUALITY_BASE_SCORE']

# Preprocess features
X = df.drop(columns=['READ_NAME'])
imputer = SimpleImputer(strategy='mean')
X_imputed = imputer.fit_transform(X)

scaler = MinMaxScaler()
X_scaled = scaler.fit_transform(X_imputed)

y_pred = xgb_model.predict(X_scaled)


# Process BAM file
samfile = pysam.AlignmentFile(input_bam, 'r')
header = samfile.header.to_dict()

# Add ROQ tag to header
header['PG'].append({
    'ID': 'ROQ',
    'PN': 'ROQ Tag',
    'CL': 'Read Overlapping Quality',
    'VN': '1.0'
})

outfile = pysam.AlignmentFile(output_bam, 'wb', header=header)

# Add ROQ to each read
roq_values = calculate_roq(y_pred)
read_map = {row['READ_NAME']: roq for row, roq in zip(df.to_dict('records'), roq_values)}

for read in samfile:
    if read.query_name in read_map:
        roq_value = read_map[read.query_name]
        read.set_tag('RQ', roq_value, value_type='f')
    outfile.write(read)

samfile.close()
outfile.close()
print(f"Updated Alignment file saved to {output_bam}")
