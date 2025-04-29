import pandas as pd
import pysam
import numpy as np


# Function definition extract doage from vcf file

def extract_ds_from_vcf(vcf_file_path,file_path):
  # Open the VCF file using pysam
  vcf_file = pysam.VariantFile(vcf_file_path)

# Extract header information
  header = vcf_file.header
  sample_names = header.samples

# Initialize a list to store dictionaries with variant information and DS values for each sample
  variant_data = []

# Iterate over variants in the VCF file
  for record in vcf_file:
  # Extract variant information
      variant_info = {
        'CHROM': record.chrom,
        'POS': record.pos,
        'ID': record.id,
        'REF': record.ref,
        'ALT': record.alts[0] if record.alts else '.'  
      }

# Extract DS values for each sample
    #sample_ds_values = {sample: record.samples[sample]['DS'] for sample in sample_names}
      sample_ds_values = {sample: float(ds) for ds in record.samples[sample]['DS'] for sample in sample_names}
# Combine variant information and DS values
      variant_info.update(sample_ds_values)
      variant_data.append(variant_info)

# Create a DataFrame with variant information and DS values
  ds_df = pd.DataFrame(variant_data)

  ds_df['chrom_pos'] = ds_df['ID'] + ':' + ds_df['REF'] + ':' + ds_df['ALT']

# PRS Calculation
  snp_weight = pd.read_table(file_path)

  result = ds_df[ds_df['chrom_pos'].isin(snp_weight['newid'])]
  weight = snp_weight[snp_weight['newid'].isin(result['chrom_pos'])]['beta_grid4']

  sum=np.dot(result['SAMPLE'].astype(float),weight)

# Return PRS per Chromosomes
  return sum



#######################
# Main Script


file_path = '/fh/fast/peters_u/GECCO_Working/mintaworking/toelisabeth/trans_prs_Nov_19.txt'
vcf_file_path = '/fh/fast/peters_u/GECCO_Working/robert_working/nextflow/FH-PRS-tool-example/prs_prototype_test/output/'

chromosomes = [f'{i}' for i in range(1, 23)]

prs=0
for chrom in chromosomes:
  ds_data_chrom = extract_ds_from_vcf(vcf_file_path + str(chrom) +'.imputed.vcf.gz', file_path)
  prs=prs+ds_data_chrom
  
# Display PRS for the sample
print(prs)