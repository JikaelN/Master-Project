import os
import pandas as pd
from io import StringIO

# Read VCF
def readVCF(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        
    #Identify header
    header = [line for line in lines if line.startswith("#")]
    
    # Find column name starting with #CHROM
    for line in header:
        
        if line.startswith("#CHROM"):
            column_names = line.strip().split("\t")  # Split columns by tab
            break
    else:
        raise ValueError("No #CHROM line found in the header!")
    
    # Extract the body (all lines not starting with #)
    body = [line for line in lines if not line.startswith("#")]
    
        # Create a DataFrame using the extracted column names
    df = pd.read_csv(
    StringIO("".join(body)),
    sep="\t",
    header=None,
    names=column_names
    )
    
    # Filter out rows with "MULTIALLELIC" in the INFO column
    if "INFO" in df.columns:
        df = df[~df["INFO"].str.contains("MULTIALLELIC", na=False)]
    
    df = df.drop(["#CHROM","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"], axis = 1)

    return df

df = readVCF("R_analysis/Data/neutral_raw/replicate_1/Seed_31995_pop_1.vcf")
print(df)

# Let him cook
#{"replicate" : {"seed" : ["df_po"p]}}

#Spare matrix > faster computation time with only 2 or 1 and nothing for 0 associated with file that keep genome position +
# information on population & replicate