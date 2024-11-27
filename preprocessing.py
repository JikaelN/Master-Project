import os
import pandas as pd
import random
import numpy as np
import matplotlib.pyplot as plt




# read_vcf and return a DF with all the data
def read_vcf(file_path):
    """
    Reads a VCF file and returns a pandas DataFrame.
    """
    with open(file_path, "r") as f:
        lines = f.readlines()
        # Extract header columns
        column_names = None
        for line in lines:
            if line.startswith("#CHROM"):
                column_names = line.strip().replace("#", "").split("\t")
                break
        # Extract data
        data = [line.strip().split("\t") for line in lines if not line.startswith("#")]
        df = pd.DataFrame(data, columns=column_names)
        #Drop MULTIALLELIC row since they mess up with sampling later
        df = df[~df.INFO.str.contains("MULTIALLELIC")]
    return df

#get common loci and return df for all the population with only common loci
def get_common_loci(list_df):
    #empty loci list that will contain list of loci for all our df
    all_loci_list = list()
    #extract loci
    for df in list_df:
        all_loci_list.append(df.POS)
    
    #get common using set.interaction and map. map convert each list to set which are next operated with intersection
    common = list(set.intersection(*map(set, all_loci_list)))
    
    #get list of df with only common loci
    final_df_list = [df[df.POS.isin(common)] for df in list_df]
    return final_df_list

#Sample 1'000 rows from our data
def th_sample(list_df, log_dir):
    '''
    Sample 1'000 rows from our data
    '''
    #get list pos (loci) from the first df
    loci_list = list(list_df[0].POS)
    # Check if enough common loci
    if len(loci_list) > 1000 :
        # Get 1000 random loci 
        
        sampled = random.sample(loci_list, 1000)
        # Sampled the same loci over all df
        sampled_df = [df[df.POS.isin(sampled)].reset_index(drop=False, inplace=False) for df in list_df]
        
        
        # Log sampled loci
        
        path_to_file = log_dir +"/"+ "log_sample.txt"
        with open(path_to_file, "w") as f:
            f.write(" ".join(sampled))
    
        return sampled_df
    else:
        print("Not enough common loci to sample. Increase the number of individuals in the population.")
        return []
    
    

# Set effect size normally distributed on all our dataset
def norm_effect_size(df_list, output_directory):
    '''
    Convert data in columns 'i0' to the last column into numeric, sum values, and apply a normally distributed effect size,
    while keeping the 'POS' column intact.
    '''
    # normal distribution
    mu = 0 #mean
    sigma = 1 # sd
    size = len(df_list[0]) #sample size
    
    #normal distribution
    dist = np.random.normal(mu, sigma, size)
    #convert it to to Df column-wise to apply it later on the whole df
    
    dist_df = pd.DataFrame({col: dist for col in df_list[0].loc[:, "i0":].columns})
    
    #save plot distribution in output dir
    output_dir = output_directory + "/"+ "normal_distribution.png"
    # Plot the distribution
    plt.figure(figsize=(8, 6))
    plt.hist(dist, bins=30, density=True, alpha=0.7, color='blue', edgecolor='black')
    plt.title("Normal Distribution (mu=0, sigma=1)")
    plt.xlabel("Value")
    plt.ylabel("Density")
    plt.grid(axis='y', alpha=0.75)
    
    plt.savefig(output_dir)
    plt.close
    
    # transform i0 -> in info into numeric, sum them, and apply the distribution
    for idx , df in enumerate(df_list):
        pos_column = df[['POS']]

        i_data = df.loc[:, 'i0':] #get all data for all individuals
        i_data_sum = i_data.map(lambda x: sum(map(int, x.split("|")))) #split by | transform number into int and sum them
        
         # Apply the effect size to the summed data
        i_data_normalized = i_data_sum * dist_df
        
        # Merge back the 'POS' column
        result_df = pd.concat([pos_column, i_data_normalized], axis=1)
        
        # Save the transformed DataFrame to the output directory
        output_path = f"{output_directory}/pop_{idx}.csv"
        result_df.to_csv(output_path, index=False)
        print(f"Transformed DataFrame {idx} saved to: {output_path}")
        
# Process VCF files for multiple samples
def process_vcf_files(directory, sample_base_dir):
    """
    Processes all VCF files in a replicate directory for 10 samples.
    """
    files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(".vcf")]
    df_list = [read_vcf(file) for file in files]
    common_df_list = get_common_loci(df_list)

    if not common_df_list:
        print(f"No valid data in {directory}. Skipping transformations.")
        return

    # Perform 10 transformations
    for i in range(1, 11):
        sample_dir = os.path.join(sample_base_dir, f"Sample_{i}", os.path.basename(directory))
        os.makedirs(sample_dir, exist_ok=True)
        print(f"Processing sample {i} for {directory}...")

        sampled_data = th_sample(common_df_list, sample_dir)
        if sampled_data:
            norm_effect_size(sampled_data, sample_dir)

    # Delete original VCF files after processing
    for file in files:
        os.remove(file)
        print(f"Deleted VCF file: {file}")

        
if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python preprocessing.py <replicate_directory> <sample_base_directory>")
        sys.exit(1)

    replicate_directory = sys.argv[1]
    sample_base_directory = sys.argv[2]

    if not os.path.isdir(replicate_directory):
        print(f"Error: {replicate_directory} is not a valid directory.")
        sys.exit(1)

    print(f"Processing directory: {replicate_directory}")
    process_vcf_files(replicate_directory, sample_base_directory)






