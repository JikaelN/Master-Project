import re
import pandas as pd

# File paths
files = {
    "T-Test": "./R_analysis/result/norm/wo_maf/within_replication_qqplot/t_test_results.txt",
    "AD Test": "./R_analysis/result/norm/wo_maf/within_replication_qqplot/ad_test_results.txt",
    "ANOVA": "./R_analysis/result/norm/wo_maf/within_replication_qqplot/anova_results.txt",
    "KS Test": "./R_analysis/result/norm/wo_maf/within_replication_qqplot/ks_test_results.txt"
}

# Function to extract p-values from text content
def extract_p_values(file_content, pattern, test_name):
    matches = re.findall(pattern, file_content)
    return [{"Test": test_name, "P-Value": float(match)} for match in matches]

# Patterns to extract p-values
patterns = {
    "T-Test": r"(?<![\w.])\d\.\d+(?=[\s\*]?$)",  # p-values like 0.03 or 0.05
    "AD Test": r"P-value\s*=\s*([\d\.]+)",
    "ANOVA": r"Pr\(>F\)\s*([\d\.]+)",
    "KS Test": r"p-value\s*=\s*([\d\.]+)"
}

# Collect p-values
all_p_values = []

for test_name, file_path in files.items():
    with open(file_path, "r") as file:
        content = file.read()
        pattern = patterns[test_name]
        all_p_values.extend(extract_p_values(content, pattern, test_name))

# Create a DataFrame
p_values_df = pd.DataFrame(all_p_values)

# Summarize p-values
summary = p_values_df.groupby("Test").describe()

# Save and display the summary
summary_path = "wo_maf_within_p_values_summary.csv"
p_values_df.to_csv(summary_path, index=False)