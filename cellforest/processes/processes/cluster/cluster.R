args <- commandArgs(trailingOnly = TRUE)

input_metadata_path <- args[1]
input_rds_path <- args[2]
output_rds_path <- args[3]
output_clusters_path <- args[4]
num_pcs <- as.numeric(args[5])
resolution <- as.numeric(args[6])
eps <- as.numeric(args[7])

r_functions_filepath <- args[8]
print("eps"); print(eps)
print("functions"); print(r_functions_filepath)
source(r_functions_filepath)

print("metadata filter")
srat <- metadata_filter_paths(input_metadata_path, input_rds_path)
print("clustering"); print(date())
srat <- find_clusters(srat, output_clusters_path, num_pcs, resolution, eps)
print("clustering DONE"); print(date())
