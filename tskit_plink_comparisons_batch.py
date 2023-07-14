import numpy as np
import tskit
import tszip
import time
import os
import sys
import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("-n", "--num_samples", type=int)

args = parser.parse_args()
sample_size = args.num_samples
print(f"input args: sample size={sample_size}")
#np.set_printoptions(threshold=sys.maxsize)
win_sizes = [10, 100]
#load a simulated tree sequence from the Quebec genealogies paper
ts = tszip.decompress(
    "/nfs_home/users/osvk/projects/tskit-ld/data/simulated_chrom_22.ts.tsz"
)

def subsample_trees(N, ts, seed=123):
    np.random.seed(seed)
    indivs = np.random.choice(ts.individuals(), N, replace=False)
    ind_nodes = [] 
    for ind in indivs:
        ind_nodes.extend(ind.nodes)
    return ts.simplify(samples=ind_nodes)

def ts_to_vcf(ts_to_conv, vcf_file, log_file):
    start = time.perf_counter()
    with gzip.open(vcf_file, "wt") as vcf_file:
        ts_to_conv.write_vcf(vcf_file)
    end = time.perf_counter()

    with open(log_file, "a") as log_file:
        log_file.write(f"Time to write VCF for {ts.num_samples} samples: {end - start:.2f} seconds\n")

def vcf_to_plink(vcf_file, plink_file, log_file):
    start = time.perf_counter()
    os.system(f"plink --vcf {vcf_file} --out data/bed/{plink_file} --double-id --vcf-half-call 'h'")
    end = time.perf_counter()
    with open(log_file, "a") as log_file:
        log_file.write(f"Time to convert VCF to PLINK for {vcf_file}: {end - start:.2f} seconds\n")

#we make a reasonable assumption here that sample_size has to be even
#so we get exactly n/2 diploid genomes
num_of_sampled_indivs = int(sample_size/2)
ts_subset = subsample_trees(num_of_sampled_indivs, ts)
ts_gene = ts_subset.keep_intervals([[49_000_000, 50_000_000]])
non_biallelic_sites = np.where(np.bincount(ts_gene.mutations_site, minlength=ts_gene.num_sites) != 1)[0]
ts_gene_biallelic = ts_gene.delete_sites(non_biallelic_sites)

ts_to_vcf(ts_gene_biallelic, f"data/vcf/ukb_chr22_N{sample_size}.vcf.gz", 
          log_file=f"logs/conversions/file_conversions_N{sample_size}.log")
vcf_to_plink(f"data/vcf/ukb_chr22_N{sample_size}.vcf.gz", 
             f"ukb_chr22_N{sample_size}",
             f"logs/conversions/file_conversions_N{sample_size}.log")

plink_file = f"data/bed/ukb_chr22_N{sample_size}"

for win_size in win_sizes:
    tskit_start = time.perf_counter()
    ldcalculator = tskit.LdCalculator(ts_gene_biallelic)
    r2_tskit = ldcalculator.r2_matrix()
    #with open(f"data/tskit_ld/ukb_chr22_N{sample_size}", "w") as out_file:
    #    print(r2_tskit, file=out_file)
    tskit_time_taken = time.perf_counter() - tskit_start

    plink_start = time.perf_counter()
    os.system(f"plink --bfile {plink_file} --r2 --ld-window-r2 0 --ld-window {win_size} --ld-window-kb 1e3 --threads 1 --out data/plink_ld/ukb_chr22_N{sample_size}")
    plink_time_taken = time.perf_counter() - plink_start
    with open(f"logs/compute_times/calc_times_N{sample_size}_w{win_size}.log", "a") as log_file:
        log_file.write(f"sample_sizes:{sample_size}\n")
        log_file.write(f"tskit:{tskit_time_taken}\n")
        log_file.write(f"plink:{plink_time_taken}\n")
