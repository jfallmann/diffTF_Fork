# Modify this Snakemake call to your needs.
# If you run the analysis on a cluster, we recommend using a cluster configuration via --cluster-config

echo "Starting Snakemake\n"

# Real run, using 2 cores
snakemake --snakefile ../../src/Snakefile --cores 2 --configfile config.json --directory .
echo "Finished.\n"
