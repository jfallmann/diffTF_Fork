# Modify this Snakemake call to your needs.
# If you run the analysis on a cluster, we recommend using a cluster configuration via --cluster-config

echo "Starting Snakemake (dry run only, not actually executing anything)\n"
# Dryrun, using 2 cores
snakemake --snakefile ../../src/Snakefile --dryrun --cores 2 --configfile config.json
echo "Finished.\n"
