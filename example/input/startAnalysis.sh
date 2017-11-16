# Modify this Snakemake call to your needs.
# By default, a dry-run is performed first. Remove the --dryrun directive to start the analysis
# If you run the analysis on a cluster, we recommend using a cluster configuration via --cluster-config

echo "Starting Snakemake\n"
# Dryrun, using 4 cores
# snakemake --snakefile ../../src/Snakefile --dryrun --cores 2 --configfile config.json

# Real run, using 4 cores
snakemake --snakefile ../../src/Snakefile --cores 2 --configfile config.json --directory .
echo "Finished.\n"
