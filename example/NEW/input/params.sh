
# Config file to use
configFile="config.json"

snakefile="/g/scb2/zaugg/carnold/Projects/PWM_Ivan/src/Snakefile"


# IGNORED NOW, AUTOMATICALLY SET TO 5 WHEN LOCAL 16 OTHERWISE
nCores=16 # Maximum number of CPUs per rule executed. See the separate cluster specification

# Enable various additional output of Snakemake such as verbose messages and printing shell commands
useVerbose=false

useConda=true

# Use a dry run for testing purposes?
dryRun=false

forceRerunAll=false

# Use --notemp for developing purposes: Ignore temp() declarations. This is useful when running only a part of the workflow,
# since temp() would lead to deletion of probably needed files by other parts of the workflow.
ignoreTemp=true

# Default false
touchOutputFiles=false

# Default: Use all rules: Set to ""
allowedRules=""

# Start from only a specific rule? If yes, name the rule here, which has to correspond to the rule name in the Snakefile
runSpecificRule=""

# Run also all downstream rules when runSpecificRule is specified? If set to false, ONLY the specified rule will be run; otherwise, all downstream rules are also triggered
runAlsoDownstreamRules=true

# Cluster options
submitToCluster=true
clusterConfig="cluster.json"
# Maximum number of simultaenous jobs
maxJobsCluster=400

# Maximum number of times a job should be reexecuted when failing
maxRestartsPerJob=2

runCustomCommand=""

# Only execute the last part: Calling Snakemake
skipSummaryAndDAG=true

# You can usually leave this untouched unless you also suffer from the "bug" that dot has a problem producing PDF files. Then change to "svg" pr "png"
workflowGraphFileType="pdf"
