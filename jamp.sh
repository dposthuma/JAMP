#!/bin/bash
# Jamp 1.0
# Arjen van Bochoven & Danielle Posthuma 1 juni 2012

VERSION=1.0.10
DATE=$(date)
PHENO_HEADER="pheno_header"
PHENO_TYPES="pheno_types"
PHENO_PERM="pheno_perm"
DEC=10 # Number of decimals of logp
PLINKOPTS= # Contains the remaining options for plink
DOSAGE= # Check if we're running with --dosage
ALLPHENO= # Is --all-pheno requested
OUT=plink # Default plink output prefix
PHENO_RES= # File listing of pheno files
JAMPLOG="log" # Jamp logfile
COMB_FILES= # Files to combine

# Usage
usage()
{
cat << EOF

JAMP | v${VERSION}

JAMP is a free, open source tool to run multivariate GWAS.
JAMP depends on PLINK.

Usage: $0 options

JAMP options:
     --jperm   : Number of permutations or permstart-permend (required)
     --out     : JAMP (and PLINK) prefix, recommended
     --jcomb   : Combine empirical P-values from multiple runs
     --help    : Show this message
     --version : Show version information

Options to be passed to PLINK:
      --bfile  : Path to data
      --assoc  : Association model
      --pheno  : Path to pheno file
      --all-pheno  : Specify all phenotypes to be used

Additional PLINK options supported in JAMP:
      --file --tfile --lfile
      --linear --logistic --trend --model --dosage --fam
      --covar --covar-number

EXAMPLES:

./jamp --bfile example --assoc --pheno pheno.txt --all-pheno \\
      --out run1 --jperm 1000

This generates one P-value for each SNP based on the joint information
from all phenotypes in pheno.txt, and based on 1000 permutations.
The created output file is called run1.jamp.empp

./jamp --bfile example --assoc --pheno pheno.txt --all-pheno \\
      --out run2 --jperm 1001-2000

Same as before but instead runs a second set of 1000 permutations,
creating output file run2.jamp.empp

./jamp --jcomb run1.jamp.empp run2.jamp.empp

This combines the information from two runs of 1000 permutations
into one empirical P

(C) 2012 Arjen van Bochoven & Danielle Posthuma, Apache License
For documentation, citation & bug-report instructions:
http://ctglab.nl/software/jamp/
EOF
exit 0;
}

# Pvals gets the P values from column P, calculates the log10
# and writes the resulting val to a file with extension .logp
function pvals {
    log "# Getting P values"
	count=0;
    for file in $PHENO_RES; do
		count=`expr $count + 1`
        awk -v col=P 'NR==1{for(i=1;i<=NF;i++){if($i==col){c=i;break}}} NR>1{val=$c=="NA"?1:$c;printf "%.'$DEC'f\n", log(val)/log(10) * -1}' $file > "${J_PREFIX}.P${count}.logp"
    done 
}

# Timer
function timer()
{
    if [[ $# -eq 0 ]]; then
        echo $(date '+%s')
    else
        local  stime=$1
        etime=$(date '+%s')

        if [[ -z "$stime" ]]; then stime=$etime; fi

        dt=$((etime - stime))
        ds=$((dt % 60))
        dm=$(((dt / 60) % 60))
        dh=$((dt / 3600))
        printf '%d:%02d:%02d' $dh $dm $ds
    fi
}

# Bail out with a message
function bail {
	log "## JAMP ERROR: $1";
	echo "Type $0 --help for instructions";
	exit 1;
}

# Logging
function log {
	echo $1 | tee -a ${J_PREFIX}.${JAMPLOG}
}

# Run plink, if DOSAGE is true, run in a loop for all pheno's
function runplink {
	
	local PHENOFILE=$1
	
	# Reset result files
	PHENO_RES=
	
	# When running --dosage, run plink in a loop with --mpheno n
	# --all-pheno does not work with --dosage in the current 
	# version of plink (v1.07)

	if [[ -z $DOSAGE ]]; then
		log "plink $PLINKOPTS --pheno $PHENOFILE";
		plink $PLINKOPTS --pheno $PHENOFILE > /dev/null || bail 'There is an error running plink'
		# Get result files from log
		PHENO_RES=$(sed -n 's/.*results to.*\[ \(.*\) \].*/\1/p' $OUT.log | tr '\n' ' ')
	else
		for (( i=1; i<=$PHENO_COUNT; i++ )); do
			log "plink ${PLINKOPTS}.P${i} --mpheno $i --pheno $PHENOFILE";
			plink ${PLINKOPTS}.P${i} --mpheno $i --pheno $PHENOFILE > /dev/null || bail 'There is an error running plink'
			# Get result files from log
			PHENO_RES="$PHENO_RES "$(sed -n 's/.*results to.*\[ \(.*\) \].*/\1/p' $OUT.P${i}.log)
		done
	fi
	
	# Check if we have results
	if [[ -z $PHENO_RES ]]; then
		bail "No result files, check your plink options"
	fi
	log "Result files: $PHENO_RES";
}

# JCOMB Combine empp files to a single file
function jcomb {
	echo "--jcomb is not yet implemented";
	exit 1;
}

# Process options, pass all unknown options to plink
while (( "$#" )); do
  case "$1" in
	--version )
		echo "JAMP version $VERSION"; exit 0;;
	--help | -? )
		usage;;
	--jperm )
		PERM=$2; shift; shift;;
	--jcomb )
		while shift; do
			COMB_FILES="$COMB_FILES $1";
		done
		jcomb;;
	--pheno )
		PHENO=$2; shift; shift;;
	--out )
		OUT=$2; shift; shift;;
	--mpheno )
		bail "You specified only one phenotype with --mpheno please use --all-pheno";;
	--sex )
		bail "PLEASE do not have --sex with JAMP, consult the manual";;
	--all-pheno )
		ALLPHENO=1; shift;;
	--dosage )
		DOSAGE=1
		PLINKOPTS="$PLINKOPTS $1"
		while shift; do
			if  [[ $1 != --* ]]; then
				PLINKOPTS="$PLINKOPTS $1";
			else
				break;
			fi
		done;;
	--* )
		PLINKOPTS="$PLINKOPTS $1"
		while shift; do
			if  [[ $1 != --* ]]; then
				PLINKOPTS="$PLINKOPTS $1";
			else
				break;
			fi
		done;;
    * ) # First argument does not have --
		PLINKOPTS="$PLINKOPTS $1"; shift;;
  esac
done

# Check if there are files with current prefix
ls ${OUT}* &> /dev/null && bail "There are already files starting with ${OUT}., please remove them first"

# Jamp prefix
J_PREFIX="${OUT}.jamp"

if [[ -z $PERM ]]; then
	bail "Please specify permutations";
fi

set - $(IFS="-"; echo $PERM)
if [[ -n $2 ]]; then
	PERMSTART=$1;
	PERMEND=$2;
else
	PERMSTART=1;
	PERMEND=$PERM;
fi

# Total number of permutations
TOTALPERM=$((PERMEND-PERMSTART+1));
if [[ $TOTALPERM < 1 ]]; then
	bail "Total permutations should be more than 0"
fi

if [[ ! -e "$PHENO" ]]; then
	bail "Can't find pheno file: $PHENO"
fi

# Get number of phenotypes from pheno file
PHENO_COUNT=$(awk 'NR==1{print NF-2}' "$PHENO");

# Plink does not process --all-pheno when running --dosage
# so we're not including --all-pheno to plinkopts when --dosage
# is present
if [[ -z $DOSAGE ]]; then
	PLINKOPTS="$PLINKOPTS --all-pheno";
fi

# Add out to the end of plinkopts, important!
PLINKOPTS="$PLINKOPTS --out $OUT";

# Check if we can write in OUT
> ${J_PREFIX}.${JAMPLOG} 2>/dev/null || bail "Cannot write to ${J_PREFIX}.${JAMPLOG}.log"

log "## JAMP V${VERSION} Started: $DATE"
log ""
log "## JAMP: Arguments"
log "# plink${PLINKOPTS}"
log "# pheno file: $PHENO"
log "# phenotypes: $PHENO_COUNT"
log "# permutations: $PERMSTART - $PERMEND"

log "## JAMP: Start"
t=$(timer)

# Run plink
log "# Running plink";
runplink $PHENO

# Get P values
pvals

# Add P values
log "# Calculating P value totals"
paste -d'+' ${J_PREFIX}.P*.logp | bc > ${J_PREFIX}.sumlogp

# Split pheno file
log "# Splitting pheno file"
cut -d" " -f1-2 $PHENO > ${J_PREFIX}.$PHENO_HEADER
cut -d" " -f3- $PHENO > ${J_PREFIX}.$PHENO_TYPES

# Split plink output, to be used later
log "# Split plink output"
awk '{print $1, $2}' $(ls $PHENO_RES | awk 'NR==1{print $1}') > ${J_PREFIX}.chr_snp

# Get number of phenotypes and the associated pvals
log '# Get pvals of each phenotype'
count=0;
for file in $PHENO_RES; do
	count=`expr $count + 1`
    awk -v col=P 'NR==1{for(i=1;i<=NF;i++){if($i==col){c=i;break}}} NR>1{print $c}' $file > "${J_PREFIX}.P${count}.pval"
done
paste -d" " ${J_PREFIX}.P*.pval | awk 'NR==1{head="NPHENO";for(i=1;i<=NF;i++){head=head " P_P" i };print head}{print '$PHENO_COUNT', $0}' > ${J_PREFIX}.npheno_pheno_x

# Main loop
log "# Starting main loop"
for (( c=$PERMSTART; c<=$PERMEND; c++ ))
do
	log "## JAMP: Permutation $c of $PERMEND";
	
	# Generate permutated file
	log "# Create perm file";
    ./permutator ${J_PREFIX}.$PHENO_HEADER ${J_PREFIX}.$PHENO_TYPES ${J_PREFIX}.$PHENO_PERM $c
	
	log "# Running plink";
	runplink ${J_PREFIX}.$PHENO_PERM

    # Get P values
	pvals
	
    # Add P values
    log "# Calculating P value totals"
	paste -d'+' ${J_PREFIX}.P*.logp | bc > ${J_PREFIX}.perm$c
done

log '## JAMP: Permutations done'

# Glue perm files together
log '# Glueing perm files'
paste -d' ' ${J_PREFIX}.perm* > ${J_PREFIX}.sumlogp_perm

# Count occurences
log "# Generating SUMLOGP, NPERMS, EMP_P"
awk 'BEGIN{c=1;while((getline n < "'${J_PREFIX}.sumlogp'") > 0){arr[c++] = n}print "SUMLOGP NPERMS EMP_P"}
{c=0;for(i=1;i<=NF;i++){if($i>=arr[FNR]){c++}}print arr[FNR], NF, c/'$TOTALPERM'}' ${J_PREFIX}.sumlogp_perm > ${J_PREFIX}.sum_perm_emp

# CHR SNP NPHENO P_Nn SUMLOGP NPERMS EMP_P
paste -d" " ${J_PREFIX}.chr_snp ${J_PREFIX}.npheno_pheno_x ${J_PREFIX}.sum_perm_emp > ${J_PREFIX}.empp

ELAPSED=$(timer $t)

log "## JAMP: Done. Processing took $ELAPSED. Output in ${J_PREFIX}.empp";
