#! /bin/bash
export CLASSPATH=$4
echo ${CLASSPATH}

echo $6

thresh="_threshold"
thresh=$5$thresh
perm="_permnet"
perm=$5$perm
rand="_randscore"
rand=$5$rand

# currently issue where distance threshold is not set
java be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser --method ensemble --input $1 --guessingparam 100.0 --stand col_norm --matrixtype count --format tab_table --resamplemethod shuffle_rows --edgethreshold 0.05 --verbosity FATAL --iterations 100 --inference mrnet --renorm --multigraph --thresholdguessing edgeNumber --topbottom --correlnonrandp none --measure1 conf --scoremergestrategy mean --minsupport 2 --measure2 supp --multicorr none --min 1 --max 3 --ensemblemethods correl_pearson/correl_spearman/sim_mutInfo/dist_bray/dist_kullbackleibler --minetdisc equalfreq --kernelwidth 0.25 --metadataattribs Kingdom/Phylum/Class/Order/Family/Genus/Species --minetmiestimator mi.shrink --metadata $2 --networkmergestrategy union --nantreatmentparam 1 --nantreatment none --filter rand  --output $thresh
java be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser --method ensemble --input $1 --stand col_norm --matrixtype count --format tab_table --resamplemethod shuffle_rows --edgethreshold 0.05 --verbosity FATAL --iterations 100 --inference mrnet --renorm --multigraph --randroutine edgeScores --correlnonrandp none --measure1 conf --scoremergestrategy mean --minsupport 2 --measure2 supp --multicorr none --min 1 --max 3 --ensemblemethods correl_pearson/correl_spearman/sim_mutInfo/dist_bray/dist_kullbackleibler --minetdisc equalfreq --kernelwidth 0.25 --metadataattribs Kingdom/Phylum/Class/Order/Family/Genus/Species --minetmiestimator mi.shrink --metadata $2 --networkmergestrategy union --nantreatmentparam 1 --nantreatment none --filter rand --output $perm --ensembleparamfile $thresh

java be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser --method ensemble --edgethreshold 0.05 --verbosity FATAL --iterations 100 --inference mrnet --pvaluemerge brown --resamplemethod bootstrap --stand col_norm --format tab_table --matrixtype count --guessingparam $6 --input $1 --multigraph --topbottom --ensemblemethods correl_pearson/correl_spearman/sim_mutInfo/dist_bray/dist_kullbackleibler --max 3 --thresholdguessing edgeNumber --nantreatment none --nantreatmentparam 1 --networkmergestrategy union --metadataattribs Kingdom/Phylum/Class/Order/Family/Genus/Species --minetmiestimator mi.shrink --minetdisc equalfreq --kernelwidth 0.25 --metadata $2 --minsupport 2 --multicorr benjaminihochberg --correlnonrandp none --scoremergestrategy mean --measure1 conf --measure2 supp --min 1 --filter rand/confidence_boot --filterparameter 10.0  --output $thresh
java be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser --method ensemble --edgethreshold 0.05 --verbosity FATAL --iterations 100 --inference mrnet --pvaluemerge brown --resamplemethod bootstrap --stand col_norm --format tab_table --matrixtype count --input $1 --multigraph --randroutine edgeScores --ensemblemethods correl_pearson/correl_spearman/sim_mutInfo/dist_bray/dist_kullbackleibler --max 3 --nantreatment none --nantreatmentparam 1 --networkmergestrategy union --metadataattribs Kingdom/Phylum/Class/Order/Family/Genus/Species --minetmiestimator mi.shrink --minetdisc equalfreq --kernelwidth 0.25 --metadata $2 --minsupport 2 --multicorr benjaminihochberg --correlnonrandp none --scoremergestrategy mean --measure1 conf --measure2 supp --min 1 --filter rand/confidence_boot --filterparameter 10.0 --output $3 --ensembleparamfile $thresh

echo "Finished a network!";
sleep 10