#!/bin/bash

set -e

mkdir -p output

# Preprocess
ml-exec-process ephys.whiten \
	--inputs timeseries:dataset/raw_data.mda \
	--outputs timeseries_out:output/pre.mda
ml-exec-process ephys.mask_out_artifacts \
	--inputs timeseries:output/pre.mda \
	--outputs timeseries_out:output/pre2.mda \
	--parameters \
		threshold:10 \
		chunk_size:2000 \

# Spike sorting
ml-exec-process ms4alg.sort \
	--inputs \
		timeseries:output/pre2.mda geom:dataset/geom.csv \
	--outputs \
		firings_out:output/firings.mda \
	--parameters \
		detect_sign:0 \
		adjacency_radius:0 \
		detect_interval:30 \
		detect_threshold:4 \
		max_num_clips_for_pca:2000 \ 

# Compute templates
ml-exec-process ephys.compute_templates \
	--inputs timeseries:dataset/raw_data.mda firings:output/firings.mda \
	--outputs templates_out:output/templates.mda \
	--parameters \
		clip_size:150

ml-exec-process ephys.compute_cluster_metrics \
    --inputs firings:output/firings.mda timeseries:dataset/raw_data.mda \
    --outputs metrics_out:output/metrics.json \
    --parameters \
        clip_size:150 \
        samplerate:30000.0 \
        refrac_msec:1.0 \


