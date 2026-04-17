# Strain-level analysis module

root := justfile_directory()

db := "/ei/projects/8/88e80936-2a5d-4f4a-afab-6f74b374c765/data/DB/MarkerG/protal/protal_r226_archaea"
map := root / "data/input/reads/strain/protal_map.tsv"
outdir := root / "data/intermediate/strainprofiling/protal051"

# Run protal 0.5.1 strain profiling (profile part)
protal051-profile start="1" end="200":
	mkdir -p {{outdir}}
	{{root}}/code/scripts/strainprofiling/protal051.sh \
		{{map}} \
		{{db}} \
		{{outdir}} \
		{{start}} \
		{{end}}

# Run protal 0.5.1 strain profiling *msa part)
protal051-msa:
	mkdir -p {{outdir}}
	{{root}}/code/scripts/strainprofiling/protal051_strain.sh \
		{{map}} \
		{{db}} \
		{{outdir}}


# Submit protal 0.5.1 strain profiling to SLURM (single batch)
submit-protal051-profile start="0" end="":
	sbatch \
		-o "{{root}}/logs/strainprofiling/protal051.out" \
		-e "{{root}}/logs/strainprofiling/protal051.err" \
		--mem=128G \
		-J strain-protal051-profile \
		-c 32 \
		-p ei-long,qib-long \
		--wrap="just strain::protal051-profile {{start}} {{end}}"

# Submit protal 0.5.1 strain profiling to SLURM (MSA)
submit-protal051-msa:
	sbatch \
		-o "{{root}}/logs/strainprofiling/protal051-msa.out" \
		-e "{{root}}/logs/strainprofiling/protal051-msa.err" \
		--mem=128G \
		-J strain-protal051-msa \
		-c 32 \
		-p ei-long,qib-long \
		--wrap="just strain::protal051-msa"

processed_dir  := root / "data/intermediate/strainprofiling/processed/protal051"
per_gene_dir   := processed_dir / "per_gene"

# Compute per-gene distance matrices for all MSA files (protal051)
per-gene-distances min_overlap="30":
	mkdir -p {{per_gene_dir}}
	for msa in {{root}}/data/intermediate/strainprofiling/protal051/strains/s__*.msa.fna.gz; do \
		species=$(basename "$msa" .msa.fna.gz); \
		partition="{{root}}/data/intermediate/strainprofiling/protal051/strains/${species}.partition.txt"; \
		[ -f "$partition" ] || continue; \
		python3 {{root}}/code/scripts/strainprofiling/msa_per_gene_distances.py \
			"$msa" "$partition" \
			--min-overlap {{min_overlap}} \
			-o {{per_gene_dir}}; \
	done

# Submit per-gene distance computation to SLURM (one job per species)
submit-per-gene-distances-parallel min_overlap="30":
	mkdir -p {{root}}/logs/strainprofiling {{per_gene_dir}}
	for msa in {{root}}/data/intermediate/strainprofiling/protal051/strains/s__*.msa.fna.gz; do \
		species=$(basename "$msa" .msa.fna.gz); \
		partition="{{root}}/data/intermediate/strainprofiling/protal051/strains/${species}.partition.txt"; \
		[ -f "$partition" ] || continue; \
		sbatch \
			-o "{{root}}/logs/strainprofiling/per_gene_${species}.out" \
			-e "{{root}}/logs/strainprofiling/per_gene_${species}.err" \
			--mem=16G -J "pgd-${species}" -c 1 \
			-p ei-short,qib-short \
			--wrap="python3 {{root}}/code/scripts/strainprofiling/msa_per_gene_distances.py '$msa' '$partition' --min-overlap {{min_overlap}} -o {{per_gene_dir}}"; \
	done

# Compute pairwise distance matrices for all MSA files (protal051)
pairwise-distances min_overlap="200":
	mkdir -p {{processed_dir}}
	for msa in {{root}}/data/intermediate/strainprofiling/protal051/strains/s__*.msa.fna.gz; do \
		species=$(basename "$msa" .msa.fna.gz); \
		python3 {{root}}/code/scripts/strainprofiling/msa_pairwise_distance.py \
			"$msa" \
			--min-overlap {{min_overlap}} \
			-o {{processed_dir}}/${species}.dist.tsv; \
	done

# Submit pairwise distance computation to SLURM
submit-pairwise-distances min_overlap="200":
	mkdir -p {{root}}/logs/strainprofiling
	sbatch \
		-o "{{root}}/logs/strainprofiling/pairwise_distances.out" \
		-e "{{root}}/logs/strainprofiling/pairwise_distances.err" \
		--mem=32G \
		-J strain-pairwise-distances \
		-c 1 \
		-p ei-medium,qib-medium \
		--wrap="just strain::pairwise-distances {{min_overlap}}"

# Submit pairwise distance computation as parallel SLURM array (one job per species)
submit-pairwise-distances-parallel min_overlap="200":
	mkdir -p {{root}}/logs/strainprofiling {{processed_dir}}
	for msa in {{root}}/data/intermediate/strainprofiling/protal051/strains/s__*.msa.fna.gz; do \
		species=$(basename "$msa" .msa.fna.gz); \
		sbatch \
			-o "{{root}}/logs/strainprofiling/pairwise_dist_${species}.out" \
			-e "{{root}}/logs/strainprofiling/pairwise_dist_${species}.err" \
			--mem=8G \
			-J "pdist-${species}" \
			-c 1 \
			-p ei-short,qib-short \
			--wrap="python3 {{root}}/code/scripts/strainprofiling/msa_pairwise_distance.py '$msa' --min-overlap {{min_overlap}} -o {{processed_dir}}/${species}.dist.tsv"; \
	done

iqtree_dir := root / "data/intermediate/strainprofiling/iqtree/protal051"

# Run IQ-TREE3 (--fast, partitioned) on MSAs that have a partition file
# Usage: just strain::iqtree-partitioned [species] [threads]
iqtree-partitioned species="" threads="AUTO":
	mkdir -p {{iqtree_dir}}/partitioned
	pixi run -e iqtree311 \
		{{root}}/code/scripts/strainprofiling/run_iqtree_partitioned.sh \
		{{root}}/data/intermediate/strainprofiling/protal051/strains \
		{{iqtree_dir}}/partitioned \
		{{threads}} "{{species}}"

# Run IQ-TREE3 (--fast, no partition) on raw MSAs
# Usage: just strain::iqtree-unpartitioned [species] [threads]
iqtree-unpartitioned species="" threads="AUTO":
	mkdir -p {{iqtree_dir}}/unpartitioned
	pixi run -e iqtree311 \
		{{root}}/code/scripts/strainprofiling/run_iqtree_unpartitioned.sh \
		{{root}}/data/intermediate/strainprofiling/protal051/strains \
		{{iqtree_dir}}/unpartitioned \
		{{threads}} "{{species}}" raw

# Run IQ-TREE3 (--fast, no partition) on processed MSAs
# Usage: just strain::iqtree-unpartitioned-processed [species] [threads]
iqtree-unpartitioned-processed species="" threads="AUTO":
	mkdir -p {{iqtree_dir}}/unpartitioned/processed
	pixi run -e iqtree311 \
		{{root}}/code/scripts/strainprofiling/run_iqtree_unpartitioned.sh \
		{{root}}/data/intermediate/strainprofiling/protal051/strains \
		{{iqtree_dir}}/unpartitioned \
		{{threads}} "{{species}}" processed

# Submit IQ-TREE3 partitioned to SLURM (one job per species)
submit-iqtree-partitioned-parallel:
	mkdir -p {{root}}/logs/strainprofiling {{iqtree_dir}}/partitioned
	for msa in {{root}}/data/intermediate/strainprofiling/protal051/strains/s__*.msa.fna.gz; do \
		species=$(basename "$msa" .msa.fna.gz); \
		partition="{{root}}/data/intermediate/strainprofiling/protal051/strains/${species}.partition.txt"; \
		[ -f "$partition" ] || continue; \
		[ -f "{{iqtree_dir}}/partitioned/${species}.treefile" ] && continue; \
		sbatch \
			-o "{{root}}/logs/strainprofiling/iqtree_part_${species}.out" \
			-e "{{root}}/logs/strainprofiling/iqtree_part_${species}.err" \
			--mem=32G -J "iqtree-part-${species}" -c 4 \
			-p ei-long,qib-long \
			--wrap="pixi run -e iqtree311 \
				{{root}}/code/scripts/strainprofiling/run_iqtree_partitioned.sh \
				{{root}}/data/intermediate/strainprofiling/protal051/strains \
				{{iqtree_dir}}/partitioned 4 ${species}"; \
	done

# Submit IQ-TREE3 unpartitioned (raw) to SLURM (one job per species)
submit-iqtree-unpartitioned-parallel:
	mkdir -p {{root}}/logs/strainprofiling {{iqtree_dir}}/unpartitioned
	for msa in {{root}}/data/intermediate/strainprofiling/protal051/strains/s__*.msa.fna.gz; do \
		[[ "$msa" == *.processed.msa.fna.gz ]] && continue; \
		species=$(basename "$msa" .msa.fna.gz); \
		[ -f "{{iqtree_dir}}/unpartitioned/${species}.treefile" ] && continue; \
		sbatch \
			-o "{{root}}/logs/strainprofiling/iqtree_nopart_${species}.out" \
			-e "{{root}}/logs/strainprofiling/iqtree_nopart_${species}.err" \
			--mem=16G -J "iqtree-nopart-${species}" -c 4 \
			-p ei-long,qib-long \
			--wrap="pixi run -e iqtree311 \
				{{root}}/code/scripts/strainprofiling/run_iqtree_unpartitioned.sh \
				{{root}}/data/intermediate/strainprofiling/protal051/strains \
				{{iqtree_dir}}/unpartitioned 4 ${species} raw"; \
	done

# Submit IQ-TREE3 unpartitioned (processed) to SLURM (one job per species)
submit-iqtree-unpartitioned-processed-parallel:
	mkdir -p {{root}}/logs/strainprofiling {{iqtree_dir}}/unpartitioned/processed
	for msa in {{root}}/data/intermediate/strainprofiling/protal051/strains/s__*.processed.msa.fna.gz; do \
		species=$(basename "$msa" .processed.msa.fna.gz); \
		[ -f "{{iqtree_dir}}/unpartitioned/processed/${species}.treefile" ] && continue; \
		sbatch \
			-o "{{root}}/logs/strainprofiling/iqtree_nopart_proc_${species}.out" \
			-e "{{root}}/logs/strainprofiling/iqtree_nopart_proc_${species}.err" \
			--mem=16G -J "iqtree-noprt-proc-${species}" -c 4 \
			-p ei-long,qib-long \
			--wrap="pixi run -e iqtree311 \
				{{root}}/code/scripts/strainprofiling/run_iqtree_unpartitioned.sh \
				{{root}}/data/intermediate/strainprofiling/protal051/strains \
				{{iqtree_dir}}/unpartitioned 4 ${species} processed"; \
	done

# Submit protal 0.5.1 strain profiling in (PORIFLES) batches of 10 (1-10, 11-20, ..., 191-200)
submit-protal051-profile-batches:
	mkdir -p {{root}}/logs/strainprofiling
	for i in $(seq 0 19); do \
		start=$(( i * 10 + 1 )); \
		end=$(( i * 10 + 10 )); \
		sbatch \
			-o "{{root}}/logs/strainprofiling/protal051_profile_${start}-${end}.out" \
			-e "{{root}}/logs/strainprofiling/protal051_profile_${start}-${end}.err" \
			--mem=128G \
			-J "strain-protal051-profile_${start}-${end}" \
			-c 32 \
			-p ei-long,qib-long \
			--wrap="just strain::protal051-profile ${start} ${end}"; \
	done
