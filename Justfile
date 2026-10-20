set shell := ["bash", "-cu"]

# Benchpro binary for the examples; set BENCHPRO to use another one (e.g. an installed copy).
benchpro := env("BENCHPRO", "target/release/benchpro")

# Build target/release/benchpro for the examples (skipped when BENCHPRO is set).
[group('examples')]
build-for-examples:
    [ -n "${BENCHPRO:-}" ] || cargo build --release -q

# Run all examples.
[group('examples')]
examples: example-minimal example-full example-gtdb-and-ncbi example-mistakes example-matrix example-other-commands

# Smallest meta file (5 columns); taxonomy inferred from the files.
[group('examples')]
example-minimal: build-for-examples
    mkdir -p examples/01-minimal/output
    # 1. Check the meta file and see the inferred taxonomy of each row.
    {{benchpro}} profile --meta examples/01-minimal/meta.tsv --validate-meta
    # 2. Run the benchmark.
    {{benchpro}} profile --meta examples/01-minimal/meta.tsv --outprefix examples/01-minimal/output/run

# Every optional column, --adjusted and --allow-alternatives.
[group('examples')]
example-full: build-for-examples
    mkdir -p examples/02-full/output
    # --adjusted: phylogeny-adjusted TP/FP/FN (needs GoldStdTree)
    # --allow-alternatives: alternative-name matching
    {{benchpro}} profile --meta examples/02-full/meta.tsv --outprefix examples/02-full/output/run --adjusted --allow-alternatives

# GTDB and NCBI tools in one run; Tool = tool + database.
[group('examples')]
example-gtdb-and-ncbi: build-for-examples
    mkdir -p examples/03-gtdb-and-ncbi/output
    # No Taxonomy column: each row's taxonomy is inferred and checked per tool.
    {{benchpro}} profile --meta examples/03-gtdb-and-ncbi/meta.tsv --outprefix examples/03-gtdb-and-ncbi/output/run

# Six common mistakes and the message each one produces.
[group('examples')]
example-mistakes: build-for-examples
    mkdir -p examples/04-mistakes/output
    # Each meta file contains one mistake; every command is expected to fail.
    @echo "=== a: required column GoldStd is missing"
    {{benchpro}} --log-level error profile --meta examples/04-mistakes/a-missing-column.tsv --outprefix examples/04-mistakes/output/a; echo "--> exit code $?"
    @echo "=== b: GTDB prediction paired with an NCBI gold standard"
    {{benchpro}} --log-level error profile --meta examples/04-mistakes/b-wrong-gold-standard.tsv --outprefix examples/04-mistakes/output/b; echo "--> exit code $?"
    @echo "=== c: one tool, two taxonomies (an NCBI file in a GTDB tool)"
    {{benchpro}} --log-level error profile --meta examples/04-mistakes/c-tool-mixes-taxonomies.tsv --outprefix examples/04-mistakes/output/c; echo "--> exit code $?"
    @echo "=== d: empty required cell (Tool)"
    {{benchpro}} --log-level error profile --meta examples/04-mistakes/d-empty-cell.tsv --outprefix examples/04-mistakes/output/d; echo "--> exit code $?"
    @echo "=== e: --adjusted without any GoldStdTree"
    {{benchpro}} --log-level error profile --meta examples/04-mistakes/e-adjusted-without-tree.tsv --outprefix examples/04-mistakes/output/e --adjusted; echo "--> exit code $?"
    @echo "=== f: output folder does not exist"
    {{benchpro}} --log-level error profile --meta examples/04-mistakes/e-adjusted-without-tree.tsv --outprefix examples/04-mistakes/no-such-folder/run; echo "--> exit code $?"

# Abundance matrices (samples as columns) paired by regex.
[group('examples')]
example-matrix: build-for-examples
    mkdir -p examples/05-matrix/output
    # ProfileRegex/GoldStdRegex capture groups pair column GAST3 with gastrointestinal_3.
    {{benchpro}} profile --meta examples/05-matrix/meta.tsv --outprefix examples/05-matrix/output/run

# The normalize, merge and strain subcommands.
[group('examples')]
example-other-commands: build-for-examples
    mkdir -p examples/06-other-commands/output
    # normalize: detect the tool that wrote a profile ...
    {{benchpro}} --log-level warning normalize -i data/profile_examples/sylph090.profile --detect
    {{benchpro}} --log-level warning normalize -i data/profile_examples/metaphlan402.profile --detect
    {{benchpro}} --log-level warning normalize -i data/profile_examples/motus404.relab --detect
    # ... and convert it to CAMI format.
    {{benchpro}} normalize -i data/profile_examples/sylph090.profile -o examples/06-other-commands/output/sylph.cami.profile --output-format cami
    # merge: per-sample profiles into one species x sample matrix; regex groups name the samples.
    {{benchpro}} merge --input data/test_data/profiles/predictions/sylph/Gastrointestinal_{0,1,2}.profile --output examples/06-other-commands/output/sylph_species_matrix.tsv --rank species --sample-regex '.*/(.*?)/Gastrointestinal_(\d+)\.profile'
    # strain: monophyly and MSA error of strain-level trees against a gold standard.
    {{benchpro}} strain --meta examples/06-other-commands/strain-samplesheet.tsv --outprefix examples/06-other-commands/output/strain

install:
    cargo install --locked --path . --root ~/.local
    install -m 755 scripts/benchpro_visualize_monophyly.R ~/.local/bin/benchpro_visualize_monophyly

sylph-debug:
    mkdir -p sandbox/output
    cargo run -- profile --meta ./data/input/meta_debug_sylph_mouse.xlsx --outprefix sandbox/output/sylph_debug
    
merge-test-r214:
    #!/usr/bin/env bash
    set -euo pipefail
    shopt -s nullglob
    cargo run -- merge \
        --input data/test_data/profiles/predictions/*/GTDB/r214/{*,*/*,*/*/*}.{tsv,bracken,profile} \
        --output data/test_data/profiles/predictions/matrix/r214.tsv \
        --sample-regex ".*/(.*?)/GTDB/r214/human/(.*?)/(\d+)_.*" \
        --sample-regex ".*/(.*?)/GTDB/r214/(.*?)/.*?_(\d+)_.*" \
        --sample-regex ".*/(.*?)/GTDB/r214/(.*?)/(\d+).profile" \
        --sample-regex ".*/(.*?)/GTDB/r214/human/(.*?)_(\d+).profile" \
        --sample-regex ".*/(.*?)/GTDB/r214/(.*?)_?(\d+).profile" \
        --test-sample-regex

merge-test-r207:
    #!/usr/bin/env bash
    set -euo pipefail
    shopt -s nullglob
    cargo run -- merge \
        --input data/test_data/profiles/predictions/*/GTDB/r207/{*,*/*,*/*/*}.{tsv,bracken,profile} \
        --output data/test_data/profiles/predictions/matrix/r207.tsv \
        --sample-regex ".*/(.*?)/GTDB/r207/human/(.*?)/(\d+)\..*" \
        --sample-regex ".*/(.*?)/GTDB/r207/(.*?)/.*?_(\d+)_.*" \
        --sample-regex ".*/(.*?)/GTDB/r207/(.*?)_?(\d+).profile" \
        --sample-regex ".*/(.*?)/GTDB/r207/(.*?)/.*?(\d+)[_.].*" \
        --test-sample-regex

merge-r207:
    #!/usr/bin/env bash
    set -euo pipefail
    shopt -s nullglob
    cargo run -- merge \
        --input data/test_data/profiles/predictions/*/GTDB/r207/{*,*/*,*/*/*}.{tsv,bracken,profile} \
        --output data/test_data/profiles/predictions/matrix/r207.tsv \
        --sample-regex ".*/(.*?)/GTDB/r207/human/(.*?)/(\d+)\..*" \
        --sample-regex ".*/(.*?)/GTDB/r207/(.*?)/.*?_(\d+)_.*" \
        --sample-regex ".*/(.*?)/GTDB/r207/(.*?)_?(\d+).profile" \
        --sample-regex ".*/(.*?)/GTDB/r207/(.*?)/.*?(\d+)[_.].*"

merge-r214:
    #!/usr/bin/env bash
    set -euo pipefail
    shopt -s nullglob
    cargo run -- merge \
        --input data/test_data/profiles/predictions/*/GTDB/r214/{*,*/*,*/*/*}.{tsv,bracken,profile} \
        --output data/test_data/profiles/predictions/matrix/r214.tsv \
        --sample-regex ".*/(.*?)/GTDB/r214/human/(.*?)/(\d+)_.*" \
        --sample-regex ".*/(.*?)/GTDB/r214/(.*?)/.*?_(\d+)_.*" \
        --sample-regex ".*/(.*?)/GTDB/r214/(.*?)/(\d+).profile" \
        --sample-regex ".*/(.*?)/GTDB/r214/human/(.*?)_(\d+).profile" \
        --sample-regex ".*/(.*?)/GTDB/r214/(.*?)_?(\d+).profile"


merge-test-ncbi:
    #!/usr/bin/env bash
    set -euo pipefail
    shopt -s nullglob
    cargo run -- merge \
        --input data/test_data/profiles/predictions/*/NCBI/{*,*/*,*/*/*}.{tsv,bracken,profile} \
        --output data/test_data/profiles/predictions/matrix/ncbi.tsv \
        --sample-regex ".*/(.*?)/NCBI/human/(.*?)/(\d+)\..*" \
        --sample-regex ".*/(.*?)/NCBI/(.*?)/.*?_(\d+)_.*" \
        --sample-regex ".*/(.*?)/NCBI/(.*?)_?(\d+).profile" \
        --sample-regex ".*/(.*?)/NCBI/(.*?)/.*?(\d+)[_.].*" \
        --test-sample-regex


merge-ncbi:
    #!/usr/bin/env bash
    set -euo pipefail
    shopt -s nullglob
    cargo run -- merge \
        --input data/test_data/profiles/predictions/*/NCBI/{*,*/*,*/*/*}.{tsv,bracken,profile} \
        --output data/test_data/profiles/predictions/matrix/ncbi.tsv \
        --sample-regex ".*/(.*?)/NCBI/human/(.*?)/(\d+)\..*" \
        --sample-regex ".*/(.*?)/NCBI/(.*?)/.*?_(\d+)_.*" \
        --sample-regex ".*/(.*?)/NCBI/(.*?)_?(\d+).profile" \
        --sample-regex ".*/(.*?)/NCBI/(.*?)/.*?(\d+)[_.].*" \
        --cami-ignore-lineage-error


msa-test:
    #!/usr/bin/env bash
    set -euo pipefail
    cargo run -- msa \
        --meta data/test_data/strain/msa/test/meta.xlsx \
        --output data/test_data/strain/msa/test/stats

strain-test1:
    #!/usr/bin/env bash
    set -euo pipefail
    cargo run -- strain \
        --meta data/test_data/strain/test1/samplesheet/samplesheet.tsv \
        --outprefix data/test_data/strain/test1/output/test1

plot-strain-test1:
    Rscript scripts/benchpro_visualize_monophyly.R \
        data/test_data/strain/test1/output/test1.genome_proximity.tsv \
        data/test_data/strain/test1/output/test1.monophyly_plot.pdf
