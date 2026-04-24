set shell := ["bash", "-cu"]

install-local:
    cargo install --path . --root ~/.local
    install -m 755 scripts/benchpro_visualize_monophyly.R ~/.local/bin/benchpro_visualize_monophyly

sylph-debug:
    mkdir -p sandbox/output
    cargo run -- --meta ./data/input/meta_debug_sylph_mouse.xlsx --outprefix sandbox/output/sylph_debug
    
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
