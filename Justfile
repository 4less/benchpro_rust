install-local:
    cargo install --path /home/fritscher/git/4less/benchpro_rust --root ~/.local

sylph-debug:
    mkdir -p sandbox/output
    cargo run -- --meta /home/fritscher/git/4less/protal_paper/data/input/meta_debug_sylph_mouse.xlsx --outprefix sandbox/output/sylph_debug
