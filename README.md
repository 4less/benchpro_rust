# Benchpro

Rust rewrite of 4less/benchpro, benchmarking tool for taxonomic profilers.

# Requirements
Make sure rust is installed.

# Run test
Clone the repository, change directories into the project and run cargo run. The executable is located at target/release/benchpro.

```
git clone
cd benchpro_rust
cargo run -- --meta data/test_data/meta/meta.xlsx --outprefix data/test_data/output/test
```
