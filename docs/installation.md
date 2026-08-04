# Installation

## With `just`

```bash
just install          # builds release and installs into $CARGO_HOME/bin (default ~/.cargo/bin)
just install-local    # installs into ~/.local/bin instead
just prefix=/usr/local install
just uninstall
```

`just install` also installs `benchpro_visualize_monophyly` (the R helper used by
`benchpro strain` plots) when `scripts/` is present, and warns if the target `bin` directory is not
on your `PATH`. `just` on its own lists every recipe; `just check` runs fmt, clippy and the tests as
CI would.


This is a CLI tool, so the goal is to build a binary and put it on your PATH. No admin rights are required.

## Option A: Build in place with Cargo

From the repo root:

```bash
cargo build --release
```

The binary will be at:

```
target/release/benchpro
```

Run it directly:

```bash
./target/release/benchpro --help
```

To make it available on your PATH without admin rights, add a user-local bin directory and symlink:

```bash
mkdir -p "$HOME/.local/bin"
ln -s "$PWD/target/release/benchpro" "$HOME/.local/bin/benchpro"
```

Then ensure `~/.local/bin` is on your PATH (add this to your shell profile if needed):

```bash
export PATH="$HOME/.local/bin:$PATH"
```

## Option B: Cargo install to a user-local prefix

You can install the binary into a user-local location:

```bash
cargo install --path . --root "$HOME/.local"
```

This installs to:

```
$HOME/.local/bin/benchpro
```

Again, ensure `~/.local/bin` is on your PATH.
