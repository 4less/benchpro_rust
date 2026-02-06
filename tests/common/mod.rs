use std::{
    env, fs,
    path::{Path, PathBuf},
    process::Command,
    time::{SystemTime, UNIX_EPOCH},
};

pub fn unique_temp_dir(prefix: &str) -> PathBuf {
    let base = env::temp_dir();
    let timestamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("SystemTime before UNIX_EPOCH")
        .as_millis();
    let pid = std::process::id();
    base.join(format!("{}_{}_{}", prefix, pid, timestamp))
}

pub fn read_file(path: &Path) -> Vec<u8> {
    fs::read(path).unwrap_or_else(|err| {
        panic!("Failed to read {}: {}", path.display(), err);
    })
}

pub fn benchpro_bin() -> PathBuf {
    if let Ok(path) = env::var("CARGO_BIN_EXE_benchpro") {
        return PathBuf::from(path);
    }

    let mut path = PathBuf::from("target/debug/benchpro");
    if cfg!(windows) {
        path.set_extension("exe");
    }

    if path.exists() {
        return path;
    }

    let status = Command::new("cargo")
        .args(["build", "--bin", "benchpro"])
        .status()
        .expect("Failed to build benchpro binary");
    assert!(status.success(), "benchpro build failed");

    if !path.exists() {
        panic!("Benchpro binary not found at {}", path.display());
    }

    path
}
