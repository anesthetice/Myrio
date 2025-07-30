fmt:
    cargo +nightly fmt

runmax:
    RUSTFLAGS="-C target-cpu=native" cargo run --release
