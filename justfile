fmt:
    cargo +nightly fmt

runmax:
    RUSTFLAGS="-C target-cpu=native" cargo run --release

flamegraph:
    RUSTFLAGS="-C force-frame-pointers=yes" cargo flamegraph --package myrio-exp --profile profiling
