fmt:
    cargo +nightly fmt

runmax:
    RUSTFLAGS="-C target-cpu=native -Awarnings" cargo run --release

test:
    RUSTFLAGS="-Awarnings" cargo test -p myrio-cli --no-default-features
    RUSTFLAGS="-Awarnings" cargo test -p myrio-exp --no-default-features
    RUSTFLAGS="-Awarnings" cargo test -p myrio-core --no-default-features

flamegraph:
    RUSTFLAGS="-C force-frame-pointers=yes" cargo flamegraph --package myrio-exp --profile profiling
