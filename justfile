fmt:
    cargo +nightly fmt

runmaxcli +ARGS:
    RUSTFLAGS="-C target-cpu=native -Awarnings" cargo run -p myrio-cli --release -- {{ARGS}}

runmaxexp:
    RUSTFLAGS="-C target-cpu=native -Awarnings" cargo run -p myrio-exp --release

test:
    RUSTFLAGS="-Awarnings" cargo test -p myrio-cli --no-default-features
    RUSTFLAGS="-Awarnings" cargo test -p myrio-exp --no-default-features
    RUSTFLAGS="-Awarnings" cargo test -p myrio-core --no-default-features

flamegraph:
    RUSTFLAGS="-C force-frame-pointers=yes" cargo flamegraph --package myrio-exp --profile profiling

regen-myrio-db:
    just runmaxcli tree new -i ./ignore/myrio-db/BOLD_Plantae_20250831_ITS.fasta -g "ITS"
    just runmaxcli tree new -i ./ignore/myrio-db/BOLD_Plantae_20250831_matK.fasta -g "matK"
    just runmaxcli tree new -i ./ignore/myrio-db/BOLD_Plantae_20250831_rbcL.fasta -g "rbcL"
    just runmaxcli tree new -i ./ignore/myrio-db/BOLD_Plantae_20250831_trnH-psbA.fasta -g "trnH-psbA"
