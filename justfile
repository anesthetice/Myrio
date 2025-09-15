fmt:
    cargo +nightly fmt

install:
    RUSTFLAGS="-C target-cpu=native" cargo install --path myrio-cli

runmaxcli +ARGS:
    RUSTFLAGS="-C target-cpu=native -Awarnings" cargo run -p myrio-cli --release -- {{ARGS}}

rundbgcli +ARGS:
    RUSTFLAGS="-Awarnings" RUST_BACKTRACE=1 cargo run -p myrio-cli -- {{ARGS}}

runmaxexp:
    RUSTFLAGS="-C target-cpu=native -Awarnings" cargo run -p myrio-exp --release

test:
    #RUSTFLAGS="-Awarnings" cargo test -p myrio-cli --no-default-features
    #RUSTFLAGS="-Awarnings" cargo test -p myrio-exp --no-default-features
    RUSTFLAGS="-Awarnings" cargo test -p myrio-core --no-default-features

flamegraph:
    RUSTFLAGS="-C force-frame-pointers=yes" cargo flamegraph --package myrio-exp --profile profiling

flamegraph-alt:
    RUSTFLAGS="-C force-frame-pointers=yes" cargo flamegraph --package myrio-cli --profile profiling -- run -i ignore/queries/Berberis_Julianae_Campus_matk_rbcL_psbA-trnH_ITS_barcode6.fastq -t ignore/myrio-db/


regen-myrio-db:
    just runmaxcli tree new -i ./ignore/myrio-db/BOLD_Plantae_20250831_ITS.fasta -g "ITS"
    just runmaxcli tree new -i ./ignore/myrio-db/BOLD_Plantae_20250831_matK.fasta -g "matK"
    just runmaxcli tree new -i ./ignore/myrio-db/BOLD_Plantae_20250831_rbcL.fasta -g "rbcL"
    just runmaxcli tree new -i ./ignore/myrio-db/BOLD_Plantae_20250831_trnH-psbA.fasta -g "trnH-psbA"

regen-myrio-db-with-precomputation:
    just runmaxcli tree new -i ./ignore/myrio-db/BOLD_Plantae_20250831_ITS.fasta -g "ITS" -k 18
    just runmaxcli tree new -i ./ignore/myrio-db/BOLD_Plantae_20250831_matK.fasta -g "matK" -k 18
    just runmaxcli tree new -i ./ignore/myrio-db/BOLD_Plantae_20250831_rbcL.fasta -g "rbcL" -k 18
    just runmaxcli tree new -i ./ignore/myrio-db/BOLD_Plantae_20250831_trnH-psbA.fasta -g "trnH-psbA" -k 18
