use std::{env, io::Error};

use clap::ValueEnum;
use clap_complete::{Shell, generate_to};

include!("src/app/cli.rs");

fn main() -> Result<(), Error> {
    let outdir = match env::var_os("OUT_DIR") {
        Some(outdir) => outdir,
        None => return Ok(()),
    };

    let mut cmd = build_cli();
    for &shell in Shell::value_variants() {
        generate_to(shell, &mut cmd, "myrio", &outdir)?;
    }

    Ok(())
}
