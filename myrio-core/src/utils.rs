#[cfg(feature = "indicatif")]
pub(crate) fn greenify(s: &str) -> String {
    console::style(s).green().to_string()
}

#[cfg(feature = "indicatif")]
pub(crate) fn simple_spinner(
    start_message: Option<String>,
    steady_tick_ms: Option<u64>,
) -> indicatif::ProgressBar {
    let spinner = indicatif::ProgressBar::new_spinner().with_style(
        indicatif::ProgressStyle::default_spinner()
            .tick_chars("⊶⊷✔")
            .template("{spinner:.green} [{elapsed_precise}] {msg}")
            .unwrap(),
    );
    if let Some(ms) = steady_tick_ms {
        spinner.enable_steady_tick(std::time::Duration::from_millis(ms));
    }
    if let Some(msg) = start_message {
        spinner.set_message(msg);
    }
    spinner
}

#[cfg(feature = "indicatif")]
pub(crate) fn simple_progressbar(
    len: usize,
    text: impl std::fmt::Display,
) -> indicatif::ProgressBar {
    indicatif::ProgressBar::new(len as u64)
        .with_style(
            indicatif::ProgressStyle::with_template(&format!(
                "{{msg}} [{{elapsed_precise}}] {{bar:40.cyan/blue}} {{pos}} {text}"
            ))
            .unwrap(),
        )
        .with_message("⋆")
        .with_finish(indicatif::ProgressFinish::WithMessage(greenify("✔").into()))
}
