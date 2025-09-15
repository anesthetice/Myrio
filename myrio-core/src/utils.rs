// Imports
use indicatif::{MultiProgress, ProgressBar};

pub fn greenify(s: &str) -> String {
    console::style(s).green().to_string()
}

pub fn simple_spinner<S: ToString>(
    start_message: Option<S>,
    steady_tick_ms: Option<u64>,
    multi: Option<&MultiProgress>,
) -> ProgressBar {
    #[cfg(feature = "progress")]
    {
        let mut spinner = indicatif::ProgressBar::new_spinner().with_style(
            indicatif::ProgressStyle::default_spinner()
                .tick_chars("⊶⊷✔")
                .template("{spinner:.green} [{elapsed_precise}] {msg}")
                .unwrap(),
        );
        if let Some(ms) = steady_tick_ms {
            spinner.enable_steady_tick(std::time::Duration::from_millis(ms));
        }
        if let Some(msg) = start_message {
            spinner.set_message(msg.to_string());
        }
        if let Some(multi) = multi {
            spinner = multi.add(spinner);
            multi.set_move_cursor(true);
        }
        spinner
    }

    #[cfg(not(feature = "progress"))]
    ProgressBar::hidden()
}

pub fn simple_progressbar(
    len: usize,
    text: impl std::fmt::Display,
    multi: Option<&MultiProgress>,
) -> indicatif::ProgressBar {
    #[cfg(feature = "progress")]
    {
        let mut pb = indicatif::ProgressBar::new(len as u64)
            .with_style(
                indicatif::ProgressStyle::with_template(&format!(
                    "{{msg}} [{{elapsed_precise}}] {{bar:40.cyan/blue}} {{pos}} {text}"
                ))
                .unwrap(),
            )
            .with_message("⋆")
            .with_finish(indicatif::ProgressFinish::WithMessage(greenify("✔").into()));
        if let Some(multi) = multi {
            pb = multi.add(pb);
            multi.set_move_cursor(true);
        }
        pb
    }

    #[cfg(not(feature = "progress"))]
    ProgressBar::hidden()
}
