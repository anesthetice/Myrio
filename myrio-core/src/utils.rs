// Imports
use console::StyledObject;
use indicatif::{MultiProgress, ProgressBar};

use crate::data::Float;

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

#[inline]
pub fn extract_singlevec<T>(vec: Vec<T>) -> T {
    if let Ok([val]) = <[_; 1]>::try_from(vec) {
        val
    } else {
        panic!("Expected vector containing a single element");
    }
}

pub fn colorizer_from_confidence_score<T>(
    conf: Float
) -> fn(console::StyledObject<T>) -> console::StyledObject<T> {
    match conf {
        0.0..0.25 => StyledObject::<T>::red,
        0.25..0.5 => StyledObject::<T>::yellow,
        0.5..0.75 => StyledObject::<T>::green,
        0.75..1.0 => StyledObject::<T>::cyan,
        _ => {
            eprintln!("Warning: got unexpected confidence score of {conf}");
            StyledObject::<T>::bold
        }
    }
}
