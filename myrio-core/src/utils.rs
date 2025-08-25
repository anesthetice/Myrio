#[cfg(feature = "indicatif")]
pub(crate) fn greenify(s: &str) -> String {
    console::style(s).green().to_string()
}
