// Modules
pub mod app;

// Imports
use crate::app::App;

fn main() -> anyhow::Result<()> {
    let app = App::load()?;
    app.run()
}
