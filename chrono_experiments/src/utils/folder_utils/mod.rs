use std::path::Path;
use std::{env, io};

pub fn check_folder(folder_path: &str) -> Result<(), io::Error> {
    let path = Path::new(folder_path);
    if path.exists() && path.is_dir() {
        println!("Folder found at path: {}", path.display());
    } else {
        eprintln!("No folder found at: {}", path.display());
    }

    Ok(())
}

pub fn check_working_dir() {
    match env::current_dir() {
        Ok(path) => println!("Current working directory: {}", path.display()),
        Err(e) => eprintln!("Failed to get current directory: {}", e),
    }
}
