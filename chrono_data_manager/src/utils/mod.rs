/// Extract GPS dataset identifier from dataset name (e.g., "gbm19670" -> 19670)
pub fn extract_gps_dataset_id(dataset: &str) -> Option<u32> {
    // The dataset string is typically just the basename without extension, e.g., "gbm19670"
    // "gbm" (3 chars) + Week (4 chars) + Day (1 char) = 5-digit identifier
    if dataset.len() >= 8 && dataset.starts_with("gbm") {
        let id_str = &dataset[3..8];
        id_str.parse::<u32>().ok()
    } else {
        None
    }
}
