mod data_loader;
mod data_manager;
mod data_path;
mod errors;
mod types;
mod utils;

// data loader
pub(crate) use crate::data_loader::load_all_satellites::load_all_satellites;
pub(crate) use crate::data_loader::load_clk::parse_clk as load_clock_data;
pub(crate) use crate::data_loader::load_data::load_data;
pub(crate) use crate::data_loader::load_g_values::load_gm_values_from_csv;
pub(crate) use crate::data_loader::load_inferometer::{
    get_zenodo_data_path, load_arb_rot_data, load_einstein_elevator_data,
};
pub(crate) use crate::data_loader::load_mhd::load_omni_data;
pub(crate) use crate::data_loader::load_pta::load_all_pta;
pub(crate) use crate::data_loader::load_pta::{get_ipta_data_path, get_nanograv_data_path};
pub(crate) use crate::data_loader::load_pulsar::discover_binary_systems;
pub(crate) use crate::data_loader::load_sp3::parse_sp3 as load_orbit_data;
pub(crate) use crate::data_loader::load_sparc::{
    list_galaxies, load_all_rotation_curves, load_galaxy_by_name,
};

// data manger
pub use crate::data_manager::DataManager;
// data paths
pub use crate::data_path::*;

// Errors
pub use crate::errors::conversion_error::ConversionError;

// Types
pub use crate::types::clock_types::ClockData;
pub use crate::types::gnss_types::GnssDataResult;
pub use crate::types::interferometry_types::*;
pub use crate::types::mhd_types::*;
pub use crate::types::observatory_types::{FleetData, ObservatoryEpoch, SatelliteState};
pub use crate::types::orbit_types::OrbitData;
pub use crate::types::pta_types::*;
pub use crate::types::satelite_types::SatId;
pub use crate::types::sparc_types::*;
pub use crate::types::strong_gravity_types::*;

// Utils
pub use crate::utils::extract_gps_dataset_id;

// 2017 t0 2018 uses the new high precision IGS14 Reference Frame
pub const YEARS: [&str; 2] = ["2017", "2018"];

// 2016 usess the OLD low precision IGNS08  Reference Frame and thus
// reduces precison of multi year calculations by at least two orders of magntide.
// Therefore, it is strongly recomneced to use the YEARS constant above to exclusde 2016.
pub const YEARS_INCL_2016: [&str; 3] = ["2016", "2017", "2018"];

/// Anomalous GPS datasets to exclude from analysis
/// Format: 5-digit identifier = WWWWD (Week + Day)
pub const ANOMALOUS_WEEKS: &[u32] = &[
    //  2016
    // E14-specific issues
    19236, 19095, 1923, 1909,
    // 2017
    // IGS14 Reference Frame Transition (Jan 1 - Feb 25, 2017)
    // IGS switched from IGS08 to IGS14 on Jan 29, 2017 (Week 1934)
    // Excluding full Jan/Feb period to avoid transition artifacts
    19300, 19301, 19302, 19303, 19304, 19305, 19306, // Week 1930 (Jan 01-07)
    19310, 19311, 19312, 19313, 19314, 19315, 19316, // Week 1931 (Jan 08-14)
    19320, 19321, 19322, 19323, 19324, 19325, 19326, // Week 1932 (Jan 15-21)
    19330, 19331, 19332, 19333, 19334, 19335, 19336, // Week 1933 (Jan 22-28)
    19340, 19341, 19342, 19343, 19344, 19345, 19346, // Week 1934 (Jan 29-Feb 04) - TRANSITION
    19350, 19351, 19352, 19353, 19354, 19355, 19356, // Week 1935 (Feb 05-11)
    19360, 19361, 19362, 19363, 19364, 19365, 19366, // Week 1936 (Feb 12-18)
    19370, 19371, 19372, 19373, 19374, 19375, 19376, // Week 1937 (Feb 19-25)
    // E14-specific anomalies only
    19625, // Aug 18, 2017 (E14 EXTREME anomaly)
    // 2018
    // E14/E18 Extreme Outliers (Full September Exclusion)
    // Week 1999 (Apr 29 - May 05, 2018) - E18 Extreme Anomaly
    19990, 19991, 19992, 19993, 19994, 19995, 19996,
    // September 2018 "Crisis" - Filtering ALL of September
    20176, // Week 2017 (Partial: Sep 1)
    20180, 20181, 20182, 20183, 20184, 20185, 20186, // Week 2018 (Sep 02 - Sep 08)
    20190, 20191, 20192, 20193, 20194, 20195, 20196, // Week 2019 (Sep 09 - Sep 15)
    20200, 20201, 20202, 20203, 20204, 20205, 20206, // Week 2020 (Sep 16 - Sep 22)
    20210, 20211, 20212, 20213, 20214, 20215, 20216, // Week 2021 (Sep 23 - Sep 29)
    20220, 20221, 20222, 20223, 20224, 20225, 20226, // Week 2022 - Rank 1 Anomaly
];
