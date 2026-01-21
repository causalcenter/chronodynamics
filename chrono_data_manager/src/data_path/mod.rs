/*
 * SPDX-License-Identifier: MIT
 * Copyright (c) "2025" . The DeepCausality Authors and Contributors. All Rights Reserved.
 */

use crate::YEARS;
use std::path::PathBuf;

/// Get the absolute path to the data input directory (data/gnss)
pub fn get_gnss_data_input_path() -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("data/gnss");
    path
}

/// Get the absolute path to the data output directory for a specific subdir
pub fn get_data_output_path(subdir: &str) -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("data/output");
    path.push(subdir);
    path
}

/// Get the absolute path to the SPARC data directory (data/zenodo/sparc)
pub fn get_sparc_data_path() -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("data/zenodo/sparc");
    path
}

/// Get datasets for a year
pub fn get_year_datasets(year: &str) -> Vec<&'static str> {
    match year {
        "2016" => vec![
            // 2016 has 371 datasets (53 weeks * 7 days)
            "gbm18770", "gbm18771", "gbm18772", "gbm18773", "gbm18774", "gbm18775", "gbm18776",
            "gbm18780", "gbm18781", "gbm18782", "gbm18783", "gbm18784", "gbm18785", "gbm18786",
            "gbm18790", "gbm18791", "gbm18792", "gbm18793", "gbm18794", "gbm18795", "gbm18796",
            "gbm18800", "gbm18801", "gbm18802", "gbm18803", "gbm18804", "gbm18805", "gbm18806",
            "gbm18810", "gbm18811", "gbm18812", "gbm18813", "gbm18814", "gbm18815", "gbm18816",
            "gbm18820", "gbm18821", "gbm18822", "gbm18823", "gbm18824", "gbm18825", "gbm18826",
            "gbm18830", "gbm18831", "gbm18832", "gbm18833", "gbm18834", "gbm18835", "gbm18836",
            "gbm18840", "gbm18841", "gbm18842", "gbm18843", "gbm18844", "gbm18845", "gbm18846",
            "gbm18850", "gbm18851", "gbm18852", "gbm18853", "gbm18854", "gbm18855", "gbm18856",
            "gbm18860", "gbm18861", "gbm18862", "gbm18863", "gbm18864", "gbm18865", "gbm18866",
            "gbm18870", "gbm18871", "gbm18872", "gbm18873", "gbm18874", "gbm18875", "gbm18876",
            "gbm18880", "gbm18881", "gbm18882", "gbm18883", "gbm18884", "gbm18885", "gbm18886",
            "gbm18890", "gbm18891", "gbm18892", "gbm18893", "gbm18894", "gbm18895", "gbm18896",
            "gbm18900", "gbm18901", "gbm18902", "gbm18903", "gbm18904", "gbm18905", "gbm18906",
            "gbm18910", "gbm18911", "gbm18912", "gbm18913", "gbm18914", "gbm18915", "gbm18916",
            "gbm18920", "gbm18921", "gbm18922", "gbm18923", "gbm18924", "gbm18925", "gbm18926",
            "gbm18930", "gbm18931", "gbm18932", "gbm18933", "gbm18934", "gbm18935", "gbm18936",
            "gbm18940", "gbm18941", "gbm18942", "gbm18943", "gbm18944", "gbm18945", "gbm18946",
            "gbm18950", "gbm18951", "gbm18952", "gbm18953", "gbm18954", "gbm18955", "gbm18956",
            "gbm18960", "gbm18961", "gbm18962", "gbm18963", "gbm18964", "gbm18965", "gbm18966",
            "gbm18970", "gbm18971", "gbm18972", "gbm18973", "gbm18974", "gbm18975", "gbm18976",
            "gbm18980", "gbm18981", "gbm18982", "gbm18983", "gbm18984", "gbm18985", "gbm18986",
            "gbm18990", "gbm18991", "gbm18992", "gbm18993", "gbm18994", "gbm18995", "gbm18996",
            "gbm19000", "gbm19001", "gbm19002", "gbm19003", "gbm19004", "gbm19005", "gbm19006",
            "gbm19010", "gbm19011", "gbm19012", "gbm19013", "gbm19014", "gbm19015", "gbm19016",
            "gbm19020", "gbm19021", "gbm19022", "gbm19023", "gbm19024", "gbm19025", "gbm19026",
            "gbm19030", "gbm19031", "gbm19032", "gbm19033", "gbm19034", "gbm19035", "gbm19036",
            "gbm19040", "gbm19041", "gbm19042", "gbm19043", "gbm19044", "gbm19045", "gbm19046",
            "gbm19050", "gbm19051", "gbm19052", "gbm19053", "gbm19054", "gbm19055", "gbm19056",
            "gbm19060", "gbm19061", "gbm19062", "gbm19063", "gbm19064", "gbm19065", "gbm19066",
            "gbm19070", "gbm19071", "gbm19072", "gbm19073", "gbm19074", "gbm19075", "gbm19076",
            "gbm19080", "gbm19081", "gbm19082", "gbm19083", "gbm19084", "gbm19085", "gbm19086",
            "gbm19090", "gbm19091", "gbm19092", "gbm19093", "gbm19094", "gbm19095", "gbm19096",
            "gbm19100", "gbm19101", "gbm19102", "gbm19103", "gbm19104", "gbm19105", "gbm19106",
            "gbm19110", "gbm19111", "gbm19112", "gbm19113", "gbm19114", "gbm19115", "gbm19116",
            "gbm19120", "gbm19121", "gbm19122", "gbm19123", "gbm19124", "gbm19125", "gbm19126",
            "gbm19130", "gbm19131", "gbm19132", "gbm19133", "gbm19134", "gbm19135", "gbm19136",
            "gbm19140", "gbm19141", "gbm19142", "gbm19143", "gbm19144", "gbm19145", "gbm19146",
            "gbm19150", "gbm19151", "gbm19152", "gbm19153", "gbm19154", "gbm19155", "gbm19156",
            "gbm19160", "gbm19161", "gbm19162", "gbm19163", "gbm19164", "gbm19165", "gbm19166",
            "gbm19170", "gbm19171", "gbm19172", "gbm19173", "gbm19174", "gbm19175", "gbm19176",
            "gbm19180", "gbm19181", "gbm19182", "gbm19183", "gbm19184", "gbm19185", "gbm19186",
            "gbm19190", "gbm19191", "gbm19192", "gbm19193", "gbm19194", "gbm19195", "gbm19196",
            "gbm19200", "gbm19201", "gbm19202", "gbm19203", "gbm19204", "gbm19205", "gbm19206",
            "gbm19210", "gbm19211", "gbm19212", "gbm19213", "gbm19214", "gbm19215", "gbm19216",
            "gbm19220", "gbm19221", "gbm19222", "gbm19223", "gbm19224", "gbm19225", "gbm19226",
            "gbm19230", "gbm19231", "gbm19232", "gbm19233", "gbm19234", "gbm19235", "gbm19236",
            "gbm19240", "gbm19241", "gbm19242", "gbm19243", "gbm19244", "gbm19245", "gbm19246",
            "gbm19250", "gbm19251", "gbm19252", "gbm19253", "gbm19254", "gbm19255", "gbm19256",
            "gbm19260", "gbm19261", "gbm19262", "gbm19263", "gbm19264", "gbm19265", "gbm19266",
            "gbm19270", "gbm19271", "gbm19272", "gbm19273", "gbm19274", "gbm19275", "gbm19276",
            "gbm19280", "gbm19281", "gbm19282", "gbm19283", "gbm19284", "gbm19285", "gbm19286",
            "gbm19290", "gbm19291", "gbm19292", "gbm19293", "gbm19294", "gbm19295", "gbm19296",
        ],
        "2017" => vec![
            // 2017 has 364 datasets (52 weeks * 7 days)
            "gbm19300", "gbm19301", "gbm19302", "gbm19303", "gbm19304", "gbm19305", "gbm19306",
            "gbm19310", "gbm19311", "gbm19312", "gbm19313", "gbm19314", "gbm19315", "gbm19316",
            "gbm19320", "gbm19321", "gbm19322", "gbm19323", "gbm19324", "gbm19325", "gbm19326",
            "gbm19330", "gbm19331", "gbm19332", "gbm19333", "gbm19334", "gbm19335", "gbm19336",
            "gbm19340", "gbm19341", "gbm19342", "gbm19343", "gbm19344", "gbm19345", "gbm19346",
            "gbm19350", "gbm19351", "gbm19352", "gbm19353", "gbm19354", "gbm19355", "gbm19356",
            "gbm19360", "gbm19361", "gbm19362", "gbm19363", "gbm19364", "gbm19365", "gbm19366",
            "gbm19370", "gbm19371", "gbm19372", "gbm19373", "gbm19374", "gbm19375", "gbm19376",
            "gbm19380", "gbm19381", "gbm19382", "gbm19383", "gbm19384", "gbm19385", "gbm19386",
            "gbm19390", "gbm19391", "gbm19392", "gbm19393", "gbm19394", "gbm19395", "gbm19396",
            "gbm19400", "gbm19401", "gbm19402", "gbm19403", "gbm19404", "gbm19405", "gbm19406",
            "gbm19410", "gbm19411", "gbm19412", "gbm19413", "gbm19414", "gbm19415", "gbm19416",
            "gbm19420", "gbm19421", "gbm19422", "gbm19423", "gbm19424", "gbm19425", "gbm19426",
            "gbm19430", "gbm19431", "gbm19432", "gbm19433", "gbm19434", "gbm19435", "gbm19436",
            "gbm19440", "gbm19441", "gbm19442", "gbm19443", "gbm19444", "gbm19445", "gbm19446",
            "gbm19450", "gbm19451", "gbm19452", "gbm19453", "gbm19454", "gbm19455", "gbm19456",
            "gbm19460", "gbm19461", "gbm19462", "gbm19463", "gbm19464", "gbm19465", "gbm19466",
            "gbm19470", "gbm19471", "gbm19472", "gbm19473", "gbm19474", "gbm19475", "gbm19476",
            "gbm19480", "gbm19481", "gbm19482", "gbm19483", "gbm19484", "gbm19485", "gbm19486",
            "gbm19490", "gbm19491", "gbm19492", "gbm19493", "gbm19494", "gbm19495", "gbm19496",
            "gbm19500", "gbm19501", "gbm19502", "gbm19503", "gbm19504", "gbm19505", "gbm19506",
            "gbm19510", "gbm19511", "gbm19512", "gbm19513", "gbm19514", "gbm19515", "gbm19516",
            "gbm19520", "gbm19521", "gbm19522", "gbm19523", "gbm19524", "gbm19525", "gbm19526",
            "gbm19530", "gbm19531", "gbm19532", "gbm19533", "gbm19534", "gbm19535", "gbm19536",
            "gbm19540", "gbm19541", "gbm19542", "gbm19543", "gbm19544", "gbm19545", "gbm19546",
            "gbm19550", "gbm19551", "gbm19552", "gbm19553", "gbm19554", "gbm19555", "gbm19556",
            "gbm19560", "gbm19561", "gbm19562", "gbm19563", "gbm19564", "gbm19565", "gbm19566",
            "gbm19570", "gbm19571", "gbm19572", "gbm19573", "gbm19574", "gbm19575", "gbm19576",
            "gbm19580", "gbm19581", "gbm19582", "gbm19583", "gbm19584", "gbm19585", "gbm19586",
            "gbm19590", "gbm19591", "gbm19592", "gbm19593", "gbm19594", "gbm19595", "gbm19596",
            "gbm19600", "gbm19601", "gbm19602", "gbm19603", "gbm19604", "gbm19605", "gbm19606",
            "gbm19610", "gbm19611", "gbm19612", "gbm19613", "gbm19614", "gbm19615", "gbm19616",
            "gbm19620", "gbm19621", "gbm19622", "gbm19623", "gbm19624", "gbm19625", "gbm19626",
            "gbm19630", "gbm19631", "gbm19632", "gbm19633", "gbm19634", "gbm19635", "gbm19636",
            "gbm19640", "gbm19641", "gbm19642", "gbm19643", "gbm19644", "gbm19645", "gbm19646",
            "gbm19650", "gbm19651", "gbm19652", "gbm19653", "gbm19654", "gbm19655", "gbm19656",
            "gbm19660", "gbm19661", "gbm19662", "gbm19663", "gbm19664", "gbm19665", "gbm19666",
            "gbm19670", "gbm19671", "gbm19672", "gbm19673", "gbm19674", "gbm19675", "gbm19676",
            "gbm19680", "gbm19681", "gbm19682", "gbm19683", "gbm19684", "gbm19685", "gbm19686",
            "gbm19690", "gbm19691", "gbm19692", "gbm19693", "gbm19694", "gbm19695", "gbm19696",
            "gbm19700", "gbm19701", "gbm19702", "gbm19703", "gbm19704", "gbm19705", "gbm19706",
            "gbm19710", "gbm19711", "gbm19712", "gbm19713", "gbm19714", "gbm19715", "gbm19716",
            "gbm19720", "gbm19721", "gbm19722", "gbm19723", "gbm19724", "gbm19725", "gbm19726",
            "gbm19730", "gbm19731", "gbm19732", "gbm19733", "gbm19734", "gbm19735", "gbm19736",
            "gbm19740", "gbm19741", "gbm19742", "gbm19743", "gbm19744", "gbm19745", "gbm19746",
            "gbm19750", "gbm19751", "gbm19752", "gbm19753", "gbm19754", "gbm19755", "gbm19756",
            "gbm19760", "gbm19761", "gbm19762", "gbm19763", "gbm19764", "gbm19765", "gbm19766",
            "gbm19770", "gbm19771", "gbm19772", "gbm19773", "gbm19774", "gbm19775", "gbm19776",
            "gbm19780", "gbm19781", "gbm19782", "gbm19783", "gbm19784", "gbm19785", "gbm19786",
            "gbm19790", "gbm19791", "gbm19792", "gbm19793", "gbm19794", "gbm19795", "gbm19796",
            "gbm19800", "gbm19801", "gbm19802", "gbm19803", "gbm19804", "gbm19805", "gbm19806",
            "gbm19810", "gbm19811", "gbm19812", "gbm19813", "gbm19814", "gbm19815", "gbm19816",
        ],
        "2018" => vec![
            // 2018 has 371 datasets (53 weeks * 7 days)
            "gbm19820", "gbm19821", "gbm19822", "gbm19823", "gbm19824", "gbm19825", "gbm19826",
            "gbm19830", "gbm19831", "gbm19832", "gbm19833", "gbm19834", "gbm19835", "gbm19836",
            "gbm19840", "gbm19841", "gbm19842", "gbm19843", "gbm19844", "gbm19845", "gbm19846",
            "gbm19850", "gbm19851", "gbm19852", "gbm19853", "gbm19854", "gbm19855", "gbm19856",
            "gbm19860", "gbm19861", "gbm19862", "gbm19863", "gbm19864", "gbm19865", "gbm19866",
            "gbm19870", "gbm19871", "gbm19872", "gbm19873", "gbm19874", "gbm19875", "gbm19876",
            "gbm19880", "gbm19881", "gbm19882", "gbm19883", "gbm19884", "gbm19885", "gbm19886",
            "gbm19890", "gbm19891", "gbm19892", "gbm19893", "gbm19894", "gbm19895", "gbm19896",
            "gbm19900", "gbm19901", "gbm19902", "gbm19903", "gbm19904", "gbm19905", "gbm19906",
            "gbm19910", "gbm19911", "gbm19912", "gbm19913", "gbm19914", "gbm19915", "gbm19916",
            "gbm19920", "gbm19921", "gbm19922", "gbm19923", "gbm19924", "gbm19925", "gbm19926",
            "gbm19930", "gbm19931", "gbm19932", "gbm19933", "gbm19934", "gbm19935", "gbm19936",
            "gbm19940", "gbm19941", "gbm19942", "gbm19943", "gbm19944", "gbm19945", "gbm19946",
            "gbm19950", "gbm19951", "gbm19952", "gbm19953", "gbm19954", "gbm19955", "gbm19956",
            "gbm19960", "gbm19961", "gbm19962", "gbm19963", "gbm19964", "gbm19965", "gbm19966",
            "gbm19970", "gbm19971", "gbm19972", "gbm19973", "gbm19974", "gbm19975", "gbm19976",
            "gbm19980", "gbm19981", "gbm19982", "gbm19983", "gbm19984", "gbm19985", "gbm19986",
            "gbm19990", "gbm19991", "gbm19992", "gbm19993", "gbm19994", "gbm19995", "gbm19996",
            "gbm20000", "gbm20001", "gbm20002", "gbm20003", "gbm20004", "gbm20005", "gbm20006",
            "gbm20010", "gbm20011", "gbm20012", "gbm20013", "gbm20014", "gbm20015", "gbm20016",
            "gbm20020", "gbm20021", "gbm20022", "gbm20023", "gbm20024", "gbm20025", "gbm20026",
            "gbm20030", "gbm20031", "gbm20032", "gbm20033", "gbm20034", "gbm20035", "gbm20036",
            "gbm20040", "gbm20041", "gbm20042", "gbm20043", "gbm20044", "gbm20045", "gbm20046",
            "gbm20050", "gbm20051", "gbm20052", "gbm20053", "gbm20054", "gbm20055", "gbm20056",
            "gbm20060", "gbm20061", "gbm20062", "gbm20063", "gbm20064", "gbm20065", "gbm20066",
            "gbm20070", "gbm20071", "gbm20072", "gbm20073", "gbm20074", "gbm20075", "gbm20076",
            "gbm20080", "gbm20081", "gbm20082", "gbm20083", "gbm20084", "gbm20085", "gbm20086",
            "gbm20090", "gbm20091", "gbm20092", "gbm20093", "gbm20094", "gbm20095", "gbm20096",
            "gbm20100", "gbm20101", "gbm20102", "gbm20103", "gbm20104", "gbm20105", "gbm20106",
            "gbm20110", "gbm20111", "gbm20112", "gbm20113", "gbm20114", "gbm20115", "gbm20116",
            "gbm20120", "gbm20121", "gbm20122", "gbm20123", "gbm20124", "gbm20125", "gbm20126",
            "gbm20130", "gbm20131", "gbm20132", "gbm20133", "gbm20134", "gbm20135", "gbm20136",
            "gbm20140", "gbm20141", "gbm20142", "gbm20143", "gbm20144", "gbm20145", "gbm20146",
            "gbm20150", "gbm20151", "gbm20152", "gbm20153", "gbm20154", "gbm20155", "gbm20156",
            "gbm20160", "gbm20161", "gbm20162", "gbm20163", "gbm20164", "gbm20165", "gbm20166",
            "gbm20170", "gbm20171", "gbm20172", "gbm20173", "gbm20174", "gbm20175", "gbm20176",
            "gbm20180", "gbm20181", "gbm20182", "gbm20183", "gbm20184", "gbm20185", "gbm20186",
            "gbm20190", "gbm20191", "gbm20192", "gbm20193", "gbm20195", "gbm20196", "gbm20200",
            "gbm20201", "gbm20202", "gbm20203", "gbm20204", "gbm20205", "gbm20206", "gbm20210",
            "gbm20211", "gbm20212", "gbm20213", "gbm20214", "gbm20215", "gbm20216", "gbm20220",
            "gbm20221", "gbm20222", "gbm20223", "gbm20224", "gbm20225", "gbm20226", "gbm20230",
            "gbm20231", "gbm20232", "gbm20233", "gbm20234", "gbm20235", "gbm20240", "gbm20241",
            "gbm20242", "gbm20243", "gbm20244", "gbm20245", "gbm20246", "gbm20250", "gbm20251",
            "gbm20252", "gbm20253", "gbm20254", "gbm20255", "gbm20256", "gbm20260", "gbm20261",
            "gbm20262", "gbm20263", "gbm20264", "gbm20265", "gbm20266", "gbm20270", "gbm20271",
            "gbm20272", "gbm20273", "gbm20274", "gbm20275", "gbm20276", "gbm20280", "gbm20281",
            "gbm20282", "gbm20283", "gbm20284", "gbm20285", "gbm20286", "gbm20290", "gbm20291",
            "gbm20292", "gbm20293", "gbm20294", "gbm20295", "gbm20296", "gbm20300", "gbm20301",
            "gbm20302", "gbm20303", "gbm20304", "gbm20305", "gbm20306", "gbm20310", "gbm20311",
            "gbm20312", "gbm20313", "gbm20314", "gbm20315", "gbm20316", "gbm20320", "gbm20321",
            "gbm20322", "gbm20323", "gbm20324", "gbm20325", "gbm20326", "gbm20330", "gbm20331",
            "gbm20332", "gbm20333", "gbm20334", "gbm20335", "gbm20336", "gbm20340", "gbm20341",
            "gbm20342", "gbm20343", "gbm20344", "gbm20345", "gbm20346",
        ],
        _ => vec![],
    }
}

/// Energy term data from gqcd_chrono_mass detailed analysis
#[derive(Debug, Clone)]
pub struct MassEnergyTerms {
    pub term_time: f64,      // c² × Δrate (clock-derived energy)
    pub term_kinetic: f64,   // 0.5 × Δv² (orbit-derived kinetic)
    pub term_potential: f64, // Δ(1/r) (orbit-derived potential)
    pub g_value: f64,        // Derived G value
    pub r_a: f64,            // Radius at point A
    pub r_b: f64,            // Radius at point B
}

/// Load energy terms from gqcd_chrono_mass detailed analysis CSV file
pub fn load_mass_data(
    path: &std::path::Path,
    data: &mut Vec<MassEnergyTerms>,
) -> std::io::Result<()> {
    use std::fs::File;
    use std::io::{BufRead, BufReader};

    let file = File::open(path)?;
    let reader = BufReader::new(file);

    for (idx, line) in reader.lines().enumerate() {
        let line = line?;
        if idx == 0 {
            continue; // Skip header
        }

        // CSV: Timestamp,R_A,R_B,Delta_H_km,G,Error_Pct,Term_Time,Term_Kinetic,Term_Potential
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() >= 9
            && let (Ok(r_a), Ok(r_b), Ok(g), Ok(term_time), Ok(term_kinetic), Ok(term_potential)) = (
                parts[1].parse::<f64>(),
                parts[2].parse::<f64>(),
                parts[4].parse::<f64>(),
                parts[6].parse::<f64>(),
                parts[7].parse::<f64>(),
                parts[8].parse::<f64>(),
            )
        {
            data.push(MassEnergyTerms {
                term_time,
                term_kinetic,
                term_potential,
                g_value: g,
                r_a,
                r_b,
            });
        }
    }

    Ok(())
}

/// Load all mass experiment data from output directory
pub fn load_all_mass_data() -> std::io::Result<Vec<MassEnergyTerms>> {
    use std::fs;

    let mut all_data = Vec::new();
    let input_base_path = get_data_output_path("gqcd_chrono_mass");

    for year in YEARS.iter() {
        let year_path = input_base_path.join(year);
        if !year_path.exists() {
            continue;
        }

        for entry in fs::read_dir(&year_path)? {
            let entry = entry?;
            let file_path = entry.path();

            if let Some(name) = file_path.file_name().and_then(|n| n.to_str())
                && name.ends_with("_detailed_analysis.csv")
            {
                load_mass_data(&file_path, &mut all_data)?;
            }
        }
    }

    Ok(all_data)
}
