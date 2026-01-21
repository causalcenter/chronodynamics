use deep_causality_num::RealField;

#[derive(Debug, Clone)]
pub struct PulsarObservation<R>
where
    R: RealField,
{
    pub name: String,
    pub observations: Vec<PulsarDataPoint<R>>,
}

#[derive(Debug, Clone)]
pub struct PulsarDataPoint<R>
where
    R: RealField,
{
    pub timestamp_unix: i64,
    pub mjd: R,
    pub error_us: R,
    pub freq_mhz: R,
    pub residual_s: R, // Currently just 0.0, placeholder for model subtraction
}
