use deep_causality_num::RealField;

/// Raman interferometry shot from arbitrary orientations dataset
#[derive(Debug, Clone)]
#[allow(dead_code)] // Fields retained for future analysis extensions
pub struct RamanShot<R>
where
    R: RealField,
{
    iteration: u32,
    timestamp: String,
    phase_rad: R,     // Phase in radians
    ratio_mean: R,    // Population ratio (interference fringe)
    ratio_std_dev: R, // Uncertainty
}

impl<R> RamanShot<R>
where
    R: RealField,
{
    pub fn new(
        iteration: u32,
        timestamp: String,
        phase_rad: R,
        ratio_mean: R,
        ratio_std_dev: R,
    ) -> Self {
        Self {
            iteration,
            timestamp,
            phase_rad,
            ratio_mean,
            ratio_std_dev,
        }
    }
}

impl<R> RamanShot<R>
where
    R: RealField + Clone,
{
    pub fn iteration(&self) -> u32 {
        self.iteration
    }

    pub fn timestamp(&self) -> &str {
        &self.timestamp
    }

    pub fn phase_rad(&self) -> R {
        self.phase_rad
    }

    pub fn ratio_mean(&self) -> R {
        self.ratio_mean
    }

    pub fn ratio_std_dev(&self) -> R {
        self.ratio_std_dev
    }
}

/// Einstein Elevator experiment data (extracted from MATLAB)
#[derive(Debug, Clone)]
#[allow(dead_code)] // Fields retained for future analysis extensions
pub struct EinsteinElevatorData<R>
where
    R: RealField,
{
    freq_data: Vec<R>,       // Frequency sweep (Hz)
    ratio_data: Vec<R>,      // Population ratios
    pe_simulated: Vec<R>,    // Simulated transition probability
    residual: Vec<R>,        // Model - Experiment
    interrogation_time_s: R, // T (seconds)
    temperature_k: R,        // Atom temperature (Kelvin)
}

impl<R> EinsteinElevatorData<R>
where
    R: RealField,
{
    pub fn new(
        freq_data: Vec<R>,
        ratio_data: Vec<R>,
        pe_simulated: Vec<R>,
        residual: Vec<R>,
        interrogation_time_s: R,
        temperature_k: R,
    ) -> Self {
        Self {
            freq_data,
            ratio_data,
            pe_simulated,
            residual,
            interrogation_time_s,
            temperature_k,
        }
    }
}

impl<R> EinsteinElevatorData<R>
where
    R: RealField + Clone,
{
    pub fn freq_data(&self) -> &Vec<R> {
        &self.freq_data
    }

    pub fn ratio_data(&self) -> &Vec<R> {
        &self.ratio_data
    }

    pub fn pe_simulated(&self) -> &Vec<R> {
        &self.pe_simulated
    }

    pub fn residual(&self) -> &Vec<R> {
        &self.residual
    }

    pub fn interrogation_time_s(&self) -> R {
        self.interrogation_time_s
    }

    pub fn temperature_k(&self) -> R {
        self.temperature_k
    }
}

/// Analysis results for Action-Phase equivalence
#[derive(Debug, Clone)]
pub struct ActionPhaseResults<R>
where
    R: RealField,
{
    sample_size: usize,
    t_field_gradient: R,     // ∇τ = g/c² (s/m)
    mean_phase_rad: R,       // Mean measured phase
    action_measured_j_s: R,  // S = ℏ × Δφ
    action_predicted_j_s: R, // S = -mc² × Δτ
    action_ratio: R,         // S_measured / S_predicted
    phase_std_dev: R,
    energy_balance_ratio: R,
}

impl<R> ActionPhaseResults<R>
where
    R: RealField,
{
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        sample_size: usize,
        t_field_gradient: R,
        mean_phase_rad: R,
        action_measured_j_s: R,
        action_predicted_j_s: R,
        action_ratio: R,
        phase_std_dev: R,
        energy_balance_ratio: R,
    ) -> Self {
        Self {
            sample_size,
            t_field_gradient,
            mean_phase_rad,
            action_measured_j_s,
            action_predicted_j_s,
            action_ratio,
            phase_std_dev,
            energy_balance_ratio,
        }
    }
}

impl<R> ActionPhaseResults<R>
where
    R: RealField + Clone,
{
    pub fn sample_size(&self) -> usize {
        self.sample_size
    }

    pub fn t_field_gradient(&self) -> R {
        self.t_field_gradient
    }

    pub fn mean_phase_rad(&self) -> R {
        self.mean_phase_rad
    }

    pub fn action_measured_j_s(&self) -> R {
        self.action_measured_j_s
    }

    pub fn action_predicted_j_s(&self) -> R {
        self.action_predicted_j_s
    }

    pub fn action_ratio(&self) -> R {
        self.action_ratio
    }

    pub fn phase_std_dev(&self) -> R {
        self.phase_std_dev
    }

    pub fn energy_balance_ratio(&self) -> R {
        self.energy_balance_ratio
    }
}
