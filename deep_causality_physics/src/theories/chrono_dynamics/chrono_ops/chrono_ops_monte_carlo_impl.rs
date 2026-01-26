use crate::{ChronoGauge, ChronoGaugeMonteCarloOps};
use deep_causality_num::RealField;
use deep_causality_topology::{RandomField, TopologyError};
use std::fmt::Debug;

impl<R> ChronoGaugeMonteCarloOps<R> for ChronoGauge<R>
where
    R: RealField
        + Clone
        + From<f64>
        + Into<f64>
        + Default
        + Debug
        + deep_causality_num::FromPrimitive
        + deep_causality_num::ToPrimitive
        + RandomField,
{
    // =========================================================================
    // Monte Carlo Methods
    // =========================================================================

    fn thermalize<Rng: deep_causality_rand::Rng>(
        &mut self,
        n_sweeps: usize,
        epsilon: R,
        rng: &mut Rng,
    ) -> Result<f64, TopologyError> {
        let mut total_acceptance = 0.0;

        for _ in 0..n_sweeps {
            let acceptance = self.try_metropolis_sweep(epsilon, rng)?;
            total_acceptance += acceptance;
        }

        if n_sweeps > 0 {
            Ok(total_acceptance / n_sweeps as f64)
        } else {
            Ok(0.0)
        }
    }

    fn measure_with_error<F, Rng>(
        &mut self,
        observable: F,
        n_measurements: usize,
        skip: usize,
        epsilon: R,
        rng: &mut Rng,
    ) -> Result<(R, R), TopologyError>
    where
        F: Fn(&Self) -> Result<R, TopologyError>,
        Rng: deep_causality_rand::Rng,
    {
        let mut sum = R::zero();
        let mut sum_sq = R::zero();
        let mut count = 0;

        // Thermalize first (optional, user should call thermalize before)
        // But we do need decorrelation steps (skip) between measurements

        for _ in 0..n_measurements {
            // Decorrelation sweeps
            for _ in 0..skip {
                self.try_metropolis_sweep(epsilon, rng)?;
            }

            // Measurement
            let val = observable(self)?;
            sum += val;
            sum_sq += val * val;
            count += 1;
        }

        if count == 0 {
            return Ok((R::zero(), R::zero()));
        }

        let n = <R as From<f64>>::from(count as f64);
        let mean = sum / n;

        // Variance = <x^2> - <x>^2
        let mean_sq = sum_sq / n;
        let variance = mean_sq - mean * mean;

        // Standard error = sqrt(Variance / N)
        // Note: This assumes independent samples. If skip is too small,
        // autocorrelation increases the error (needs jackknife/bootstrap).
        let error = RealField::sqrt(variance / <R as From<f64>>::from((count as f64).abs()));

        Ok((mean, error))
    }

    fn auto_tune_metropolis<Rng: deep_causality_rand::Rng>(
        &mut self,
        initial_epsilon: R,
        target_acceptance: f64,
        max_steps: usize,
        rng: &mut Rng,
    ) -> Result<R, TopologyError> {
        let mut epsilon = initial_epsilon;
        let checks = 50; // Number of sweeps per check

        for _ in 0..max_steps {
            let mut total_acc = 0.0;
            for _ in 0..checks {
                total_acc += self.try_metropolis_sweep(epsilon, rng)?;
            }
            let rate = total_acc / checks as f64;

            // Simple adaptive logic
            // Ideal rate is usually around 0.5 for metropolis
            if (rate - target_acceptance).abs() < 0.05 {
                return Ok(epsilon);
            }

            if rate < target_acceptance {
                // Acceptance too low -> steps too big -> decrease epsilon
                let factor = <R as From<f64>>::from(0.9);
                epsilon *= factor;
            } else {
                // Acceptance too high -> steps too small -> increase epsilon
                let factor = <R as From<f64>>::from(1.1);
                epsilon *= factor;
            }
        }

        Ok(epsilon)
    }
}
