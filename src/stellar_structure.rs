use crate::rk4::{Derivatives, RK4Solver, State};
use std::fs::File;
use std::io::Write;
use std::io::Error as IoError;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum StellarError {
    #[error("IO error: {0}")]
    Io(#[from] IoError),
    #[error("Math error: {0}")]
    Math(String)
}

#[derive(Debug, Clone)]
pub struct PhysicalConstants {
    pub k: f64,
    pub gamma: f64,
    pub g: f64,
}

#[derive(Debug, Clone)]
pub struct IntegrationParams {
    pub log_dr: f64,
    pub surface_pressure_threshold: f64
}

impl Default for IntegrationParams {
    fn default() -> Self {
        Self {
            log_dr: 0.01,
            surface_pressure_threshold: 1.1
        }
    }
}

pub struct StellarStructure {
    constants: PhysicalConstants
}

impl StellarStructure {
    pub fn new(constants: PhysicalConstants) -> Self {
        Self { constants }
    }
}

impl Derivatives<f64> for StellarStructure {
    fn derivatives(&self, log_r: f64, state: &State<f64>) -> State<f64> {
        let epsilon: f64 = 1e-10; // Regularization Term

        let log_m: f64 = state.values[0];
        let log_p: f64 = state.values[1];

        let base: f64 = 10.0;
        let log_rho: f64 = (log_p - self.constants.k.log10()) / self.constants.gamma;

        let dlogm_dlogr: f64 = 4.0 * std::f64::consts::PI * base.powf(log_r + epsilon).powi(3) * base.powf(log_rho + epsilon) / base.powf(log_m + epsilon);
        let dlogp_dlogr: f64 = - self.constants.g * base.powf(log_m + epsilon) * base.powf(log_rho + epsilon) / (base.powf(log_p + epsilon) * base.powf(log_r + epsilon));
    

        State {
            values: vec![dlogm_dlogr, dlogp_dlogr]
        }
    }
}

pub struct StellarModel {
    pub structure: StellarStructure,
    pub r_start: f64,
    pub r_end: f64,
    pub rho_c: f64,
    params: IntegrationParams
}

impl StellarModel {
    pub fn new(constants: PhysicalConstants, r_start: f64, r_end: f64, rho_c: f64, params: Option<IntegrationParams>) -> Self {
        StellarModel {
            structure: StellarStructure::new(constants),
            r_start,
            r_end,
            rho_c,
            params: params.unwrap_or_default()
        }
    }

    pub fn run(&self, output_file: &str) -> Result<(f64, f64), StellarError> {
        let base: f64 = 10.0;
        let log_r_start: f64 = self.r_start.log10();
        let log_r_end: f64 = self.r_end.log10();
        let _dr: f64 = base.powf(self.params.log_dr);

        let log_rho_c: f64 = self.rho_c.log10();

        let fraction: f64 = 4.0 / 3.0;
        let log_m0: f64 = fraction.log10() + std::f64::consts::PI.log10() + 3.0 * self.r_start.log10() + self.rho_c.log10();
        let log_p0: f64 = self.structure.constants.k.log10() + (self.structure.constants.gamma * log_rho_c);

        let mut state: State<f64> = State { values: vec![log_m0, log_p0] };

        let solver: RK4Solver<'_, f64, StellarStructure> = RK4Solver::new(&self.structure, self.params.log_dr);

        let mut file = File::create(output_file)?;
        writeln!(file, "r,m,P")?;
    
        let mut log_r: f64 = log_r_start;
    
        let log_pressure_threshold: f64 = self.params.surface_pressure_threshold.log10();

        while log_r < log_r_end {
            self.write_state_to_file(&mut file, log_r, &state)?;
    
            if self.is_terminationn_condition_met(log_r, &state, log_pressure_threshold) {
                break;
            }
    
            state = solver.step(log_r, &state);
            log_r += self.params.log_dr;
        }
        Ok((base.powf(state.values[0]), base.powf(state.values[1])))
    }

    fn write_state_to_file(&self, file: &mut File, log_r: f64, state: &State<f64>) -> Result<(), StellarError> {
        let base: f64 = 10.0;
        writeln!(file, "{},{},{}", base.powf(log_r), base.powf(state.values[0]), base.powf(state.values[1]))?;
        Ok(())
    }

    fn is_terminationn_condition_met(&self, log_r: f64, state: &State<f64>, log_pressure_threshold: f64) -> bool {
        let base: f64 = 10.0;
        if state.values[1].is_nan() {
            println!("Surface or Instability Detected at r = {:.2e}, m = {:.2e}, P = {:.2e}", base.powf(log_r), base.powf(state.values[0]), base.powf(state.values[1]));
            return true;
        }
        if state.values[1] < log_pressure_threshold {
            println!("Surface Detected At r = {:.2e}, m = {:.2e}, P = {:.2e}", base.powf(log_r), base.powf(state.values[0]), base.powf(state.values[1]));
            return true;
        }
        false
    } 

}