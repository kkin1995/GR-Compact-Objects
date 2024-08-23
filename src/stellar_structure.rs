//! Stellar Structure and Modelling Library
//! 
//! This module provides the necessary components to model the structure of a 
//! compact body, such as a Neutron Star, using the Runge-Kutta 4th Order (RK4)
//! integration method. The module includes physical constants, integration parameters,
//! and functionality to calculate the mass and pressure profiles of a star.
use crate::rk4::{Derivatives, RK4Solver, State};
use std::fs::File;
use std::io::Write;
use std::f64::consts::PI;
use std::io::Error as IoError;
use thiserror::Error;

/// Custom Error Type
/// 
/// This enum represents the possible errors that can occur during the execution
/// of the stellar structure integration, including I/O errors and custom mathematical
/// errors.
#[derive(Error, Debug)]
pub enum StellarError {
    /// Error related to input/output operations.
    #[error("IO error: {0}")]
    Io(#[from] IoError),

    /// Custom error for mathematical issues encountered during integration.
    #[error("Math error: {0}")]
    Math(String)
}

/// Physical constants required for the stellar structure equations.
/// 
/// This struct holds the essential physical constants such as the proportionality
/// constant (`k`), the polytropic index (`gamma`), and the gravitational constant (`g`).
#[derive(Debug, Clone)]
pub struct PhysicalConstants {
    pub k: f64,
    pub gamma: f64,
    pub g: f64,
}

/// Parameters that control the integration process.
/// 
/// This struct defines parameters such as the logarithmic step size (`log_dr`)
/// and the surface pressure threshold (`surface_pressure_threshold`).
#[derive(Debug, Clone)]
pub struct IntegrationParams {
    pub log_dr: f64,
    pub surface_pressure_threshold: f64
}

impl Default for IntegrationParams {
    /// Provides default integration parameters.
    /// 
    /// Returns a default step size of `0.01` for `log_dr` and a surface pressure
    /// threshold of `1.1`.
    fn default() -> Self {
        Self {
            log_dr: 0.01,
            surface_pressure_threshold: 1.1
        }
    }
}

/// Represents the stellar structure and provides methods to compute derivatives.
/// 
/// This struct contains the physical constants and implements the `Derivatives<f64>` trait
/// to calculate the derivatives of logarithmic mass and pressure with respect to the logarithmic radius.
pub struct StellarStructure {
    constants: PhysicalConstants
}

impl StellarStructure {
    /// Creates a new `StellarStructure` instance.
    /// 
    /// # Arguments
    /// 
    /// * `constants` - Physical constants required for the stellar structure equations.
    /// 
    /// # Example
    /// 
    /// ```
    /// let constants = PhysicalConstants{ k: 1.0e13, gamma: 5.0 / 3.0, g: 6.67430e-8 };
    /// let structure = StellarStructure::new(constants);
    /// ```
    pub fn new(constants: PhysicalConstants) -> Self {
        Self { constants }
    }
}

impl Derivatives<f64> for StellarStructure {
    /// Computes the derivatives of logarithmic mass and pressure with respect to the logarithmic radius.
    /// 
    /// This function calculates the logarithmic derivaives needed for the RK4 integration method.
    /// 
    /// # Arguments
    /// 
    /// * `log_r` - The logarithm of the current radius.
    /// * `state` - The current state, containing the logarithmic mass and pressure.
    /// 
    /// # Returns
    /// 
    /// A `State<f64>` containing the derivatives of the logarithmic mass and pressure.
    fn derivatives(&self, log_r: f64, state: &State<f64>) -> State<f64> {
        let epsilon: f64 = 1e-10; // Regularization Term
        let relative_epsilon = |x: f64| epsilon * (1.0 + x.abs());
        let log_m: f64 = state.values[0];
        let log_p: f64 = state.values[1];

        let base: f64 = 10.0;
        let log_rho: f64 = (log_p - self.constants.k.log10()) / self.constants.gamma;
        
        if log_r.is_nan() || log_m.is_nan() || log_p.is_nan() || log_rho.is_nan() ||
        log_r.is_infinite() || log_m.is_infinite() || log_p.is_infinite() || log_rho.is_infinite() {
            return State { values: vec![f64::NAN, f64::NAN] };
        }

        let safe_log_m = log_m + relative_epsilon(log_m);
        let safe_log_p = log_p + relative_epsilon(log_p);
        let safe_log_r = log_r + relative_epsilon(log_r);
        let safe_log_rho = log_rho + relative_epsilon(log_rho);

        if safe_log_r.is_sign_negative() || safe_log_m.is_sign_negative() || safe_log_p.is_sign_negative() || safe_log_rho.is_sign_negative() {
            return State { values: vec![f64::NAN, f64::NAN] };
        }

        // let dlogm_dlogr: f64 = 4.0 * std::f64::consts::PI * (safe_log_r * base.ln()).exp().powi(3) * (safe_log_rho * base.ln()).exp() / (safe_log_m * base.ln()).exp();
        // let dlogp_dlogr: f64 = - self.constants.g * (safe_log_m * base.ln()).exp() * (safe_log_rho * base.ln()).exp() / ((safe_log_p * base.ln()).exp() * (safe_log_r * base.ln()).exp());

        let dlogm_dlogr: f64 = (4.0 * PI * base.powf(safe_log_r).powf(3.0) * base.powf(safe_log_rho)) / base.powf(safe_log_m);
        let dlogp_dlogr: f64 = - (self.constants.g * base.powf(safe_log_m) * base.powf(safe_log_rho)) / (base.powf(safe_log_p) * base.powf(safe_log_r));
        
        // if (safe_log_r * 100.0).round() / 100.0 != (safe_log_r * 100.0).floor() / 100.0 {
        //     println!("Debug:\t\tr(cm)\t m(g)\tP(dyne/cm^2)\trho(g/cm^3)\tdm/dr\tdP/dr");
        //     println!("Debug\t\t{:.2e}\t{:.2e}\t{:.2e}\t{:.2e}\t{:.2e}\t{:.2e}",
        //         base.powf(safe_log_r),
        //         base.powf(safe_log_m),
        //         base.powf(safe_log_p),
        //         base.powf(safe_log_rho),
        //         dlogm_dlogr * base.powf(safe_log_m) / base.powf(safe_log_r), // dm/dr
        //         dlogp_dlogr * base.powf(safe_log_p) / base.powf(safe_log_r)  // dP/dr
        // );
        //     println!("Log Debug:\t\tlog(m)\tlog(p)\tlog(r)\tlog(rho)\tdlogm_dlogr\tdlogp_dlogr");
        //     println!("Log Debug:\t\t{:.2e}\t{:.2e}\t{:.2e}\t{:.2e}\t{:.2e}\t{:.2e}", safe_log_m, safe_log_p, safe_log_r, safe_log_rho, dlogm_dlogr, dlogp_dlogr);
        //     println!("----------------------------------------");
        // }

        if !dlogm_dlogr.is_finite() || !dlogp_dlogr.is_finite() ||
        dlogm_dlogr.abs() > 1e30 || dlogm_dlogr > 1e30 {
            println!("Warning: Extreme Derivate Values Detected At log_r = {:2e}", log_r);
            return State { values: vec![f64::NAN, f64::NAN] };
        }

        State {
            values: vec![dlogm_dlogr, dlogp_dlogr]
        }
    }
}

/// Encapsulates the stellar model, including the integration and output functionalities.
/// 
/// The `StellarModel` struct allows you to initialize a stellar model with specific parameters,
/// run the integration process, and output the results to a file.
pub struct StellarModel {
    pub structure: StellarStructure,
    pub r_start: f64,
    pub r_end: f64,
    pub rho_c: f64,
    params: IntegrationParams
}

impl StellarModel {
    /// Creates a new `StellarModel` instance.
    /// 
    /// # Arguments
    /// 
    /// * - `constants` -  The physical constants used in the model.
    /// * - `r_start` - The starting radius for the integration.
    /// * - `r_end` - The ending radius for the integration.
    /// * - `rho_c` - The central density of the compact body.
    /// * - `params` - Optional integration parameters. Defaults to `IntegrationParameters::default()` if `None`.
    /// 
    /// # Example
    /// 
    /// ```
    /// let constants = PhysicalConstants { k: 1.0e13, gamma: 5.0 / 3.0, g: 6.67430e-8 };
    /// let params = IntegrationParams::default();
    /// let model = StellarModel::new(constants, 1.0e-6, 1.0e6, 1.0e14, Some(params))l
    /// ```
    pub fn new(constants: PhysicalConstants, r_start: f64, r_end: f64, rho_c: f64, params: Option<IntegrationParams>) -> Self {
        StellarModel {
            structure: StellarStructure::new(constants),
            r_start,
            r_end,
            rho_c,
            params: params.unwrap_or_default()
        }
    }

    /// Runs the stellar structure integration and writes the results to a file.
    /// 
    /// This function integrates the stellar structure equatios from `r_start` to `r_end` and
    /// writes the radius, mass, and pressure profiles to the specified output file.
    /// 
    /// # Arguments
    /// 
    /// * `output_file` - The path to the file where the results will be written.
    /// 
    /// # Returns
    /// 
    /// A tuple containing the final radius, mass, and pressure at the star's surface.
    /// 
    /// # Errors
    /// 
    /// Returns a `StellarError` if an I/O error occurs or if there is a mathematical error during integration.
    pub fn run(&self, output_file: &str) -> Result<(f64, f64, f64), StellarError> {
        let base: f64 = 10.0;
        let log_r_start: f64 = self.r_start.log10();
        let log_r_end: f64 = self.r_end.log10();
        let mut log_dr: f64 = self.params.log_dr;
        let _dr: f64 = base.powf(log_dr);

        let log_rho_c: f64 = self.rho_c.log10();

        let fraction: f64 = 4.0 / 3.0;
        let log_m0: f64 = fraction.log10() + std::f64::consts::PI.log10() + 3.0 * self.r_start.log10() + self.rho_c.log10();
        let log_p0: f64 = self.structure.constants.k.log10() + (self.structure.constants.gamma * log_rho_c);

        let mut state: State<f64> = State { values: vec![log_m0, log_p0] };

        let solver: RK4Solver<'_, f64, StellarStructure> = RK4Solver::new(&self.structure);

        let mut file = File::create(output_file)?;
        writeln!(file, "r,m,P")?;
    
        let mut log_r: f64 = log_r_start;
    
        let _log_pressure_threshold: f64 = self.params.surface_pressure_threshold.log10();

        let mut last_valid_state = state.clone();
        let mut last_valid_log_r: f64 = log_r;

        while log_r < log_r_end {
            self.write_state_to_file(&mut file, log_r, &state)?;
    
            if self.is_termination_condition_met(log_r, &state, _log_pressure_threshold) {
                state = last_valid_state;
                log_r = last_valid_log_r;
                break;
            }
            
            last_valid_state = state.clone();
            last_valid_log_r = log_r;
            
            let next_state: State<f64> = solver.step(log_r, &state, self.params.log_dr);

            let log_rho: f64 = (state.values[1] - self.structure.constants.k.log10()) / self.structure.constants.gamma;
            let next_log_rho: f64 = (next_state.values[1] - self.structure.constants.k.log10()) / self.structure.constants.gamma;
            let log_rho_diff: f64 = (next_log_rho - log_rho).abs();

            if log_rho_diff > 1e-2 {
                log_dr /= 2.0;
            } else if log_rho_diff < 1e-4 {
                log_dr *= 2.0;
            }
            
            state = next_state;
            log_r += self.params.log_dr;
        }
        Ok((base.powf(log_r), base.powf(state.values[0]), base.powf(state.values[1])))
    }

    /// Writes the current state (radius, mass, pressure) to the output file.
    /// 
    /// This helper function formats and writes the current radius, mass, and pressure
    /// to the specified file.
    /// 
    /// # Arguments
    /// 
    /// * `file` - A mutable reference to the file where the state will be written.
    /// * `log_r` - The logarithm of the current radius.
    /// * `state` - The current state, containing the logarithmic mass and pressure.
    /// 
    /// # Returns
    /// 
    /// Returns `Ok(())` on success or a `StellarError` if an I/O error occurs. 
    fn write_state_to_file(&self, file: &mut File, log_r: f64, state: &State<f64>) -> Result<(), StellarError> {
        let base: f64 = 10.0;
        writeln!(file, "{},{},{}", base.powf(log_r), base.powf(state.values[0]), base.powf(state.values[1]))?;
        Ok(())
    }

    /// Checks whether the integration should terminate.
    /// 
    /// This function checks whether the integration process should terminate, either due to detecting
    /// the star's surface or encountering instability in the solution.
    /// 
    /// # Arguments
    /// 
    /// * `log_r` - The logarithm of the current radius.
    /// * `state` - THe current state, containing the logarithmic mass and pressure.
    /// * `_log_pressure_threshold` - Optional. The log of the surface pressure threshold.
    /// 
    /// # Returns
    /// 
    /// Returns `true` if a `NaN` value or the surface (according to `_log_surface_threshold`) is detected.
    /// Returns `false` otherwise.
    fn is_termination_condition_met(&self, log_r: f64, state: &State<f64>, _log_pressure_threshold: f64) -> bool {
        let base: f64 = 10.0;
        if state.values[0].is_nan() || state.values[1].is_nan() {
            println!("Instability Detected at r = {:.2e}, m = {:.2e}, P = {:.2e}", base.powf(log_r), base.powf(state.values[0]), base.powf(state.values[1]));
            return true;
        }
        if state.values[1] < _log_pressure_threshold {
            println!("Surface Detected At r = {:.2e}, m = {:.2e}, P = {:.2e}", base.powf(log_r), base.powf(state.values[0]), base.powf(state.values[1]));
            return true;
        }
        false
    } 

}