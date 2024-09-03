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
/// # Fields
/// 
/// * `k`: The proportionality constannt in the polytropic equation of state.
/// * `gamma`: The polytropic index, related to the compressibility of the stellar material.
/// * `g`: The gravitational constant.
/// 
/// # Physical Assumptions
/// 
/// The model assumes:
/// - Spherical symmetry of the star.
/// - The material follows a polytropic equation of state: P = k * rho^gamma.
/// 
/// # Parameter Choices
/// 
/// - `gamma` is often chosen as 5/3 for a non-relativistic degenerate gas, or 4/3 for an
/// ultra-relativistic gas.
/// - `g` is the standard gravitational constant: 6.6743e-8 dyn cm^2 g^-2.
#[derive(Debug, Clone)]
pub struct PhysicalConstants {
    pub k: f64,
    pub gamma: f64,
    pub g: f64,
}

/// Parameters that control the integration process.
/// 
/// This struct defines parameters such as the logarithmic step size (`log_dr`).
#[derive(Debug, Clone)]
pub struct IntegrationParams {
    pub log_dr: f64
}

impl Default for IntegrationParams {
    /// Provides default integration parameters.
    /// 
    /// Returns a default step size of `0.01` for `log_dr`.
    fn default() -> Self {
        Self {
            log_dr: 0.01
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
    params: IntegrationParams,
    surface_pressure: f64,
    surface_density: f64
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
    /// * - `surface_pressure` - The pressure at the surface of the star. Typically calculated from the polytropic EoS.
    /// * - `surface_density` - The density at the surface of the star. Used as a threshold to stop integration.
    /// 
    /// # Example
    /// 
    /// ```
    /// let constants = PhysicalConstants { k: 1.0e13, gamma: 5.0 / 3.0, g: 6.67430e-8 };
    /// let params = IntegrationParams::default();
    /// let model = StellarModel::new(constants, 1.0e-6, 1.0e6, 1.0e14, Some(params))l
    /// ```
    pub fn new(constants: PhysicalConstants, r_start: f64, r_end: f64, rho_c: f64, params: Option<IntegrationParams>, surface_pressure: f64, surface_density: f64) -> Self {
        StellarModel {
            structure: StellarStructure::new(constants),
            r_start,
            r_end,
            rho_c,
            params: params.unwrap_or_default(),
            surface_pressure,
            surface_density
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
        let _log_r_end: f64 = self.r_end.log10();
        let log_dr: f64 = self.params.log_dr;
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

        let mut last_valid_state = state.clone();
        let mut last_valid_log_r: f64 = log_r;

        // while log_r < log_r_end {
        while !self.is_dm_dr_negligible(log_r, &state) {
            self.write_state_to_file(&mut file, log_r, &state)?;

            if self.is_termination_condition_met( log_r, &state) {
                state = last_valid_state;
                log_r = last_valid_log_r;
                break;
            }
            
            last_valid_state = state.clone();
            last_valid_log_r = log_r;
            
            let next_state: State<f64> = solver.step(log_r, &state, log_dr);

            // let log_rho: f64 = (state.values[1] - self.structure.constants.k.log10()) / self.structure.constants.gamma;
            // let next_log_rho: f64 = (next_state.values[1] - self.structure.constants.k.log10()) / self.structure.constants.gamma;
            // let log_rho_diff: f64 = (next_log_rho - log_rho).abs();

            // if log_rho_diff > 1e-2 {
            //     log_dr /= 2.0;
            // } else if log_rho_diff < 1e-4 {
            //     log_dr *= 2.0;
            // }
            
            state = next_state;
            log_r += log_dr;
        }

        if self.is_dm_dr_negligible(log_r, &state) {
            println!("dm/dr is negligible");
        }

        let r_1 = base.powf(log_r);
        let m_1 = base.powf(state.values[0]);
        let p_1 = base.powf(state.values[1]);

        let r_2: f64 = self.find_surface_radius(r_1, p_1, m_1)?;

        Ok((r_2, m_1, self.surface_pressure))
    }

    /// Checks if the mass continuity equation (dm/dr) is negligible.
    /// 
    /// This function determines whether the change in mass with respect to radius
    /// has become negligible, indicating that we're approaching the star's surface.
    /// It also checks if the density has fallen below the defined surface density.
    /// 
    /// # Arguments
    /// 
    /// * `log_r` - The logarithm (base 10) of the radial coordinate.
    /// * `state` - A `State<f64>` containing the logarithms of the mass and pressure.
    /// 
    /// # Returns
    /// 
    /// Returns `true` if either:
    /// - The absolute value of dm/dr is less than 1e-6, or
    /// - The current density is lower than or equal to the defined surface density.
    fn is_dm_dr_negligible(&self, log_r: f64, state: &State<f64>) -> bool{
        let derivatives: State<f64> = self.structure.derivatives(log_r, state);
        let log_rho: f64 = (state.values[1] - self.structure.constants.k.log10()) / self.structure.constants.gamma;
        let rho: f64 = 10f64.powf(log_rho);

        derivatives.values[0].abs() < 1e-6 || rho <= self.surface_density
    }

    /// Calculates the total radius of the star using the bisection method.
    /// 
    /// This function is called after dm/dr becomes negligible, to precisely locate
    /// the star's surface where the density equals the defined surface density.
    /// 
    /// # Arguments
    /// 
    /// * `r_1` - Radius where dm/dr became negligible.
    /// * `p_1` - Pressure at `r_1`.
    /// * `m_1` - Mass at `r_1`, considered the total mass of the star.
    /// 
    /// # Returns
    /// 
    /// Returns the calculated surface radius of the star.
    /// Returns `StellarError` if an error occurs during pressure integration.
    /// 
    /// # Algorithm
    /// 
    /// Uses the bisection method to find the radius where the density equals
    /// the defined surface density. The search starts from `r_1` and extends
    /// outward by 20% of `r_1`.
    /// 
    /// # Notes
    /// 
    /// The tolerance for both the bisection method and the surface density
    /// comparison is set to 1e-6 * r_1.
    fn find_surface_radius(&self, r_1: f64, p_1: f64, m_1: f64) -> Result<f64, StellarError> {
        let mut a: f64 = r_1;
        let mut b: f64 = r_1 * 1.2;
        let tolerance: f64 = 1e-6 * r_1;

        while b - a > tolerance {
            let c: f64 = (a + b) / 2.0;
            let p_c: f64 = self.integrate_pressure(r_1, c, p_1, m_1)?;

            let rho_c: f64 = (p_c / self.structure.constants.k).powf(1.0 / self.structure.constants.gamma);
            if (rho_c - self.surface_density).abs() < tolerance {
                return Ok(c);
            }

            let p_a: f64 = self.integrate_pressure(r_1, a, p_1, m_1)?;
            let rho_a: f64 = (p_a / self.structure.constants.k).powf(1.0 / self.structure.constants.gamma);
            if (rho_c - self.surface_density) * (rho_a - self.surface_density) > 0.0 {
                a = c;
            } else {
                b = c;
            }
        }

        Ok((a+b) / 2.0)
    }

    /// Integrates the pressure equation from a starting radius to an ending radius.
    /// 
    /// This function is used when dm/dr has become negligible, allowing us to
    /// integrate only the pressure equation while assuming constant mass.
    /// 
    /// # Arguments
    /// 
    /// * `r_start` - Starting radius for the integration.
    /// * `r_end` - Ending radius for the integration.
    /// * `p_start` - Pressure at `r_start`.
    /// * `m` - Total mass of the star (constant in this region).
    /// 
    /// # Returns
    /// 
    /// Returns the calculated pressure at `r_end`.
    /// Returns `StellarError` if an error occurs during integration.
    /// 
    /// # Method
    /// 
    /// Uses a simple Euler method with 1000 steps to integrate dP/dr from
    /// `r_start` to `r_end`. The step size is calculated as (r_end - r_start) / 1000.
    /// 
    /// # Notes
    /// 
    /// This method assumes that the mass is constant in the region of integration,
    /// which is valid when dm/dr is negligible near the star's surface.
    fn integrate_pressure(&self, r_start: f64, r_end: f64, p_start: f64, m: f64) -> Result<f64, StellarError> {
        let n_steps: i32 = 1000;
        let dr: f64 = (r_end - r_start) / n_steps as f64;

        let mut p: f64 = p_start;
        let mut r: f64 = r_start;

        for _ in 0..n_steps {
            let rho: f64 = (p / self.structure.constants.k).powf(1.0 / self.structure.constants.gamma);
            let dp_dr: f64 = -self.structure.constants.g * m * rho / (r * r);
            p += dp_dr * dr;
            r += dr;
        }

        Ok(p)

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
    fn is_termination_condition_met(&self, log_r: f64, state: &State<f64>) -> bool {
        let base: f64 = 10.0;
        if state.values[0].is_nan() || state.values[1].is_nan() {
            println!("Instability Detected at r = {:.2e}, m = {:.2e}, P = {:.2e}", base.powf(log_r), base.powf(state.values[0]), base.powf(state.values[1]));
            return true;
        }

        false
    } 

}