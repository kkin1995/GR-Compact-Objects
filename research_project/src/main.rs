use gr_compact_objects::rk4::{Derivatives, RK4Solver, State};
use std::io::{Result, Write};
use std::fs::File;
use std::os::macos::raw::stat;

struct StellarStructure {
    k: f64,
    gamma: f64,
    g: f64,
}

impl Derivatives<f64> for StellarStructure {
    fn derivatives(&self, log_r: f64, state: &State<f64>) -> State<f64> {
        let epsilon: f64 = 1e-10; // Regularization Term

        // let m: f64 = state.values[0];
        // let p: f64 = state.values[1];
        let log_m: f64 = state.values[0];
        let log_p: f64 = state.values[1];

        let base: f64 = 10.0;
        // let rho: f64 = base.powf((p.log10() - self.k.log10()) / self.gamma);
        let log_rho: f64 = (log_p - self.k.log10()) / self.gamma;
        // let rho: f64 = (p / self.k).powf(1.0 / self.gamma);

        let dlogm_dlogr: f64 = 4.0 * std::f64::consts::PI * base.powf(log_r).powi(3) * base.powf(log_rho) / base.powf(log_m);
        let dlogp_dlogr: f64 = - self.g * base.powf(log_m) * base.powf(log_rho) / (base.powf(log_p) * base.powf(log_r));
        // let dm_dr: f64 = 4.0 * std::f64::consts::PI * r * r * rho;
        // let dp_dr: f64 = - self.g * m * rho / (r * r);
        
        // println!("Debug: r: {:.2e} log_r: {:.2e} m: {:.2e} log_m: {:.2e} P: {:2e} log_p: {:.2e} rho: {:2e} log_rho: {:.2e} dlogm_dlogr: {:2e} dlogp_dlogr: {:2e}",
        //                 base.powf(log_r), log_r, base.powf(log_m), log_m, base.powf(log_p), log_p, base.powf(log_rho), log_rho, dlogm_dlogr, dlogp_dlogr);

        State {
            values: vec![dlogm_dlogr, dlogp_dlogr]
        }
    }
}

fn main() -> Result<()> {
    let stellar_structure = StellarStructure {
        k: 5.3802e9, // CGS
        gamma: 5.0 / 3.0, // Non-Relativistic Fermion Gas
        g: 6.67430e-8, // Gravitational Constant // CGS
    };

    let r_start: f64 = 100.0;
    let log_r_start: f64 = r_start.log10();
    let r_end: f64 = 1e6;
    let log_r_end: f64 = r_end.log10();
    let dr: f64 = 1.0;

    let rho_c: f64 = 1e12;
    let log_rho_c: f64 = rho_c.log10();
    let m0: f64 = 4.0 * std::f64::consts::PI * r_start.powf(3.0) * rho_c / 3.0;
    let fraction: f64 = 4.0 / 3.0;
    let log_m0: f64 = fraction.log10() + std::f64::consts::PI.log10() + 3.0 * r_start.log10() + rho_c.log10();
    let p0: f64 = stellar_structure.k * rho_c.powf(stellar_structure.gamma);
    let log_p0: f64 = stellar_structure.k.log10() + (stellar_structure.gamma * log_rho_c);

    let mut state: State<f64> = State { values: vec![log_m0, log_p0] };

    let solver: RK4Solver<f64, StellarStructure> = RK4Solver::new(stellar_structure, dr);

    // let mut file = File::create("stellar_structure.csv")?;

    // writeln!(file, "r,m,P")?;

    let mut r: f64 = r_start;
    let mut log_r: f64 = log_r_start;

    let pressure_threshold: f64 = 1e-10;
    let log_pressure_threshold: f64 = pressure_threshold.log10();

    let base: f64 = 10.0;
    while r < r_end {
        // writeln!(file, "{},{},{}", r, base.powf(state.values[0]), base.powf(state.values[1]))?;
        println!("r: {:.2}, m: {:.6e}, P: {:.6e}", r, base.powf(state.values[0]), base.powf(state.values[1]));
        state = solver.step(log_r, &state);

        r += dr;
        log_r = r.log10();

        if state.values[1] < log_pressure_threshold {
            println!("Surface Detected At r = {:.2e}, m = {:.2e}, P = {:.2e}", r, base.powf(state.values[0]), base.powf(state.values[1]));
            break;
        }
    }

    println!("Final Mass: {}", state.values[0]);
    println!("Final Pressure: {}", state.values[1]);

    Ok(())
}