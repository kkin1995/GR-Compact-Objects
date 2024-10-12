use gr_compact_objects::stellar_structure::{EquationOfState, IntegrationParams, ModelType, PhysicalConstants, StellarModel};
use std::fs::File;
use std::io::Write;
use std::io::Error as IoError;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum SimulationError {
    #[error("IO error: {0}")]
    Io(#[from] IoError),
    #[error("Stellar model error: {0}")]
    StellarModelError(#[from] gr_compact_objects::stellar_structure::StellarError),
}
fn create_neutron_star_model(k: f64, gamma: f64, model_type: ModelType, r_start: f64, r_end: f64, rho_c: f64, log_dr: f64, surface_density: f64) -> StellarModel {
    let g: f64 = 6.67430e-8; // gravitational constant in CGS
    let c: f64 = 2.99792458e10; // speed of light in CGS
    let surface_pressure: f64 = k * surface_density.powf(gamma);

    let eos = EquationOfState {
        k: k,
        gamma: gamma
    };

    let constants = PhysicalConstants {
        g,
        c
    };

    let integration_params = IntegrationParams {
        log_dr,
        model_type
    };

    StellarModel::new(
        eos,
        constants,
        1.0,
        1e12,
        rho_c,
        Some(integration_params),
        surface_pressure,
        surface_density
    )
}

fn run_compact_star_simulation(k: f64, gamma: f64, model_type: &ModelType, r_start: f64, r_end: f64, rho_c: f64, log_dr: f64, surface_density: f64, output_file: &str) -> Result<(f64, f64, f64), gr_compact_objects::stellar_structure::StellarError> {
    let model: StellarModel = create_neutron_star_model(k, gamma, *model_type, r_start, r_end, rho_c, log_dr, surface_density);
    model.run(output_file)
}

fn get_gas_parameters(gas_type: &str, mu_e: f64) -> (f64, f64, f64, f64, ModelType) {
    match gas_type {
        "non_relativistic_electron_gas" => {
            // Non - Relativistic Electrons
            (1.0, 5.0, 1.0036e13 / mu_e.powf(5.0 / 3.0), 5.0 / 3.0, ModelType::Newton)
        },

        "relativistic_electron_gas" => {
            (7.0, 12.0, 1.2435e15 / mu_e.powf(4.0 / 3.0), 4.0 / 3.0, ModelType::TOV)
        },

        "non_relativistic_neutron_gas" => {
            (9.0, 14.0, 5.3802e9, 5.0 / 3.0, ModelType::Newton)
        },

        "relativistic_neutron_gas" => {
            (16.0, 20.0, 1.2293e15, 4.0 / 3.0, ModelType::TOV)
        },

        _ => panic!("Unknown Gas Type: {}", gas_type)
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n: i32 = 6; // Replace this with the desired number of values
    let gas_type: &str = "relativistic_neutron_gas";
    let mu_e: f64 = 2.0;
    let r_start: f64 = 0.001;
    let r_end: f64 = 1e12;
    let log_dr: f64 = 0.001;
    let surface_density: f64 = 1.0;

    let (start, end, k, gamma, model_type) = get_gas_parameters(gas_type, mu_e);

    let central_densities: Vec<f64> = (0..n)
        .map(|i: i32| {
            let t: f64 = i as f64 / (n - 1) as f64; // Linear interpolation factor
            10f64.powf(start + t * (end - start))
        })
        .collect();
    
    let start_density: f64 = central_densities[0];
    let end_density: f64 = central_densities[central_densities.len() - 1];
    let summary_file: String = format!("/Users/karankinariwala/Dropbox/KARAN/2 Areas/Education/PhD/gr_compact_objects/research_project/data/{gas_type}_summary_{start_density:.2e}_to_{end_density:.2e}.csv");
    let mut file = File::create(&summary_file)?;
    writeln!(file, "Central_Density (g/cm^3),Mass (g),Radius (cm),Pressure (dyn/cm^2)")?;

    let error_log_file = format!("/Users/karankinariwala/Dropbox/KARAN/2 Areas/Education/PhD/gr_compact_objects/research_project/data/{gas_type}_error_log.txt");
    let mut error_log = File::create(&error_log_file)?;

    for (i, &central_density) in central_densities.iter().enumerate() {

        let output_file: String = format!("/Users/karankinariwala/Dropbox/KARAN/2 Areas/Education/PhD/gr_compact_objects/research_project/data/{gas_type}_profile_central_density_{central_density:.2e}.csv");
        
        match run_compact_star_simulation(k, gamma, &model_type, r_start, r_end, central_density, log_dr, surface_density, &output_file) {
            Ok((final_radius, final_mass, final_pressure)) => {
                writeln!(file, "{},{},{},{}", central_density, final_mass, final_radius, final_pressure)?;

                println!("--------------------------------");
                println!("Running simulation {}/{}", i+1, central_densities.len());
                println!("Star Type: {}", gas_type);
                println!("Central Density: {:.4e} g cm^-3", central_density);
                println!("Final Mass: {:.4e} g (approx. {:.2e} solar masses)", final_mass, final_mass / 1.989e33);
                println!("Final Radius: {:.4e} cm (approx. {:.2e} solar radii)", final_radius, final_radius / 69.634e9);
                println!("Surface Pressure: {:.4e} dyn / cm^2", final_pressure);
                println!("Simulation Complete");
                println!("--------------------------------");
            }

            Err(e) => {
                let error_message = format!("Error in simulation for density {:.4e}: {:?}\n", central_density, e);
                eprintln!("{}", error_message);
                writeln!(error_log, "{}", error_message)?;
                continue;
            }
        }  
    }
    Ok(())
}

