use gr_compact_objects::stellar_structure::{PhysicalConstants, IntegrationParams, StellarModel};
use std::fs::File;
use std::io::Write;

struct CompactStarParams {
    gamma: f64,
    k: f64,
    central_density: f64
}

impl CompactStarParams {
    fn new(gamma: f64, k: f64, central_density: f64) -> Self {
        Self { gamma, k, central_density }
    }
}

fn create_neutron_star_model(params: CompactStarParams) -> StellarModel {
    let constants = PhysicalConstants {
        k: params.k,
        gamma: params.gamma,
        g: 6.67430e-8 // gravitational constant in CGS
    };

    let integration_params = IntegrationParams {
        log_dr: 0.0001,
        surface_pressure_threshold: 1.1
    };

    StellarModel::new(
        constants,
        0.01,
        1e12,
        params.central_density,
        Some(integration_params)
    )
}

fn run_compact_star_simulation(params: CompactStarParams, output_file: &str) -> Result<(f64, f64, f64), gr_compact_objects::stellar_structure::StellarError> {
    let model: StellarModel = create_neutron_star_model(params);
    model.run(output_file)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let n: i32 = 11; // Replace this with the desired number of values
    let gas_type: &str = "non_relativistic_electron_gas";

    let (mut start, mut end, mut k, mut gamma): (f64, f64, f64, f64) = (0.0, 0.0, 0.0, 0.0);

    if gas_type == "non_relativistic_electron_gas" {
        // Non - Relativistic Electrons
        start = 3.0;
        end = 12.0; //5.0;
        k = 1.0036e13;
        gamma = 5.0 / 3.0;
    } else if gas_type == "relativistic_electron_gas" {
        // Relativisitic Electrons
        start = 7.0;
        end = 9.0;
        k = 1.2435e15;
        gamma = 4.0 / 3.0;
    } else if gas_type == "non_relativistic_neutron_gas" {
        // Non - Relativistic Neutrons
        start = 12.0;
        end = 14.0;
        k = 5.3802e9;
        gamma = 5.0 / 3.0;
    } else if gas_type == "relativistic_neutron_gas" {
        // Relativistic Neutrons
        start = 16.0;
        end = 18.0;
        k = 1.2293e15;
        gamma = 4.0 / 3.0;
    }

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
    writeln!(file, "Central_Density,Mass,Radius,Pressure")?;

    for (i, &central_density) in central_densities.iter().enumerate() {
        let params: CompactStarParams = CompactStarParams::new(
            gamma,
            k,
            central_density
        );

        let output_file: String = format!("/Users/karankinariwala/Dropbox/KARAN/2 Areas/Education/PhD/gr_compact_objects/research_project/data/{gas_type}_profile_central_density_{central_density:.2e}.csv");
        
        match run_compact_star_simulation(params, &output_file) {
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
                eprintln!("Error in simulation for density {}: {:?}", central_density, e);
                continue;
            }
        }  
    }
    Ok(())
}

