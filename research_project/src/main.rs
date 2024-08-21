use gr_compact_objects::stellar_structure::{PhysicalConstants, IntegrationParams, StellarModel};
use std::fs::File;
use std::io::Write;

struct NeutronStarParams {
    gamma: f64,
    k: f64,
    central_density: f64
}

impl NeutronStarParams {
    fn new(gamma: f64, k: f64, central_density: f64) -> Self {
        Self { gamma, k, central_density }
    }
}

fn create_neutron_star_model(params: NeutronStarParams) -> StellarModel {
    let constants = PhysicalConstants {
        k: params.k,
        gamma: params.gamma,
        g: 6.67430e-8 // gravitational constant in CGS
    };

    let integration_params = IntegrationParams {
        log_dr: 0.001,
        surface_pressure_threshold: 1.01
    };

    StellarModel::new(
        constants,
        1.0,
        1e7,
        params.central_density,
        Some(integration_params)
    )
}

fn run_neutron_star_simulation(params: NeutronStarParams, output_file: &str) -> Result<(f64, f64, f64), gr_compact_objects::stellar_structure::StellarError> {
    let model = create_neutron_star_model(params);
    model.run(output_file)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let central_densities: Vec<f64> = (140..=150).map(|i| 10f64.powf(i as f64 / 10.0)).collect();
    let mut file = File::create("central_density_enumeration.csv")?;
    writeln!(file, "Central_Density,Mass,Radius,Pressure")?;

    for (i, &central_density) in central_densities.iter().enumerate() {
        println!("Running simulation {}/{}", i+1, central_densities.len());
        let params = NeutronStarParams::new(
            5.0 / 3.0,
            5.3802e9,
            central_density
        );

        let output_file: String = format!("neutron_star_profile_central_density_{central_density:.2e}.csv");
        
        match run_neutron_star_simulation(params, &output_file) {
            Ok((final_radius, final_mass, final_pressure)) => {
                writeln!(file, "{},{},{},{}", central_density, final_mass, final_radius, final_pressure)?;

                println!("--------------------------------");
                println!("Neutron Star Simulation Complete");
                println!("--------------------------------");
                println!("Central Density: {:.4e} g cm^-3", central_density);
                println!("Final Mass: {:.4e} g (approx. {:.2e} solar masses)", final_mass, final_mass / 1.989e33);
                println!("Final Radius: {:.4e} cm (approx. {:.2e} solar radii)", final_radius, final_radius / 69.634e9);
                println!("Surface Pressure: {:.4e} dyn / cm^2", final_pressure);
            }

            Err(e) => {
                eprintln!("Error in simulation for density {}: {:?}", central_density, e);
                continue;
            }
        }



        
    }

    Ok(())
}

