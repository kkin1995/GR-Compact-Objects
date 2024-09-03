// To Solve the Differential Equation of a Free-Falling Body.
// y'' = -g

use gr_compact_objects::rk4::{Derivatives, RK4Solver, State};
use std::fs::File;
use std::io::{Result, Write};

const G: f64 = 9.81; // ms-1

struct FreeFall;

impl Derivatives<f64> for FreeFall {
    fn derivatives(&self, _t: f64, state: &State<f64>) -> State<f64> {
        State {
            values: vec![state.values[1], -G],
        }
    }
}

fn main() -> Result<()> {
    let y0: f64 = 100.0; // meters
    let p0: f64 = 1.0; // Dropped from rest
    let mut state = State {
        values: vec![y0, p0],
    };
    let t_start = 0.0;
    let dt = 0.1;
    let t_end = 10.0;

    let solver = RK4Solver::new(FreeFall, dt);

    let mut file = File::create("free_fall_data.txt")?;

    writeln!(file, "t,y,p")?;

    let mut t = t_start;
    while t < t_end && state.values[0] > 0.0 {
        writeln!(file, "{},{},{}", t, state.values[0], state.values[1])?;
        state = solver.step(t, &state);
        t += dt;
    }

    Ok(())
}
