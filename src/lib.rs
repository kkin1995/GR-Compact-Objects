//! A crate for simulating compact objects like black holes, neutron stars, and white dwarfs
//! using general relativity.
//!
//! This crate provides different modules for simulating anything from the internal structure of 
//! compact objects to particle orbits around said compact objects. This project is currently in active
//! development and may change significantly from time to time.
//! 
//!
//! # Features
//!
//! - Generic implementation of the RK4 method
//! - Flexible `State` and `Derivatives` abstractions for various ODE systems
//! - Implementation of stellar structure simulations using Newtonian mechanics.
//!
//! # Example
//!
//! Here's a simple example of solving a harmonic oscillator using the RK4 solver:
//!
//! ```
//! use gr_compact_objects::(RKSolver, State, Derivatives);
//! struct HarmonicOscillator;
//!
//! impl Derivatives<f64> for HarmonicOscillator {
//! fn derivatives(&self, _t: f64, state: &State<f64>) -> State<f64> {
//!     State { values: vec![state.values[1], -state.values[0]] }
//!     }
//! }
//!
//! fn main() {
//!     let solver = RK4Solver::new(HarmonnicOscillator, 0.1);
//!     let mut state = State { values: vec![1.0, 0.0] };
//!     let mut t = 0.0;
//!
//!     for _ in 0..100 {
//!         state = solver.step(t, &state);
//!         t += 0.1;
//!         println!("t: {}, x: {}, v: {}", t, state.values[0], state.values[1]);
//!     }
//! }
//!
//! ```
//!
//! # Modules
//!
//! - [`rk4`]: Contains the implementation of the RK4 solver.
//! - [`stellar_structure`]: Contains the module to simulate the internal structure of compact objects for a given metric.
pub mod rk4;
pub mod adaptive_rk4;
pub mod stellar_structure;
