//! Runge - Kutta 4th Order (RK4) solver for ordinary differential equations.
//!
//! This module provides a generic implementation of the RK4 method,
//! suitable for solving a wide range of initial value problems.

/// Represents the state of a system at a given time.
///
/// The state can contain any number of variables, each of type `T`.
pub struct State<T> {
    /// The values of the state variables.
    pub values: Vec<T>,
}

/// Trait for systems of differential equations.
///
/// Implementors of this trait define how to calculate derivatives
/// for a given state and time.
pub trait Derivatives<T> {
    /// Calculates the derivatives of the state variables.
    ///
    /// # Arguments
    ///
    /// * `t` - The current time
    /// * `state` - The current state of the system
    ///
    /// # Returns
    ///
    /// A new `State<T>` containing the calculated derivatives.
    fn derivatives(&self, t: T, state: &State<T>) -> State<T>;
}

/// RK4 solver for systems of ordinary differential equations.
pub struct RK4Solver<'a, T, D: Derivatives<T>> {
    derivatives: &'a D,
    dt: T,
}

impl<'a, T, D> RK4Solver<'a, T, D>
where
    T: Copy
        + std::ops::Add<Output = T>
        + std::ops::Mul<Output = T>
        + std::ops::Div<Output = T>
        + From<f64>,
    D: Derivatives<T>,
{
    /// Creates a new RK4Solver.
    ///
    /// # Arguments
    ///
    /// * `derivatives` - An implementation of the `Derivatives` trait
    /// * `dt` - The time step for the solver
    ///
    /// # Returns
    ///
    /// A new `RK4Solver` instance.
    pub fn new(derivatives: &'a D, dt: T) -> Self {
        RK4Solver { derivatives, dt }
    }

    /// Performs one step of the RK4 method.
    ///
    /// # Arguments
    ///
    /// * `t` - The current time
    /// * `state` - The cureent state of the system
    ///
    /// # Returns
    ///
    /// A new `State<T>` representing the system after one time step.
    ///
    /// # Example
    ///
    /// ```
    /// # use gr_compact_objects::(RK4Solver, State, Derivatives);
    /// # struct MyDerivatives;
    /// # impl Derivatives<f64> for MyDerivatives {
    ///     fn derivatives(&self, _t: f64, state: &State<f64>) -> State<f64> {
    ///         State { values: vec![-state.values[1], state.values[0]] }
    ///     }
    /// }
    /// let solver = RK4Solver::new(MyDerivatives, 0.1);
    /// let state = State { values: vec![1.0, 0.0] };
    /// let next_state = solver.step(0.0, &state);
    /// ```
    pub fn step(&self, t: T, state: &State<T>) -> State<T> {
        let k1 = self.derivatives.derivatives(t, state);
        let k2 = self.derivatives.derivatives(
            t + self.dt * T::from(0.5),
            &State {
                values: state
                    .values
                    .iter()
                    .zip(k1.values.iter())
                    .map(|(&s, &k)| s + k * self.dt * T::from(0.5))
                    .collect(),
            },
        );
        let k3 = self.derivatives.derivatives(
            t + self.dt * T::from(0.5),
            &State {
                values: state
                    .values
                    .iter()
                    .zip(k2.values.iter())
                    .map(|(&s, &k)| s + k * self.dt * T::from(0.5))
                    .collect(),
            },
        );
        let k4 = self.derivatives.derivatives(
            t + self.dt,
            &State {
                values: state
                    .values
                    .iter()
                    .zip(k3.values.iter())
                    .map(|(&s, &k)| s + k * self.dt)
                    .collect(),
            },
        );

        State {
            values: state
                .values
                .iter()
                .zip(k1.values.iter())
                .zip(k2.values.iter())
                .zip(k3.values.iter())
                .zip(k4.values.iter())
                .map(|((((s, &k1), &k2), &k3), &k4)| {
                    *s + (self.dt / T::from(6.0))
                        * (k1 + k2 * T::from(2.0) + k3 * T::from(2.0) + k4)
                })
                .collect(),
        }
    }
}
