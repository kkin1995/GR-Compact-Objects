use std::clone::Clone;
use std::marker::PhantomData;

// Define a custom Absolute trait for calculating the absolute value.
pub trait Absolute {
    fn abs(self) -> Self;
}

// Implement the Absolute trait for f64.
impl Absolute for f64 {
    fn abs(self) -> Self {
        f64::abs(self)
    }
}

#[derive(Clone)]
pub struct State<T> {
    pub values: Vec<T>,
}

pub trait Derivatives<T> {
    fn derivatives(&self, t: T, state: &State<T>) -> State<T>;
}

pub struct AdaptiveRK4Solver<'a, T, D: Derivatives<T>> {
    derivatives: &'a D,
    _marker: PhantomData<T>,
}

impl<'a, T, D> AdaptiveRK4Solver<'a, T, D>
where
    T: Copy
        + std::ops::Add<Output = T>
        + std::ops::Mul<Output = T>
        + std::ops::Div<Output = T>
        + std::ops::Sub<Output = T>
        + PartialOrd
        + From<f64>
        + Absolute,
    D: Derivatives<T>,
{
    pub fn new(derivatives: &'a D) -> Self {
        AdaptiveRK4Solver { derivatives, _marker: PhantomData }
    }

    pub fn adaptive_step(&self, t: T, state: &State<T>, adaptive_dt: T, tolerance: T) -> (State<T>, T) {
        let full_step = self.step(t, state, adaptive_dt, T::from(1.0));

        let half_step_1 = self.step(t, state, adaptive_dt * T::from(0.5), T::from(1.0));
        let half_step_2 = self.step(t + adaptive_dt * T::from(0.5), &half_step_1, adaptive_dt * T::from(0.5), T::from(1.0));

        let error: T = full_step
            .values
            .iter()
            .zip(half_step_2.values.iter())
            .map(|(s_full, s_half)| (*s_full - *s_half).abs())
            .fold(T::from(0.0), |acc, e| acc + e);

        let mut new_adaptive_dt = adaptive_dt;

        if error > tolerance {
            new_adaptive_dt = adaptive_dt * T::from(0.5);
        } else if error < tolerance * T::from(0.1) {
            new_adaptive_dt = adaptive_dt * T::from(2.0);
        }

        (full_step, new_adaptive_dt)
    }

    pub fn step(&self, t: T, state: &State<T>, adaptive_dt: T, adaptive_tuner: T) -> State<T> {
        let dt  = adaptive_tuner * adaptive_dt;
        
        let k1 = self.derivatives.derivatives(t, state);
        let k2 = self.derivatives.derivatives(
            t + dt * T::from(0.5),
            &State {
                values: state
                    .values
                    .iter()
                    .zip(k1.values.iter())
                    .map(|(&s, &k)| s + k * dt * T::from(0.5))
                    .collect(),
            },
        );
        let k3 = self.derivatives.derivatives(
            t + dt * T::from(0.5),
            &State {
                values: state
                    .values
                    .iter()
                    .zip(k2.values.iter())
                    .map(|(&s, &k)| s + k * dt * T::from(0.5))
                    .collect(),
            },
        );
        let k4 = self.derivatives.derivatives(
            t + dt,
            &State {
                values: state
                    .values
                    .iter()
                    .zip(k3.values.iter())
                    .map(|(&s, &k)| s + k * dt)
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
                    *s + (dt / T::from(6.0))
                        * (k1 + k2 * T::from(2.0) + k3 * T::from(2.0) + k4)
                })
                .collect(),
        }
    }
}