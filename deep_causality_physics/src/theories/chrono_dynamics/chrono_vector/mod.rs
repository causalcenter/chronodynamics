//! Concrete vector type for chrono-gauge curvature and source operations.

use deep_causality_num::{Field, RealField};
use std::ops::{Add, Mul, Sub};

/// A concrete vector type for curvature and source operations.
///
/// This type is used as the concrete implementation behind the HKT
/// `RiemannMap` trait, which requires specific vector types for
/// tensor contraction operations.
///
/// # Safety Contract
///
/// When using `ChronoGaugeWitness::curvature()` or `::scatter()`,
/// the caller MUST pass `ChronoVector<T>` for all type parameters A, B, C, D.
/// Passing any other type results in undefined behavior.
#[derive(Debug, Clone, PartialEq)]
pub struct ChronoVector<T> {
    /// Vector components.
    pub data: Vec<T>,
}

impl<T> ChronoVector<T>
where
    T: Field + Copy,
{
    /// Creates a new chrono vector from a slice.
    #[inline]
    pub fn new(data: &[T]) -> Self {
        Self {
            data: data.to_vec(),
        }
    }

    /// Creates a zero vector of given dimension.
    #[inline]
    pub fn zeros(dim: usize) -> Self {
        Self {
            data: vec![T::zero(); dim],
        }
    }

    /// Creates a basis vector e_i.
    #[inline]
    pub fn basis(dim: usize, i: usize) -> Self {
        let mut data = vec![T::zero(); dim];
        if i < dim {
            data[i] = T::one();
        }
        Self { data }
    }

    /// Returns the dimension.
    #[inline]
    pub fn dim(&self) -> usize {
        self.data.len()
    }

    /// Returns a slice of the data.
    #[inline]
    pub fn as_slice(&self) -> &[T] {
        &self.data
    }

    /// Returns a mutable slice of the data.
    #[inline]
    pub fn as_mut_slice(&mut self) -> &mut [T] {
        &mut self.data
    }

    /// Computes the dot product with another vector.
    #[inline]
    pub fn dot(&self, other: &Self) -> T {
        self.data
            .iter()
            .zip(other.data.iter())
            .map(|(&a, &b)| a * b)
            .fold(T::zero(), |acc, x| acc + x)
    }

    /// Computes the squared magnitude.
    #[inline]
    pub fn magnitude_squared(&self) -> T {
        self.dot(self)
    }
}

impl<T: RealField> ChronoVector<T> {
    /// Computes the magnitude (Euclidean norm).
    #[inline]
    pub fn magnitude(&self) -> T {
        self.magnitude_squared().sqrt()
    }

    /// Returns a normalized copy of this vector.
    #[inline]
    pub fn normalized(&self) -> Self {
        let mag = self.magnitude();
        if mag == T::zero() {
            self.clone()
        } else {
            let inv_mag = T::one() / mag;
            Self {
                data: self.data.iter().map(|&x| x * inv_mag).collect(),
            }
        }
    }
}

impl<T> From<Vec<T>> for ChronoVector<T> {
    fn from(data: Vec<T>) -> Self {
        Self { data }
    }
}

impl<T> From<ChronoVector<T>> for Vec<T> {
    fn from(v: ChronoVector<T>) -> Self {
        v.data
    }
}

impl<T: Field + Copy> Add for ChronoVector<T> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let data: Vec<T> = self
            .data
            .iter()
            .zip(rhs.data.iter())
            .map(|(&a, &b)| a + b)
            .collect();
        Self { data }
    }
}

impl<T: Field + Copy> Sub for ChronoVector<T> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let data: Vec<T> = self
            .data
            .iter()
            .zip(rhs.data.iter())
            .map(|(&a, &b)| a - b)
            .collect();
        Self { data }
    }
}

impl<T: Field + Copy> Mul<T> for ChronoVector<T> {
    type Output = Self;

    fn mul(self, scalar: T) -> Self::Output {
        Self {
            data: self.data.iter().map(|&x| x * scalar).collect(),
        }
    }
}
