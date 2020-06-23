use thiserror::Error;
#[derive(Error, Debug)]
#[error("DomainError: {value:?} is out of range {min:?} : {max:?}")]
pub struct DomainError {
    pub value: f64,
    pub min: f64,
    pub max: f64
}
