use pyo3::prelude::*;
use pyo3::types::PyList;
use pyo3::PyAny;

/// fortran dependencies
mod constants;
mod SRtoolkit;

/// Example retrieving constant from constants file can be improved by using macros
#[pyfunction]
fn get_pi() -> PyResult<f64> {
    Ok(constants::constants::CLIGHT)
}
/// functions from SRtoolkit
#[pyfunction]
fn bofg(py: Python, arg: &PyAny) -> PyResult<PyObject> {
    if let Ok(single_value) = arg.extract::<f64>() {
        // Call the scalar version of bofg
        Ok(SRtoolkit::srtoolkit::bofg_s(single_value).into_py(py))
    } else if let Ok(vec) = arg.extract::<Vec<f64>>() {
        // Call the vector version of bofg
        Ok(SRtoolkit::srtoolkit::bofg_v(&vec).into_py(py))
    } else {
        Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>(
            "Argument must be a float or a list of floats.",
        ))
    }
}


/// A Python module implemented in Rust.
#[pymodule]
fn pyparamo(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(get_pi, m)?)?;
    m.add_function(wrap_pyfunction!(bofg, m)?)?;
    Ok(())
}
