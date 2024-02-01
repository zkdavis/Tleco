use pyo3::prelude::*;
use pyo3::types::PyList;
use pyo3::PyAny;
use ndarray::ArrayView1;

/// fortran dependencies
mod constants;
mod SRtoolkit;
mod distribs;
mod misc;
mod specialf;

/// Example retrieving constant from constants file can be improved by using macros
#[pyfunction]
fn get_pi() -> PyResult<f64> {
    Ok(constants::CLIGHT)
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

///functions from distribs
#[pyfunction]
fn eq_59_park1995(t: f64, g: Vec<f64>) -> PyResult<Vec<f64>> {
    let g_array = ArrayView1::from(&g);
    let result = distribs::eq_59_park1995(t, g_array);
    Ok(result.to_vec())
}

/// A Python module implemented in Rust.
#[pymodule]
fn pyparamo(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(get_pi, m)?)?;
    m.add_function(wrap_pyfunction!(bofg, m)?)?;
    m.add_function(wrap_pyfunction!(eq_59_park1995, m)?)?;
    Ok(())
}
