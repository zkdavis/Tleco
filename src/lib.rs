use pyo3::prelude::*;
use pyo3::types::PyList;
use pyo3::PyAny;
use ndarray::{Array1, ArrayView1};

/// fortran dependencies
mod constants;
mod SRtoolkit;
mod distribs;
mod misc;
mod specialf;
mod radiation;
mod pwl_integ;

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

#[pyfunction]
fn fp_findif_difu(dt_in: f64, g: Vec<f64>, nin: Vec<f64>, gdot_in: Vec<f64>,
              din: Vec<f64>, qin: Vec<f64>, tesc_in: f64, tlc: f64,
              check_params: Option<bool>) -> PyResult<Vec<f64>>{
    let g_array = ArrayView1::from(&g).to_owned();
    let nin_array = ArrayView1::from(&nin).to_owned();
    let gdot_in_array = ArrayView1::from(&gdot_in).to_owned();
    let din_array = ArrayView1::from(&din).to_owned();
    let qin_array = ArrayView1::from(&qin).to_owned();
    let result = distribs::fp_findif_difu(dt_in, &g_array,&nin_array,&gdot_in_array,&din_array,&qin_array,tesc_in,tlc,check_params);
    Ok(result.to_vec())
    }


#[pyfunction]
fn syn_emissivity_full(freqs: Vec<f64>, g: Vec<f64>, n: Vec<f64>, b: f64, with_abs: bool) -> PyResult<(Vec<f64>, Vec<f64>)> {
    let freqs_arr = Array1::from(freqs);
    let g_arr = Array1::from(g);
    let n_arr = Array1::from(n);

    let (jmbs, ambs) = radiation::syn_emissivity_full(&freqs_arr, &g_arr, &n_arr, b, with_abs);
    Ok((jmbs.to_vec(), ambs.to_vec()))
}


#[pyfunction]
pub fn rad_trans_blob(R: f64, jnu: Vec<f64>, anu: Vec<f64>) -> PyResult<Vec<f64>> {
    let jnu_arr = Array1::from_vec(jnu);
    let anu_arr = Array1::from_vec(anu);
    let result =  radiation::rad_trans_blob(R, &jnu_arr, &anu_arr);

    Ok(result.to_vec())
}


#[pyfunction]
fn ic_iso_powlaw_full(freqs: Vec<f64>, inu: Vec<f64>, g: Vec<f64>, n: Vec<f64>) -> PyResult<Vec<f64>> {
    let freqs_array = Array1::from_vec(freqs);
    let inu_array = Array1::from_vec(inu);
    let g_array = Array1::from_vec(g);
    let n_array = Array1::from_vec(n);

    let result = radiation::ic_iso_powlaw_full(&freqs_array, &inu_array, &g_array, &n_array);

    Ok(result.to_vec())
}


#[pyfunction]
fn ic_iso_monochrom_full(freqs: Vec<f64>, uext: f64,nuext: f64, n: Vec<f64>, g: Vec<f64>) -> PyResult<Vec<f64>> {
    let freqs_array = Array1::from_vec(freqs);
    let g_array = Array1::from_vec(g);
    let n_array = Array1::from_vec(n);
    let result = radiation::ic_iso_monochrom_full(&freqs_array, uext, nuext, &g_array, &n_array);
    Ok(result.to_vec())
}



/// A Python module implemented in Rust.
#[pymodule]
fn pyparamo(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(get_pi, m)?)?;
    m.add_function(wrap_pyfunction!(bofg, m)?)?;
    m.add_function(wrap_pyfunction!(eq_59_park1995, m)?)?;
    m.add_function(wrap_pyfunction!(fp_findif_difu, m)?)?;
    m.add_function(wrap_pyfunction!(syn_emissivity_full, m)?)?;
    m.add_function(wrap_pyfunction!(rad_trans_blob, m)?)?;
    m.add_function(wrap_pyfunction!(ic_iso_powlaw_full, m)?)?;
    m.add_function(wrap_pyfunction!(ic_iso_monochrom_full, m)?)?;
    Ok(())
}
