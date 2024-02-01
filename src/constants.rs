// Rust module for constants

// Mathematical constants
pub const PI: f64 = 3.141592653589793238462643;
pub const HALF_PI: f64 = 1.570796326794896619231322;
pub const TWO_PI: f64 = 6.283185307179586476925287;
pub const FOUR_PI: f64 = 12.56637061435917295385057;
pub const SQRT_PI: f64 = 1.772453850905516027298167;
pub const SQRT_2: f64 = 1.414213562373095048801689;
pub const SQRT_3: f64 = 1.732050807568877293527446;
pub const LZERO: f64 = -460.517018598809136803598291; // log(1e-200)

// Physical constants
pub const CLIGHT: f64 = 2.99792458e10; // [cm / s]
pub const MASS_P: f64 = 1.67262192369e-24; // [g]
pub const MASS_E: f64 = 9.1093837015e-28; // [g]
pub const ECHARGE: f64 = 4.80320467299766e-10; // [cm^(3/2) g^(1/2) / s]
pub const SIGMAT: f64 = 6.6524587321000005e-25; // [cm^2]
pub const HPLANCK: f64 = 6.62607015e-27; // [erg s]
pub const HBAR: f64 = 1.0545718176461565e-27; // [erg s]
pub const KBOLTZ: f64 = 1.380649e-16; // [erg / K]
pub const SIGMASB: f64 = 5.6703744191844314e-5; // [erg / cm^2 / K^4 / s]
pub const GGRAV: f64 = 6.67430e-8; // [cm^3 / g / s^2]
pub const EVOLT: f64 = 1.602176634e-12; // [erg]
pub const NUCONST: f64 = 2.799248987233304e6; // eCharge / (2 * pi * m_e * cLight)
pub const JMBCONST: f64 = 6.6645698196351e-30; // sqrt(3) * eCharge^2 / (2 * cLight)
pub const AMBCONST: f64 = 3.6580794255805012e-3; // sqrt(3) * eCharge^2 / (4 * m_e * cLight)
pub const BCRITICAL: f64 = 4.414005218694872e13; // m_e^2 * cLight^3 / (eCharge * hbar)
pub const ENERGY_E: f64 = 8.187105776823886e-7; // m_e cLight^2
pub const ENERGY_P: f64 = 1.5032776159851257e-3; // m_p cLight^2
pub const MEC2_H: f64 = 1.235589963807414e20; // m_e c^2 / h
pub const H_MEC2: f64 = 6.53745089363765e-21; // h / m_e c^2

// Astronomical constants
pub const ASTRO_UNIT: f64 = 1495978707e13; // cm
pub const SOLAR_MASS: f64 = 1.9884754153381438e33; // g
pub const SOLAR_RADIUS: f64 = 6.957e10; // cm
pub const SOLAR_LUM: f64 = 3.828e33; // erg / s
pub const PARSEC: f64 = 3.085677581467192e18; // cm
pub const JANSKY: f64 = 1e-23; // erg / s cm^2 Hz
pub const LIGHTYR: f64 = 9.46073e17;
