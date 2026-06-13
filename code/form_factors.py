"""
Atomic form factors for X-ray scattering
Matches MATLAB CalcAtomicfQ.m with option to use xraylib
"""

import numpy as np

# Try to import xraylib, fall back to built-in if not available
try:
    import xraylib
    HAS_XRAYLIB = True
except ImportError:
    HAS_XRAYLIB = False
    print("Warning: xraylib not available, using built-in Cromer-Mann coefficients")


# Module-level flag for form factor message
_FORM_FACTOR_PRINTED = False

class FormFactorCalculator:
    """
    Calculate atomic form factors f(Q)
    Matches MATLAB CalcAtomicfQ.m
    """
    
    # Cromer-Mann coefficients from International Tables for Crystallography, Vol. C
    # Format: [a1, a2, a3, a4, b1, b2, b3, b4, c]
    CROMER_MANN = {
        'Au': {  # Gold (Z=79)
            'a': [33.4718, 18.6655, 12.8373, 5.8908],
            'b': [2.1702, 0.2172, 16.5796, 64.9567],
            'c': 8.1302
        },
        'Te': {  # Tellurium (Z=52)
            'a': [19.9644, 19.0138, 6.1449, 2.5239],
            'b': [4.8174, 0.4209, 28.5284, 70.8403],
            'c': 4.3520
        },
        'Si': {  # Silicon (Z=14)
            'a': [6.2915, 3.0353, 1.9891, 1.5410],
            'b': [2.4386, 32.3337, 0.6785, 81.6937],
            'c': 1.1407
        },
        'Ge': {  # Germanium (Z=32)
            'a': [16.0816, 6.3747, 3.7068, 3.6830],
            'b': [2.8509, 0.2516, 11.4468, 54.7625],
            'c': 2.1313
        }
    }
    
    def __init__(self, use_xraylib=True):
        """
        Initialize form factor calculator

        Parameters:
        -----------
        use_xraylib : bool
            Use xraylib if available (default True)
        """
        self.use_xraylib = use_xraylib and HAS_XRAYLIB

        global _FORM_FACTOR_PRINTED
        if not _FORM_FACTOR_PRINTED:
            if self.use_xraylib:
                print("Using xraylib for form factors")
            else:
                print("Using built-in Cromer-Mann coefficients")
            _FORM_FACTOR_PRINTED = True
    
    def calc_f(self, Q, symbol, scale=1.0):
        """
        Calculate atomic form factor f0(Q)
        
        Matches MATLAB: f = CalcAtomicfQ(Q, symbol, 'scale', scale)
        
        Parameters:
        -----------
        Q : float or ndarray
            Scattering vector magnitude. 
            By default expects Q in units of sin(theta)/lambda [Angstrom^-1]
            If Q is passed as 4*pi*sin(theta)/lambda, set scale=4*pi
        symbol : str
            Element symbol ('Au', 'Te', 'Ge')
        scale : float
            Scale factor for Q (default 1.0)
            If Q is in units of 4π*sin(θ)/λ, use scale=4π
        
        Returns:
        --------
        f : float or ndarray
            Atomic form factor f0(s) where s = sin(θ)/λ
        """
        
        # Standardize Q to sin(theta)/lambda units
        s = np.asarray(Q) / scale
        
        if self.use_xraylib:
            return self._calc_f_xraylib(s, symbol)
        else:
            return self._calc_f_cromer_mann(s, symbol)
    
    def _calc_f_xraylib(self, s, symbol):
        """
        Calculate form factor using xraylib
        
        xraylib expects q in units of Angstrom^-1 where the momentum transfer is:
        Q = 4π*sin(θ)/λ = 4π*q_xraylib
        
        So if we have s = sin(θ)/λ, then q_xraylib = s/(4π)
        """
        
        # Get atomic number
        Z = xraylib.SymbolToAtomicNumber(symbol)
        
        # xraylib uses q in units where Q_total = 4π*q_xraylib
        # We have s = sin(θ)/λ, so need to divide by 4π
        # But wait - let's check the xraylib documentation...
        # Actually xraylib FF_Rayl expects momentum transfer / (4π)
        # So q_xraylib = |Q|/(4π) where |Q| = 4π*sin(θ)/λ
        # Therefore q_xraylib = sin(θ)/λ = s
        
        # NO WAIT - checking xraylib source, FF_Rayl expects:
        # q = sin(θ)/λ (in Angstrom^-1)
        # which is exactly what we have as 's'
        
        # Let me verify: xraylib documentation says FF_Rayl(Z, q)
        # where q is in Angstrom^-1
        # The scattering vector is Q = (4π/λ) sin(θ)
        # So Q = 4π * s where s = sin(θ)/λ
        # But xraylib wants just s = sin(θ)/λ
        
        # Calculate form factor
        if np.isscalar(s):
            f = xraylib.FF_Rayl(Z, s)
        else:
            f = np.array([xraylib.FF_Rayl(Z, s_val) for s_val in s])
        
        return f
    
    def _calc_f_cromer_mann(self, s, symbol):
        """
        Calculate form factor using Cromer-Mann coefficients
        
        f0(s) = c + sum_{i=1}^4 a_i * exp(-b_i * s^2)
        where s = sin(θ)/λ in Angstrom^-1
        """
        
        if symbol not in self.CROMER_MANN:
            raise ValueError(f"Cromer-Mann coefficients not available for {symbol}")
        
        coeffs = self.CROMER_MANN[symbol]
        a = coeffs['a']
        b = coeffs['b']
        c = coeffs['c']
        
        # Calculate f0(s)
        s_array = np.asarray(s)
        f = c * np.ones_like(s_array)
        
        for i in range(4):
            f += a[i] * np.exp(-b[i] * s_array**2)
        
        return f
    
    def add_element(self, symbol, a, b, c):
        """
        Add custom Cromer-Mann coefficients
        
        Parameters:
        -----------
        symbol : str
            Element symbol
        a : list of 4 floats
            a coefficients
        b : list of 4 floats
            b coefficients
        c : float
            constant term
        """
        if len(a) != 4 or len(b) != 4:
            raise ValueError("Must provide 4 a and 4 b coefficients")
        
        self.CROMER_MANN[symbol] = {
            'a': list(a),
            'b': list(b),
            'c': float(c)
        }
        
        print(f"Added Cromer-Mann coefficients for {symbol}")


def calc_form_factor(Q, symbol, energy_keV, scale=1.0, use_xraylib=True):
    """
    Full atomic form factor f(Q, E) = f0(Q) + f'(E) + i*Fii(E), where Fii is
    xraylib's anomalous-scattering imaginary part (Fi/Fii).

    A note on the sign of the imaginary part, since this is easy to get
    backwards: xraylib's Fii(Z,E) is the *negative* of the physically
    positive (absorptive) f''(E) -- i.e. -Fii(Z,E) is the quantity that
    matches the optical theorem applied to xraylib's own photoelectric
    cross section CS_Photo (verified in tests/test_form_factors.py). Using
    Fii directly (as done here, NOT -Fii) is nonetheless correct for this
    codebase, because code/ixs.py's structure factor uses the phase
    convention exp(-2*pi*i*Q.r) (the opposite of the crystallographic
    standard exp(+2*pi*i*Q.r)). Writing F_std for the structure factor built
    with the standard convention (f''>0, exp(+2*pi*i*Q.r)) and F_ours for
    ours (Fii<0, exp(-2*pi*i*Q.r)):

        F_ours = sum_j (f0_j+f'_j - i*f''_j) * exp(-i*theta_j)
               = sum_j conj(f0_j+f'_j + i*f''_j) * conj(exp(i*theta_j))
               = conj(F_std)

    so |F_ours|^2 == |F_std|^2 termwise, for any Q and any number of atoms.
    The two sign flips (xraylib's Fii vs. f'', and our phase convention vs.
    the crystallographic standard) cancel exactly, giving the physically
    correct (absorptive, not "gain") |F(Q)|^2 used for IXS_stokes,
    IXS_antistokes and the elastic |F|^2. Do not "fix" this by negating Fii
    without also flipping the phase convention in code/ixs.py.

    Parameters:
    -----------
    Q : float or ndarray
        Scattering vector magnitude (see `scale`)
    symbol : str
        Element symbol
    energy_keV : float
        Photon energy in keV, for the energy-dependent dispersion
        corrections f'(E) and Fii(E)
    scale : float
        Scale factor for Q (see FormFactorCalculator.calc_f)
    use_xraylib : bool
        Use xraylib if available (default True)

    Returns:
    --------
    f : complex or float
        f0(Q) + f'(E) + i*Fii(E) if xraylib is available, otherwise the
        real f0(Q) from the built-in Cromer-Mann coefficients (f'=f''=0).
    """
    calculator = FormFactorCalculator(use_xraylib=use_xraylib)
    f0 = calculator.calc_f(Q, symbol, scale=scale)

    if calculator.use_xraylib:
        Z = xraylib.SymbolToAtomicNumber(symbol)
        f1 = xraylib.Fi(Z, energy_keV)
        f2 = xraylib.Fii(Z, energy_keV)
        return f0 + f1 + 1j * f2

    return f0


# Convenience function matching MATLAB interface
def CalcAtomicfQ(Q, symbol, scale=1.0, use_xraylib=True):
    """
    Calculate atomic form factor - matches MATLAB CalcAtomicfQ.m
    
    Parameters:
    -----------
    Q : float or ndarray
        Scattering vector magnitude
    symbol : str
        Element symbol
    scale : float
        Scale factor (use 4π if Q is in units of 4π*sin(θ)/λ)
    use_xraylib : bool
        Use xraylib if available
    
    Returns:
    --------
    f : float or ndarray
        Form factor
    """
    calculator = FormFactorCalculator(use_xraylib=use_xraylib)
    return calculator.calc_f(Q, symbol, scale=scale)


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    
    print("=" * 60)
    print("Testing Form Factor Calculator")
    print("=" * 60)
    
    # Test at a single Q value
    # NOTE: Different conventions:
    # - Cromer-Mann and xraylib both use s = sin(θ)/λ
    # - Scattering vector magnitude is Q = 4π*sin(θ)/λ = 4π*s
    
    s_test = 0.5  # sin(theta)/lambda in Angstrom^-1
    Q_test = 4 * np.pi * s_test  # Full scattering vector
    
    print(f"\nTesting at s = sin(θ)/λ = {s_test:.3f} Å⁻¹")
    print(f"           Q = 4π·s      = {Q_test:.3f} Å⁻¹")
    print()
    
    for symbol in ['Au', 'Te']:
        # Both should use s directly
        f_builtin = CalcAtomicfQ(s_test, symbol, scale=1.0, use_xraylib=False)
        print(f"{symbol}: f = {f_builtin:.3f} (Cromer-Mann)")
        
        if HAS_XRAYLIB:
            f_xraylib = CalcAtomicfQ(s_test, symbol, scale=1.0, use_xraylib=True)
            print(f"{symbol}: f = {f_xraylib:.3f} (xraylib)")
            print(f"     Difference: {abs(f_builtin - f_xraylib):.3f}")
        print()
    
    # Also test if we pass Q (full scattering vector) with scale=4π
    print("\nTesting with Q = 4π·sin(θ)/λ and scale=4π:")
    for symbol in ['Au', 'Te']:
        f_builtin = CalcAtomicfQ(Q_test, symbol, scale=4*np.pi, use_xraylib=False)
        print(f"{symbol}: f = {f_builtin:.3f} (should match above)")
    
    # Plot form factors vs s = sin(theta)/lambda
    print("\nPlotting form factors vs sin(θ)/λ...")
    
    s_range = np.linspace(0, 2, 100)  # sin(theta)/lambda in Angstrom^-1
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Plot vs s
    for symbol in ['Au', 'Te']:
        f_cm = CalcAtomicfQ(s_range, symbol, use_xraylib=False)
        ax1.plot(s_range, f_cm, '-', label=f'{symbol} (Cromer-Mann)', linewidth=2)
        
        if HAS_XRAYLIB:
            f_xr = CalcAtomicfQ(s_range, symbol, use_xraylib=True)
            ax1.plot(s_range, f_xr, '--', label=f'{symbol} (xraylib)', linewidth=2, alpha=0.7)
    
    ax1.set_xlabel('s = sin(θ)/λ (Å⁻¹)', fontsize=12)
    ax1.set_ylabel('f(s)', fontsize=12)
    ax1.set_title('Atomic Form Factors vs sin(θ)/λ', fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot vs Q = 4π*s
    Q_range = 4 * np.pi * s_range
    for symbol in ['Au', 'Te']:
        f_cm = CalcAtomicfQ(Q_range, symbol, scale=4*np.pi, use_xraylib=False)
        ax2.plot(Q_range, f_cm, '-', label=f'{symbol}', linewidth=2)
    
    ax2.set_xlabel('Q = 4π·sin(θ)/λ (Å⁻¹)', fontsize=12)
    ax2.set_ylabel('f(Q)', fontsize=12)
    ax2.set_title('Atomic Form Factors vs |Q|', fontsize=14)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('form_factors_test.png', dpi=150)
    print("Saved: form_factors_test.png")
    
    plt.show(block=False)
    
    print("\n" + "=" * 60)
    print("Form factor calculator ready!")
    print("=" * 60)
