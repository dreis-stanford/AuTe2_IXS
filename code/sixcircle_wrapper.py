"""
Wrapper for sixcircle - handles all initialization issues
"""

import sys
import os
import config

class SixCircleWrapper:
    """Safe wrapper for sixcircle functionality"""
    
    def __init__(self):
        """Initialize - suppress config file errors"""
        
        # Change to a temp directory to avoid config file issues
        original_dir = os.getcwd()
        temp_dir = os.path.join(config.PROJECT_ROOT, '.sixcircle_temp')
        os.makedirs(temp_dir, exist_ok=True)
        os.chdir(temp_dir)
        
        try:
            # Import sixcircle
            import sixcircle
            self.sixcircle = sixcircle
            
            # Get functions
            from sixcircle import (
                setlat, setlambda, setfrozen, freeze,
                or0, or1, or_swap, or_check, setaz,
                ca, br, mv, wh,
                setlm, setlm_clear,
                ca_a, wmab,
                pa, UB,
                save, load
            )
            
            # Store them
            self.setlat = setlat
            self.setlambda = setlambda
            self.setfrozen = setfrozen
            self.freeze = freeze
            self.or0 = or0
            self.or1 = or1
            self.or_swap = or_swap
            self.or_check = or_check
            self.setaz = setaz
            self.ca = ca
            self.br = br
            self.mv = mv
            self.wh = wh
            self.setlm = setlm
            self.setlm_clear = setlm_clear
            self.ca_a = ca_a
            self.wmab = wmab
            self.pa = pa
            self.UB = UB
            self.save = save
            self.load = load
            
            self._initialized = False
            
            print("✓ Sixcircle wrapper created")
            print("  Call setup_for_ixs() to initialize")
            
        finally:
            # Go back to original directory
            os.chdir(original_dir)
    
    def setup_for_ixs(self, photon_energy_kev=21.747):
        """
        Setup IXS configuration
        
        Parameters:
        -----------
        photon_energy_kev : float
            Photon energy in keV
        """
        wavelength = 12.398 / photon_energy_kev
        
        self.setlambda(wavelength)
        self.setaz(1, 0, 0)
        self.setfrozen(456)
        
        self._initialized = True
        
        print(f"\n✓ Configured for IXS:")
        print(f"  Energy: {photon_energy_kev:.3f} keV")
        print(f"  Wavelength: {wavelength:.6f} Å")
        print(f"  Frozen: mu, gam, omega (456)")
    
    def setup_crystal(self, a, b, c, alpha=90, beta=90, gamma=90, name="Sample"):
        """Setup crystal lattice"""
        if not self._initialized:
            raise RuntimeError("Call setup_for_ixs() first!")
        
        self.setlat(a, b, c, alpha, beta, gamma, name)
        print(f"\n✓ Crystal '{name}' configured:")
        print(f"  a={a:.3f} Å, b={b:.3f} Å, c={c:.3f} Å")
        print(f"  α={alpha:.1f}°, β={beta:.1f}°, γ={gamma:.1f}°")

def init_sixcircle(photon_energy_kev=21.747):
    """Quick initialization"""
    sc = SixCircleWrapper()
    sc.setup_for_ixs(photon_energy_kev)
    return sc

if __name__ == "__main__":
    print("Testing sixcircle wrapper...")
    print("=" * 60)
    
    sc = init_sixcircle(21.747)
    
    print("\n" + "=" * 60)
    print("Testing with cubic crystal...")
    print("=" * 60)
    sc.setup_crystal(5.0, 5.0, 5.0, name="TestCrystal")
    
    print("\n" + "=" * 60)
    print("Calculating (1,0,0)...")
    print("=" * 60)
    sc.ca(1, 0, 0)
    
    print("\n" + "=" * 60)
    print("SUCCESS!")
    print("=" * 60)
