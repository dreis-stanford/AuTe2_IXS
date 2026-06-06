"""
FCC crystal structure coordinate transformations
Converts between conventional (cubic) and primitive (rhombohedral) bases

For Silicon and other FCC/diamond structures
"""

import numpy as np


class FCCStructure:
    """
    FCC (Face-Centered Cubic) structure transformations
    
    Conventional cell: cubic, a=b=c, α=β=γ=90°
    Primitive cell: rhombohedral, α=β=γ=60°
    """
    
    # Transformation matrices for reciprocal space (k-vectors)
    
    # Conventional to Primitive: k_prim = T * k_conv
    T_CONV_TO_PRIM = 0.5 * np.array([
        [0.0, 1.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 0.0]
    ])
    
    # Primitive to Conventional: k_conv = T_inv * k_prim
    T_PRIM_TO_CONV = np.array([
        [-1.0,  1.0,  1.0],
        [ 1.0, -1.0,  1.0],
        [ 1.0,  1.0, -1.0]
    ])
    
    @staticmethod
    def conv2prim_k(k_conv):
        """
        Convert wave vectors from conventional to primitive basis for FCC
        
        Parameters:
        -----------
        k_conv : ndarray
            Wave vector(s) in conventional cubic basis (units of 2π/a)
            Can be shape (3,), (N, 3), or (3, N)
        
        Returns:
        --------
        k_prim : ndarray
            Wave vector(s) in primitive reciprocal basis (same shape as input)
        
        Examples:
            X point [0, 1, 0] -> [0.5, 0, 0.5]
            X point [1, 0, 0] -> [0, 0.5, 0.5]
            L point [0.5, 0.5, 0.5] -> [0.5, 0.5, 0.5]
        """
        k_conv = np.asarray(k_conv)
        
        # Single vector
        if k_conv.shape == (3,):
            return FCCStructure.T_CONV_TO_PRIM @ k_conv
        
        # 2D array
        elif len(k_conv.shape) == 2:
            # Assume row format (N, 3) as default
            # This is more common in NumPy
            if k_conv.shape[1] == 3:
                return k_conv @ FCCStructure.T_CONV_TO_PRIM.T
            # Column format (3, N)
            elif k_conv.shape[0] == 3:
                return FCCStructure.T_CONV_TO_PRIM @ k_conv
            else:
                raise ValueError(f"Invalid shape: {k_conv.shape}")
        
        else:
            raise ValueError(f"Invalid shape: {k_conv.shape}")
    
    @staticmethod
    def prim2conv_k(k_prim):
        """
        Convert wave vectors from primitive to conventional basis for FCC
        
        Parameters:
        -----------
        k_prim : ndarray
            Wave vector(s) in primitive reciprocal basis
            Can be shape (3,), (N, 3), or (3, N)
        
        Returns:
        --------
        k_conv : ndarray
            Wave vector(s) in conventional cubic basis (units of 2π/a)
        
        Examples:
            [0.5, 0, 0.5] -> [0, 1, 0]
            [0, 0.5, 0.5] -> [1, 0, 0]
            [0.5, 0.5, 0.5] -> [0.5, 0.5, 0.5]
        """
        k_prim = np.asarray(k_prim)
        
        # Single vector
        if k_prim.shape == (3,):
            return FCCStructure.T_PRIM_TO_CONV @ k_prim
        
        # 2D array  
        elif len(k_prim.shape) == 2:
            # Assume row format (N, 3) as default
            if k_prim.shape[1] == 3:
                return k_prim @ FCCStructure.T_PRIM_TO_CONV.T
            # Column format (3, N)
            elif k_prim.shape[0] == 3:
                return FCCStructure.T_PRIM_TO_CONV @ k_prim
            else:
                raise ValueError(f"Invalid shape: {k_prim.shape}")
        
        else:
            raise ValueError(f"Invalid shape: {k_prim.shape}")


# Convenience functions matching MATLAB interface
def fcc_conv2prim_k(k_conv):
    """Convert from conventional to primitive coordinates for FCC"""
    return FCCStructure.conv2prim_k(k_conv)


def fcc_prim2conv_k(k_prim):
    """Convert from primitive to conventional coordinates for FCC"""
    return FCCStructure.prim2conv_k(k_prim)


if __name__ == "__main__":
    print("=" * 70)
    print("Testing FCC Coordinate Transformations")
    print("=" * 70)
    
    # Test 1: High-symmetry points
    print("\nTest 1: High-symmetry points")
    print("-" * 70)
    
    test_cases = [
        ("Γ", [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
        ("X", [0.0, 1.0, 0.0], [0.5, 0.0, 0.5]),
        ("X'", [1.0, 0.0, 0.0], [0.0, 0.5, 0.5]),
        ("L", [0.5, 0.5, 0.5], [0.5, 0.5, 0.5]),
    ]
    
    for label, k_conv_exp, k_prim_exp in test_cases:
        k_conv = np.array(k_conv_exp)
        k_prim = np.array(k_prim_exp)
        
        result_prim = fcc_conv2prim_k(k_conv)
        result_conv = fcc_prim2conv_k(k_prim)
        
        print(f"\n{label}: {k_conv} -> {result_prim}")
        assert np.allclose(result_prim, k_prim), f"Failed for {label}"
        assert np.allclose(result_conv, k_conv), f"Round-trip failed for {label}"
        print(f"  ✓ Pass")
    
    # Test 2: Array format
    print("\n" + "=" * 70)
    print("Test 2: Array of vectors")
    print("-" * 70)
    
    k_array = np.array([
        [0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.5, 0.5, 0.5],
        [1.0, 0.0, 0.0],
    ])
    
    expected = np.array([
        [0.0, 0.0, 0.0],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.5],
        [0.0, 0.5, 0.5],
    ])
    
    result = fcc_conv2prim_k(k_array)
    back = fcc_prim2conv_k(result)
    
    print("\nConventional -> Primitive:")
    for i in range(len(k_array)):
        print(f"  {k_array[i]} -> {result[i]}")
    
    assert np.allclose(result, expected), "Array transformation failed"
    assert np.allclose(back, k_array), "Round-trip failed"
    print("\n✓ Array transformation correct")
    print("✓ Round-trip verified")
    
    print("\n" + "=" * 70)
    print("All tests passed!")
    print("=" * 70)
