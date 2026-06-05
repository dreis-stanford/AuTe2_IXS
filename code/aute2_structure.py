"""
AuTe2 crystal structure and coordinate transformations
C2/m base-centered monoclinic structure

Matches MATLAB aute2_conv2prim_k.m and aute2_prim2conv_k.m
"""

import numpy as np


class AuTe2Structure:
    """
    Handle coordinate transformations for AuTe2 (C2/m monoclinic)
    
    AuTe2 has a base-centered monoclinic structure (C2/m)
    - Conventional cell: base-centered
    - Primitive cell: smaller, non-centered
    """
    
    # Transformation matrices
    # These are fixed for C2/m base-centered monoclinic
    
    # Conventional to Primitive (for k-vectors/reciprocal space)
    # [H, K, L]_conventional = T * [h, k, l]_primitive
    T_CONV_TO_PRIM = np.array([
        [ 0.5,  0.5,  0.0],
        [-0.5,  0.5,  0.0],
        [ 0.0,  0.0,  1.0]
    ])
    
    # Primitive to Conventional (for k-vectors/reciprocal space)
    # [h, k, l]_primitive = T_inv * [H, K, L]_conventional
    T_PRIM_TO_CONV = np.array([
        [ 1.0, -1.0,  0.0],
        [ 1.0,  1.0,  0.0],
        [ 0.0,  0.0,  1.0]
    ])
    
    @staticmethod
    def conv2prim_k(q_conv):
        """
        Convert reciprocal space vector from conventional to primitive
        
        Matches MATLAB: aute2_conv2prim_k.m
        
        For k-vectors (reciprocal space):
        [h, k, l]_prim = T * [H, K, L]_conv
        
        where:
        h = (H + K) / 2
        k = (-H + K) / 2
        l = L
        
        Parameters:
        -----------
        q_conv : ndarray
            Q-vector in conventional coordinates
            Can be shape (3,), (N, 3), or (3, N)
        
        Returns:
        --------
        q_prim : ndarray
            Q-vector in primitive coordinates (same shape as input)
        """
        q_conv = np.asarray(q_conv)
        original_shape = q_conv.shape
        
        # Handle different input shapes
        if q_conv.shape == (3,):
            # Single vector as column
            q_prim = AuTe2Structure.T_CONV_TO_PRIM @ q_conv
        
        elif len(q_conv.shape) == 2:
            if q_conv.shape[0] == 3:
                # Shape (3, N) - each column is a vector
                q_prim = AuTe2Structure.T_CONV_TO_PRIM @ q_conv
            elif q_conv.shape[1] == 3:
                # Shape (N, 3) - each row is a vector
                q_prim = q_conv @ AuTe2Structure.T_CONV_TO_PRIM.T
            else:
                raise ValueError(f"Invalid shape: {q_conv.shape}. Expected (3,), (N,3), or (3,N)")
        else:
            raise ValueError(f"Invalid shape: {q_conv.shape}. Expected (3,), (N,3), or (3,N)")
        
        return q_prim
    
    @staticmethod
    def prim2conv_k(q_prim):
        """
        Convert reciprocal space vector from primitive to conventional
        
        Matches MATLAB: aute2_prim2conv_k.m
        
        For k-vectors (reciprocal space):
        [H, K, L]_conv = T_inv * [h, k, l]_prim
        
        where:
        H = h - k
        K = h + k
        L = l
        
        Parameters:
        -----------
        q_prim : ndarray
            Q-vector in primitive coordinates
            Can be shape (3,), (N, 3), or (3, N)
        
        Returns:
        --------
        q_conv : ndarray
            Q-vector in conventional coordinates (same shape as input)
        """
        q_prim = np.asarray(q_prim)
        original_shape = q_prim.shape
        
        # Handle different input shapes
        if q_prim.shape == (3,):
            # Single vector as column
            q_conv = AuTe2Structure.T_PRIM_TO_CONV @ q_prim
        
        elif len(q_prim.shape) == 2:
            if q_prim.shape[0] == 3:
                # Shape (3, N) - each column is a vector
                q_conv = AuTe2Structure.T_PRIM_TO_CONV @ q_prim
            elif q_prim.shape[1] == 3:
                # Shape (N, 3) - each row is a vector
                q_conv = q_prim @ AuTe2Structure.T_PRIM_TO_CONV.T
            else:
                raise ValueError(f"Invalid shape: {q_prim.shape}. Expected (3,), (N,3), or (3,N)")
        else:
            raise ValueError(f"Invalid shape: {q_prim.shape}. Expected (3,), (N,3), or (3,N)")
        
        return q_conv
    
    @staticmethod
    def verify_transformations():
        """Verify that transformations are inverses"""
        
        # T_inv should be the inverse of T
        identity = AuTe2Structure.T_CONV_TO_PRIM @ AuTe2Structure.T_PRIM_TO_CONV
        
        if not np.allclose(identity, np.eye(3)):
            raise RuntimeError("Transformation matrices are not inverses!")
        
        print("✓ Transformation matrices verified as inverses")
        
        # Test round-trip
        q_test = np.array([0.5, 0.3, 0.2])
        q_conv = AuTe2Structure.prim2conv_k(q_test)
        q_back = AuTe2Structure.conv2prim_k(q_conv)
        
        if not np.allclose(q_test, q_back):
            raise RuntimeError("Round-trip transformation failed!")
        
        print("✓ Round-trip transformation verified")
        
        return True


# Convenience functions matching MATLAB interface
def aute2_conv2prim_k(q_conv):
    """
    Convert from conventional to primitive coordinates
    Matches MATLAB aute2_conv2prim_k.m
    """
    return AuTe2Structure.conv2prim_k(q_conv)


def aute2_prim2conv_k(q_prim):
    """
    Convert from primitive to conventional coordinates
    Matches MATLAB aute2_prim2conv_k.m
    """
    return AuTe2Structure.prim2conv_k(q_prim)


if __name__ == "__main__":
    print("=" * 70)
    print("Testing AuTe2 Coordinate Transformations")
    print("=" * 70)
    
    # Verify transformation matrices
    print("\nVerifying transformation matrices...")
    AuTe2Structure.verify_transformations()
    
    print("\n" + "=" * 70)
    print("Test 1: Single vectors")
    print("=" * 70)
    
    # Test cases from MATLAB
    test_cases = [
        ("Γ", np.array([0.0, 0.0, 0.0])),
        ("X", np.array([0.5, 0.0, 0.0])),
        ("M", np.array([0.5, 0.5, 0.0])),
        ("Z", np.array([0.0, 0.0, 0.5])),
    ]
    
    for label, q_prim in test_cases:
        q_conv = aute2_prim2conv_k(q_prim)
        q_back = aute2_conv2prim_k(q_conv)
        
        print(f"\n{label} point:")
        print(f"  Primitive:    [{q_prim[0]:.3f}, {q_prim[1]:.3f}, {q_prim[2]:.3f}]")
        print(f"  Conventional: [{q_conv[0]:.3f}, {q_conv[1]:.3f}, {q_conv[2]:.3f}]")
        print(f"  Round-trip:   [{q_back[0]:.3f}, {q_back[1]:.3f}, {q_back[2]:.3f}]")
        
        if not np.allclose(q_prim, q_back):
            print("  ⚠ Round-trip mismatch!")
        else:
            print("  ✓ Round-trip OK")
    
    print("\n" + "=" * 70)
    print("Test 2: Array of vectors (N, 3) format")
    print("=" * 70)
    
    # Array of vectors
    q_prim_array = np.array([
        [0.0, 0.0, 0.0],  # Γ
        [0.5, 0.0, 0.0],  # X
        [0.5, 0.5, 0.5],  # M
    ])
    
    q_conv_array = aute2_prim2conv_k(q_prim_array)
    q_back_array = aute2_conv2prim_k(q_conv_array)
    
    print("\nPrimitive coordinates:")
    print(q_prim_array)
    print("\nConventional coordinates:")
    print(q_conv_array)
    print("\nRound-trip:")
    print(q_back_array)
    
    if np.allclose(q_prim_array, q_back_array):
        print("\n✓ Array transformation verified")
    else:
        print("\n⚠ Array transformation failed!")
    
    print("\n" + "=" * 70)
    print("Test 3: Column vectors (3, N) format")
    print("=" * 70)
    
    q_prim_cols = q_prim_array.T  # Shape (3, N)
    q_conv_cols = aute2_prim2conv_k(q_prim_cols)
    q_back_cols = aute2_conv2prim_k(q_conv_cols)
    
    print(f"\nInput shape: {q_prim_cols.shape}")
    print(f"Output shape: {q_conv_cols.shape}")
    
    if np.allclose(q_prim_cols, q_back_cols):
        print("✓ Column format transformation verified")
    else:
        print("⚠ Column format transformation failed!")
    
    print("\n" + "=" * 70)
    print("Test 4: Example from your MATLAB code")
    print("=" * 70)
    
    # From your interactive code: Q_conv = [0.5, 0.5, 0] -> Q_prim = [0.5, 0, 0]
    Q_conv_example = np.array([0.5, 0.5, 0.0])
    Q_prim_example = aute2_conv2prim_k(Q_conv_example)
    
    print(f"\nConventional: {Q_conv_example}")
    print(f"Primitive:    {Q_prim_example}")
    print(f"Expected:     [0.5, 0.0, 0.0]")
    
    if np.allclose(Q_prim_example, [0.5, 0.0, 0.0]):
        print("✓ Matches MATLAB example!")
    else:
        print("⚠ Does not match MATLAB!")
    
    print("\n" + "=" * 70)
    print("All coordinate transformation tests complete!")
    print("=" * 70)
