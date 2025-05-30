from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple, Final
import numpy as np
from sympy import sqrt

@dataclass(frozen=True)
class GLVSpecs:
    """
    Gallant-Lambert-Vanstone endomorphism parameters.
    
    This class implements the GLV method for faster scalar multiplication
    on elliptic curves that admit an efficient endomorphism.
    
    Attributes:
        is_enabled: Whether GLV optimization is enabled
        lambda_param: Eigenvalue of the endomorphism
        constant_b: First decomposition constant
        constant_c: Second decomposition constant
    """
    is_enabled: Final[bool] = False
    lambda_param: Final[int] = 0
    constant_b: Final[int] = 0
    constant_c: Final[int] = 0
    
    def __post_init__(self) -> None:
        """Validate GLV parameters."""
        if self.is_enabled and not self._validate_parameters():
            raise ValueError("Invalid GLV parameters")
    
    def _validate_parameters(self) -> bool:
        """
        Validate GLV parameters.
        
        Returns:
            bool: True if parameters are valid
        """
        return (
            self.lambda_param != 0 and
            self.constant_b != 0 and
            self.constant_c != 0
        )
    
    def extended_euclidean_algorithm(
        self, 
        n: int, 
        lam: int
    ) -> List[Tuple[int, int, int]]:
        """
        Compute extended Euclidean algorithm sequence.
        
        Args:
            n: Curve order
            lam: Lambda parameter
            
        Returns:
            List[Tuple[int, int, int]]: Sequence of (s, t, r) values
        """
        if n <= 0 or lam <= 0:
            raise ValueError("Inputs must be positive")
            
        s0, t0, r0 = 1, 0, n
        s1, t1, r1 = 0, 1, lam
        sequence = [(s0, t0, r0), (s1, t1, r1)]
        
        while r1 != 0:
            q = r0 // r1
            r2 = r0 - q * r1
            s2 = s0 - q * s1
            t2 = t0 - q * t1
            sequence.append((s2, t2, r2))
            s0, t0, r0 = s1, t1, r1
            s1, t1, r1 = s2, t2, r2
            
        return sequence[:-1]
    
    def find_short_vectors(
        self, 
        n: int, 
        lam: int
    ) -> Tuple[Tuple[int, int], Tuple[int, int]]:
        """
        Find short vectors for scalar decomposition.
        
        Args:
            n: Curve order
            lam: Lambda parameter
            
        Returns:
            Tuple[Tuple[int, int], Tuple[int, int]]: Two shortest vectors
        """
        sequence = self.extended_euclidean_algorithm(n, lam)
        sqrt_n = int(sqrt(n))
        
        # Find largest index m where r_m >= sqrt(n)
        m = max(i for i, (_, _, r) in enumerate(sequence) if r >= sqrt_n)
        
        # Get components for v1
        rm_plus1, tm_plus1 = sequence[m + 1][2], sequence[m + 1][1]
        v1 = (rm_plus1, -tm_plus1)
        
        # Get components for v2
        if m + 2 < len(sequence):
            rm_plus2, tm_plus2 = sequence[m + 2][2], sequence[m + 2][1]
        else:
            rm_plus2, tm_plus2 = float('inf'), float('inf')
            
        v2_candidates = [
            (sequence[m][2], -sequence[m][1]),
            (rm_plus2, -tm_plus2)
        ]
        
        # Choose shorter vector
        v2 = min(v2_candidates, key=lambda v: v[0]**2 + v[1]**2)
        
        return v1, v2
    
    def decompose_scalar(
        self, 
        k: int, 
        n: int
    ) -> Tuple[int, int]:
        """
        Decompose scalar for faster multiplication.
        
        Args:
            k: Scalar to decompose
            n: Curve order
            
        Returns:
            Tuple[int, int]: Decomposed scalar components
            
        Raises:
            ValueError: If GLV is not enabled
        """
        if not self.is_enabled:
            return k, 0
            
        v1, v2 = self.find_short_vectors(n, self.lambda_param)
        v, b1, b2 = self._find_closest_lattice_point(k, v1, v2, n)
        
        k1 = (k - v[0]) % n
        k2 = (-v[1]) % n
        
        return k1, k2
    
    def _find_closest_lattice_point(
        self, 
        k: int, 
        v1: Tuple[int, int], 
        v2: Tuple[int, int], 
        n: int
    ) -> Tuple[Tuple[int, int], int, int]:
        """
        Find closest lattice point for scalar decomposition.
        
        Args:
            k: Scalar value
            v1: First basis vector
            v2: Second basis vector
            n: Curve order
            
        Returns:
            Tuple[Tuple[int, int], int, int]: Lattice point and coefficients
        """
        # Compute beta values
        beta1, beta2 = self._compute_beta(k, v1, v2, n)
        
        # Round to nearest integers
        b1 = round(beta1)
        b2 = round(beta2)
        
        # Compute lattice point
        v = (
            b1 * v1[0] + b2 * v2[0],
            b1 * v1[1] + b2 * v2[1]
        )
        
        return v, b1, b2
    
    def _compute_beta(
        self, 
        k: int, 
        v1: Tuple[int, int], 
        v2: Tuple[int, int], 
        n: int
    ) -> Tuple[float, float]:
        """
        Compute beta values for lattice point computation.
        
        Args:
            k: Scalar value
            v1: First basis vector
            v2: Second basis vector
            n: Curve order
            
        Returns:
            Tuple[float, float]: Beta values
        """
        det = v1[0] * v2[1] - v1[1] * v2[0]
        if det == 0:
            return 0.0, 0.0
            
        # Use numpy for better numerical stability
        A_inv = np.array([
            [v2[1], -v2[0]],
            [-v1[1], v1[0]]
        ]) / det
        
        beta = A_inv @ np.array([k, 0])
        return float(beta[0]), float(beta[1])

