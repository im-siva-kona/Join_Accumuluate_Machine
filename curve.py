from __future__ import annotations

import math
import hashlib
from dataclasses import dataclass
from typing import List, ClassVar, Final
from jam.ring_vrf.curve.glv import GLVSpecs

@dataclass(frozen=True)
class Curve:
    """
    Base implementation of an elliptic curve.
    
    This class provides the core functionality for elliptic curve operations,
    particularly focused on hash-to-curve operations as specified in the 
    IETF draft-irtf-cfrg-hash-to-curve.
    
    Attributes:
        PRIME_FIELD: The prime field characteristic
        ORDER: The order of the curve
        GENERATOR_X: X-coordinate of the generator point
        GENERATOR_Y: Y-coordinate of the generator point
        COFACTOR: The cofactor of the curve
        glv: GLV optimization parameters
        Z: The Z parameter for the curve
    """
    
    # Curve Parameters
    PRIME_FIELD: Final[int]
    ORDER: Final[int]
    GENERATOR_X: Final[int]
    GENERATOR_Y: Final[int]
    COFACTOR: Final[int]
    glv: GLVSpecs
    Z: Final[int]
    
    # Hash to Curve Parameters
    M: ClassVar[int] = 1  # Degree of field extension
    K: ClassVar[int] = 128  # Security parameter
    
    # Suite String Parameters
    SUITE_STRING: ClassVar[str] = "Bandersnatch_SHA-512_ELL2"
    DST: ClassVar[bytes] = f"ECVRF_Bandersnatch_XMD:SHA-512_ELL2_RO_{SUITE_STRING}".encode()
    
    def __post_init__(self) -> None:
        """Validate curve parameters after initialization."""
        if not self._validate_parameters():
            raise ValueError("Invalid curve parameters")
    
    def _validate_parameters(self) -> bool:
        """
        Validate the curve parameters.
        
        Returns:
            bool: True if parameters are valid, False otherwise
        """
        return (
            self.PRIME_FIELD > 0 and
            self.ORDER > 0 and
            0 <= self.GENERATOR_X < self.PRIME_FIELD and
            0 <= self.GENERATOR_Y < self.PRIME_FIELD and
            self.COFACTOR > 0
        )
    
    @property
    def L(self) -> int:
        """
        Calculate the length parameter for hash-to-field operations.
        
        Returns:
            int: The calculated L parameter
        """
        return math.ceil((math.ceil(math.log2(self.PRIME_FIELD)) + self.K) / 8)
    
    def hash_to_field(self, msg: bytes, count: int) -> List[int]:
        """
        Hash an arbitrary string to one or more field elements.
        
        Args:
            msg: The message to hash
            count: Number of field elements to generate
            
        Returns:
            List[int]: List of field elements
            
        Raises:
            ValueError: If count is negative or msg is None
        """
        if count < 0:
            raise ValueError("Count must be non-negative")
        if msg is None:
            raise ValueError("Message cannot be None")
            
        len_in_bytes = count * self.M * self.L
        uniform_bytes = self.expand_message_xmd(msg, len_in_bytes)
        
        u_values: List[int] = []
        for i in range(count):
            for j in range(self.M):
                elm_offset = self.L * (j + i * self.M)
                tv = uniform_bytes[elm_offset:elm_offset + self.L]
                e_j = int.from_bytes(tv, 'big') % self.PRIME_FIELD
                u_values.append(e_j)
                
        return u_values
    
    def expand_message_xmd(self, msg: bytes, len_in_bytes: int) -> bytes:
        """
        Expand a message using XMD (eXpandable Message Digest).
        
        Args:
            msg: The message to expand
            len_in_bytes: Desired length of the output in bytes
            
        Returns:
            bytes: The expanded message
            
        Raises:
            ValueError: If the input parameters are invalid
        """
        hash_fn = hashlib.sha512
        b_in_bytes = hash_fn().digest_size
        ell = math.ceil(len_in_bytes / b_in_bytes)
        
        if ell > 255 or len_in_bytes > 65535 or len(self.DST) > 255:
            raise ValueError("Invalid input size parameters")

        DST_prime = self.DST + self.I2OSP(len(self.DST), 1)

        Z_pad = self.I2OSP(0, self.L)

        l_i_b_str = self.I2OSP(len_in_bytes, 2)

        msg_prime = Z_pad + msg + l_i_b_str + self.I2OSP(0, 1) + DST_prime

        b_0 = hash_fn(msg_prime).digest()

        b_1 = hash_fn(b_0 + self.I2OSP(1, 1) + DST_prime).digest()

        b_values = [b_1]
        for i in range(2, ell + 1):
            b_i = hash_fn(self.strxor(b_0, b_values[-1]) + self.I2OSP(i, 1) + DST_prime).digest()
            b_values.append(b_i)

        uniform_bytes = b''.join(b_values)

        return uniform_bytes[:len_in_bytes]

    def mod_inverse(self, val: int) -> int:
        """
        Compute modular multiplicative inverse.
        
        Args:
            val: Value to invert
            field: Modulus
            
        Returns:
            int: Modular inverse
            
        Raises:
            ValueError: If inverse doesn't exist
        """
        if pow(val, self.PRIME_FIELD - 1, self.PRIME_FIELD) != 1:
            raise ValueError("No inverse exists")
        return pow(val, self.PRIME_FIELD - 2, self.PRIME_FIELD)
    
    @staticmethod
    def CMOV(a:int, b:int, cond:int)->int:
        """Constant-time conditional move: if cond is True, return b; else return a."""
        return b if cond else a

    @staticmethod
    def sgn0(x:int)->int:
        """Return the sign of x: 1 if odd, 0 if even."""
        return x % 2

    def find_z_ell2(self) -> int:
        return 5

    def is_square(self, val:int)->bool:
        if val == 0:
            return True
        return pow(val, (self.PRIME_FIELD - 1) // 2, self.PRIME_FIELD) == 1

    def mod_sqrt(self, val: int) -> int:
        """
        Compute the square root modulo prime field.
        
        Args:
            val: Value to compute square root of
            
        Returns:
            int: Square root of val modulo prime field
            
        Raises:
            ValueError: If no square root exists
        """
        if val == 0:
            return 0
            
        if not self.is_square(val):
            raise ValueError("No square root exists")
            
        # Tonelli-Shanks algorithm
        q = self.PRIME_FIELD - 1
        s = 0
        while q % 2 == 0:
            q //= 2
            s += 1
            
        if s == 1:
            return pow(val, (self.PRIME_FIELD + 1) // 4, self.PRIME_FIELD)
            
        # Find quadratic non-residue
        z = 2
        while self.is_square(z):
            z += 1
            
        m = s
        c = pow(z, q, self.PRIME_FIELD)
        t = pow(val, q, self.PRIME_FIELD)
        r = pow(val, (q + 1) // 2, self.PRIME_FIELD)
        
        while t != 1:
            i = 0
            temp = t
            while temp != 1:
                temp = (temp * temp) % self.PRIME_FIELD
                i += 1
                if i == m:
                    raise ValueError("No square root exists")
                    
            b = pow(c, 1 << (m - i - 1), self.PRIME_FIELD)
            m = i
            c = (b * b) % self.PRIME_FIELD
            t = (t * c) % self.PRIME_FIELD
            r = (r * b) % self.PRIME_FIELD
            
        return r

    @staticmethod
    def sha512(data: bytes) -> bytes:
        """Calculate SHA-512 hash"""
        return hashlib.sha512(data).digest()


    @staticmethod
    def I2OSP(value:int, length:int)->bytes:
        if value >= 256 ** length:
            raise ValueError("integer too large")
        return value.to_bytes(length, 'big')


    @staticmethod
    def OS2IP(octets:bytearray)->int:
        return int.from_bytes(octets, 'big')


    @staticmethod
    def strxor(s1:bytes, s2:bytes)->bytes:
        return bytes(a ^ b for a, b in zip(s1, s2))



