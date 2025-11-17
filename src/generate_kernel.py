#!/usr/bin/env python3
"""
Python function to write NumPy arrays to binary files in Fortran-compatible format.

Binary file format:
- Header: 2 int32 values (dimensions n1, n2)
- Data: n1 x n2 double precision (float64) values in column-major order
"""
import numpy as np
import struct


def write_binary_for_fortran(filepath, data):
    """
    Write a NumPy array to a binary file in Fortran-compatible format.
    
    The binary file will contain:
    - Header: 2 int32 values representing the array dimensions
    - Data: Double precision (float64) values in column-major (Fortran) order
    
    Parameters:
    -----------
    filepath : str
        Path to the output binary file
    data : numpy.ndarray
        2D NumPy array to write. Will be converted to float64 if needed.
        
    Returns:
    --------
    None
    
    Raises:
    -------
    ValueError : If data is not a 2D array
    
    Example:
    --------
    >>> import numpy as np
    >>> data = np.random.randn(10, 15)
    >>> write_binary_for_fortran('output.bin', data)
    """
    # Ensure data is a NumPy array
    data = np.asarray(data)
    
    # Check that data is 2D
    if data.ndim != 2:
        raise ValueError(f"Data must be 2D array, got {data.ndim}D array")
    
    # Get dimensions
    n1, n2 = data.shape
    
    # Convert to float64 (double precision) if not already
    data_f64 = data.astype(np.float64, copy=False)
    
    # Write to binary file
    with open(filepath, 'wb') as f:
        # Write header (2 int32 values)
        # Using native byte order (system default)
        f.write(struct.pack('ii', n1, n2))
        
        # Write data in column-major (Fortran) order
        # Transpose to convert from row-major (C) to column-major (Fortran)
        f.write(data_f64.T.tobytes())
    
    print(f"Successfully wrote binary file: {filepath}")
    print(f"  Dimensions: {n1} x {n2}")
    print(f"  File size: {2*4 + n1*n2*8} bytes")


def read_binary_from_fortran(filepath):
    """
    Read a binary file in Fortran format back into a NumPy array.
    
    This is useful for testing and verification.
    
    Parameters:
    -----------
    filepath : str
        Path to the binary file
        
    Returns:
    --------
    data : numpy.ndarray
        2D NumPy array with the data from the file
    """
    with open(filepath, 'rb') as f:
        # Read header (2 int32 values)
        n1, n2 = struct.unpack('ii', f.read(8))
        
        # Read data
        data_bytes = f.read(n1 * n2 * 8)
        
        # Convert to NumPy array and reshape
        # Data is in column-major order, so read with Fortran order
        data = np.frombuffer(data_bytes, dtype=np.float64).reshape((n2, n1)).T
    
    print(f"Successfully read binary file: {filepath}")
    print(f"  Dimensions: {n1} x {n2}")
    
    return data


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    r = np.linspace(-8, 8, 2048) ** 2
    r = np.add.outer(r,r)
    k = np.exp(-r)
    plt.imshow(k)
    plt.show()
    plt.imshow(k - 0.5*np.exp(-r/2))
    plt.show()
    write_binary_for_fortran("kernel1.bin", k)
    write_binary_for_fortran("kernel2.bin", k - 0.5*np.exp(-r/2))
    write_binary_for_fortran("kernel3.bin", k*0)
