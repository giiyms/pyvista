using System;

namespace PyVista.Core.Utilities;

/// <summary>
/// Static utility class for 3D transformation matrix operations.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.core.utilities.transformations</c>
/// module. All methods are pure matrix math with no VTK dependency.
/// </para>
/// </summary>
public static class Transformations
{
    /// <summary>
    /// Returns a new 4×4 identity matrix.
    /// </summary>
    /// <returns>A 4×4 identity matrix.</returns>
    public static double[,] Identity4x4()
    {
        var m = new double[4, 4];
        m[0, 0] = 1.0;
        m[1, 1] = 1.0;
        m[2, 2] = 1.0;
        m[3, 3] = 1.0;
        return m;
    }

    /// <summary>
    /// Computes a 4×4 rotation matrix for rotation about an arbitrary axis by the
    /// given angle using Rodrigues' rotation formula.
    /// <para>
    /// The rotation is counterclockwise when facing the direction of
    /// <paramref name="axis"/>.
    /// </para>
    /// </summary>
    /// <param name="axis">
    /// The direction vector of the rotation axis. It need not be a unit vector,
    /// but it must not be a zero vector.
    /// </param>
    /// <param name="angleDegrees">The rotation angle in degrees.</param>
    /// <param name="point">
    /// Optional origin of the rotation axis. When <c>null</c>, the rotation axis
    /// passes through the coordinate origin <c>(0, 0, 0)</c>.
    /// </param>
    /// <returns>A 4×4 rotation matrix.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="axis"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when <paramref name="axis"/> is a zero vector or does not have 3 elements.
    /// </exception>
    public static double[,] AxisAngleRotation(double[] axis, double angleDegrees, double[]? point = null)
    {
        ArgumentNullException.ThrowIfNull(axis);
        if (axis.Length != 3)
        {
            throw new ArgumentException("Axis must have exactly 3 elements.", nameof(axis));
        }

        double angle = angleDegrees * Math.PI / 180.0;

        // Return identity for zero rotation
        if (angle % (2.0 * Math.PI) == 0.0)
        {
            return Identity4x4();
        }

        // Normalize axis
        double norm = Math.Sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
        if (norm < 1e-15)
        {
            throw new ArgumentException("Cannot rotate around a zero vector axis.", nameof(axis));
        }

        double nx = axis[0] / norm;
        double ny = axis[1] / norm;
        double nz = axis[2] / norm;

        // Build the skew-symmetric cross-product matrix K
        // K = [[0, -nz, ny], [nz, 0, -nx], [-ny, nx, 0]]
        double sinA = Math.Sin(angle);
        double cosA = Math.Cos(angle);

        // Round for special angles (multiples of 90°)
        if (angle % (Math.PI / 2.0) == 0.0)
        {
            sinA = Math.Round(sinA);
            cosA = Math.Round(cosA);
        }

        // R = I + sin(a)*K + (1 - cos(a))*K^2  (Rodrigues' formula)
        double oneMinusCos = 1.0 - cosA;

        var result = Identity4x4();

        result[0, 0] = cosA + oneMinusCos * nx * nx;
        result[0, 1] = oneMinusCos * nx * ny - sinA * nz;
        result[0, 2] = oneMinusCos * nx * nz + sinA * ny;

        result[1, 0] = oneMinusCos * ny * nx + sinA * nz;
        result[1, 1] = cosA + oneMinusCos * ny * ny;
        result[1, 2] = oneMinusCos * ny * nz - sinA * nx;

        result[2, 0] = oneMinusCos * nz * nx - sinA * ny;
        result[2, 1] = oneMinusCos * nz * ny + sinA * nx;
        result[2, 2] = cosA + oneMinusCos * nz * nz;

        // If a point is provided, add the translation component:
        // b = point - R * point
        if (point != null)
        {
            if (point.Length != 3)
            {
                throw new ArgumentException("Point must have exactly 3 elements.", nameof(point));
            }

            double px = point[0], py = point[1], pz = point[2];
            double rpx = result[0, 0] * px + result[0, 1] * py + result[0, 2] * pz;
            double rpy = result[1, 0] * px + result[1, 1] * py + result[1, 2] * pz;
            double rpz = result[2, 0] * px + result[2, 1] * py + result[2, 2] * pz;

            result[0, 3] = px - rpx;
            result[1, 3] = py - rpy;
            result[2, 3] = pz - rpz;
        }

        return result;
    }

    /// <summary>
    /// Computes a 4×4 reflection matrix for reflection across a plane defined by
    /// the given normal vector.
    /// <para>
    /// The reflection plane passes through <paramref name="point"/> (or the origin
    /// when <paramref name="point"/> is <c>null</c>).
    /// </para>
    /// </summary>
    /// <param name="normal">
    /// The normal vector of the reflection plane. It need not be a unit vector,
    /// but it must not be a zero vector.
    /// </param>
    /// <param name="point">
    /// Optional point on the reflection plane. When <c>null</c>, the plane passes
    /// through the origin.
    /// </param>
    /// <returns>A 4×4 reflection matrix.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="normal"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when <paramref name="normal"/> is a zero vector or does not have 3 elements.
    /// </exception>
    public static double[,] ReflectionMatrix(double[] normal, double[]? point = null)
    {
        ArgumentNullException.ThrowIfNull(normal);
        if (normal.Length != 3)
        {
            throw new ArgumentException("Normal must have exactly 3 elements.", nameof(normal));
        }

        double norm = Math.Sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
        if (norm < 1e-15)
        {
            throw new ArgumentException("Cannot reflect across a zero normal vector.", nameof(normal));
        }

        double nx = normal[0] / norm;
        double ny = normal[1] / norm;
        double nz = normal[2] / norm;

        // R = I - 2 * n * n^T
        var result = Identity4x4();
        result[0, 0] = 1.0 - 2.0 * nx * nx;
        result[0, 1] = -2.0 * nx * ny;
        result[0, 2] = -2.0 * nx * nz;

        result[1, 0] = -2.0 * ny * nx;
        result[1, 1] = 1.0 - 2.0 * ny * ny;
        result[1, 2] = -2.0 * ny * nz;

        result[2, 0] = -2.0 * nz * nx;
        result[2, 1] = -2.0 * nz * ny;
        result[2, 2] = 1.0 - 2.0 * nz * nz;

        // b = point - R * point
        if (point != null)
        {
            if (point.Length != 3)
            {
                throw new ArgumentException("Point must have exactly 3 elements.", nameof(point));
            }

            double px = point[0], py = point[1], pz = point[2];
            double rpx = result[0, 0] * px + result[0, 1] * py + result[0, 2] * pz;
            double rpy = result[1, 0] * px + result[1, 1] * py + result[1, 2] * pz;
            double rpz = result[2, 0] * px + result[2, 1] * py + result[2, 2] * pz;

            result[0, 3] = px - rpx;
            result[1, 3] = py - rpy;
            result[2, 3] = pz - rpz;
        }

        return result;
    }

    /// <summary>
    /// Applies a 4×4 transformation matrix to an array of 3D points.
    /// <para>
    /// Points are stored as a flat array of length <c>N * 3</c>, where each
    /// consecutive triple <c>(x, y, z)</c> represents one point.
    /// </para>
    /// </summary>
    /// <param name="transformation">A 4×4 homogeneous transformation matrix.</param>
    /// <param name="points">
    /// A flat array of 3D points with length <c>N * 3</c>.
    /// </param>
    /// <returns>
    /// A new flat array of the same length containing the transformed points.
    /// </returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="transformation"/> or <paramref name="points"/>
    /// is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when <paramref name="points"/> length is not a multiple of 3, or
    /// when <paramref name="transformation"/> is not 4×4.
    /// </exception>
    public static double[] ApplyTransformationToPoints(double[,] transformation, double[] points)
    {
        ArgumentNullException.ThrowIfNull(transformation);
        ArgumentNullException.ThrowIfNull(points);

        if (transformation.GetLength(0) != 4 || transformation.GetLength(1) != 4)
        {
            throw new ArgumentException("Transformation must be a 4x4 matrix.", nameof(transformation));
        }

        if (points.Length % 3 != 0)
        {
            throw new ArgumentException("Points array length must be a multiple of 3.", nameof(points));
        }

        int n = points.Length / 3;
        var result = new double[points.Length];

        for (int i = 0; i < n; i++)
        {
            int offset = i * 3;
            double x = points[offset];
            double y = points[offset + 1];
            double z = points[offset + 2];

            result[offset] = transformation[0, 0] * x + transformation[0, 1] * y + transformation[0, 2] * z + transformation[0, 3];
            result[offset + 1] = transformation[1, 0] * x + transformation[1, 1] * y + transformation[1, 2] * z + transformation[1, 3];
            result[offset + 2] = transformation[2, 0] * x + transformation[2, 1] * y + transformation[2, 2] * z + transformation[2, 3];
        }

        return result;
    }

    /// <summary>
    /// Multiplies two 4×4 matrices and returns the product.
    /// <para>
    /// Computes <c>A × B</c>.
    /// </para>
    /// </summary>
    /// <param name="a">Left-hand 4×4 matrix.</param>
    /// <param name="b">Right-hand 4×4 matrix.</param>
    /// <returns>The 4×4 product matrix.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="a"/> or <paramref name="b"/> is <c>null</c>.
    /// </exception>
    public static double[,] MultiplyMatrices(double[,] a, double[,] b)
    {
        ArgumentNullException.ThrowIfNull(a);
        ArgumentNullException.ThrowIfNull(b);

        var result = new double[4, 4];
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                double sum = 0;
                for (int k = 0; k < 4; k++)
                {
                    sum += a[i, k] * b[k, j];
                }

                result[i, j] = sum;
            }
        }

        return result;
    }

    /// <summary>
    /// Computes the inverse of a 4×4 matrix using cofactor expansion.
    /// </summary>
    /// <param name="matrix">The 4×4 matrix to invert.</param>
    /// <returns>The 4×4 inverse matrix.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="matrix"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="InvalidOperationException">
    /// Thrown when the matrix is singular (determinant is zero).
    /// </exception>
    public static double[,] InvertMatrix(double[,] matrix)
    {
        ArgumentNullException.ThrowIfNull(matrix);

        // Flatten to a local array for easier indexing
        double[] m = new double[16];
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                m[i * 4 + j] = matrix[i, j];
            }
        }

        double[] inv = new double[16];

        inv[0] = m[5] * m[10] * m[15] - m[5] * m[11] * m[14] -
                 m[9] * m[6] * m[15] + m[9] * m[7] * m[14] +
                 m[13] * m[6] * m[11] - m[13] * m[7] * m[10];

        inv[4] = -m[4] * m[10] * m[15] + m[4] * m[11] * m[14] +
                  m[8] * m[6] * m[15] - m[8] * m[7] * m[14] -
                  m[12] * m[6] * m[11] + m[12] * m[7] * m[10];

        inv[8] = m[4] * m[9] * m[15] - m[4] * m[11] * m[13] -
                 m[8] * m[5] * m[15] + m[8] * m[7] * m[13] +
                 m[12] * m[5] * m[11] - m[12] * m[7] * m[9];

        inv[12] = -m[4] * m[9] * m[14] + m[4] * m[10] * m[13] +
                   m[8] * m[5] * m[14] - m[8] * m[6] * m[13] -
                   m[12] * m[5] * m[10] + m[12] * m[6] * m[9];

        inv[1] = -m[1] * m[10] * m[15] + m[1] * m[11] * m[14] +
                  m[9] * m[2] * m[15] - m[9] * m[3] * m[14] -
                  m[13] * m[2] * m[11] + m[13] * m[3] * m[10];

        inv[5] = m[0] * m[10] * m[15] - m[0] * m[11] * m[14] -
                 m[8] * m[2] * m[15] + m[8] * m[3] * m[14] +
                 m[12] * m[2] * m[11] - m[12] * m[3] * m[10];

        inv[9] = -m[0] * m[9] * m[15] + m[0] * m[11] * m[13] +
                  m[8] * m[1] * m[15] - m[8] * m[3] * m[13] -
                  m[12] * m[1] * m[11] + m[12] * m[3] * m[9];

        inv[13] = m[0] * m[9] * m[14] - m[0] * m[10] * m[13] -
                  m[8] * m[1] * m[14] + m[8] * m[2] * m[13] +
                  m[12] * m[1] * m[10] - m[12] * m[2] * m[9];

        inv[2] = m[1] * m[6] * m[15] - m[1] * m[7] * m[14] -
                 m[5] * m[2] * m[15] + m[5] * m[3] * m[14] +
                 m[13] * m[2] * m[7] - m[13] * m[3] * m[6];

        inv[6] = -m[0] * m[6] * m[15] + m[0] * m[7] * m[14] +
                  m[4] * m[2] * m[15] - m[4] * m[3] * m[14] -
                  m[12] * m[2] * m[7] + m[12] * m[3] * m[6];

        inv[10] = m[0] * m[5] * m[15] - m[0] * m[7] * m[13] -
                  m[4] * m[1] * m[15] + m[4] * m[3] * m[13] +
                  m[12] * m[1] * m[7] - m[12] * m[3] * m[5];

        inv[14] = -m[0] * m[5] * m[14] + m[0] * m[6] * m[13] +
                   m[4] * m[1] * m[14] - m[4] * m[2] * m[13] -
                   m[12] * m[1] * m[6] + m[12] * m[2] * m[5];

        inv[3] = -m[1] * m[6] * m[11] + m[1] * m[7] * m[10] +
                  m[5] * m[2] * m[11] - m[5] * m[3] * m[10] -
                  m[9] * m[2] * m[7] + m[9] * m[3] * m[6];

        inv[7] = m[0] * m[6] * m[11] - m[0] * m[7] * m[10] -
                 m[4] * m[2] * m[11] + m[4] * m[3] * m[10] +
                 m[8] * m[2] * m[7] - m[8] * m[3] * m[6];

        inv[11] = -m[0] * m[5] * m[11] + m[0] * m[7] * m[9] +
                   m[4] * m[1] * m[11] - m[4] * m[3] * m[9] -
                   m[8] * m[1] * m[7] + m[8] * m[3] * m[5];

        inv[15] = m[0] * m[5] * m[10] - m[0] * m[6] * m[9] -
                  m[4] * m[1] * m[10] + m[4] * m[2] * m[9] +
                  m[8] * m[1] * m[6] - m[8] * m[2] * m[5];

        double det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];
        if (Math.Abs(det) < 1e-15)
        {
            throw new InvalidOperationException("Matrix is singular and cannot be inverted.");
        }

        double invDet = 1.0 / det;
        var result = new double[4, 4];
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                result[i, j] = inv[i * 4 + j] * invDet;
            }
        }

        return result;
    }

    /// <summary>
    /// Computes the determinant of a 4×4 matrix.
    /// </summary>
    /// <param name="matrix">The 4×4 matrix.</param>
    /// <returns>The determinant value.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="matrix"/> is <c>null</c>.
    /// </exception>
    public static double Determinant(double[,] matrix)
    {
        ArgumentNullException.ThrowIfNull(matrix);
        double[] m = new double[16];
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                m[i * 4 + j] = matrix[i, j];
            }
        }

        double c0 = m[5] * m[10] * m[15] - m[5] * m[11] * m[14] -
                     m[9] * m[6] * m[15] + m[9] * m[7] * m[14] +
                     m[13] * m[6] * m[11] - m[13] * m[7] * m[10];

        double c4 = -m[4] * m[10] * m[15] + m[4] * m[11] * m[14] +
                      m[8] * m[6] * m[15] - m[8] * m[7] * m[14] -
                      m[12] * m[6] * m[11] + m[12] * m[7] * m[10];

        double c8 = m[4] * m[9] * m[15] - m[4] * m[11] * m[13] -
                     m[8] * m[5] * m[15] + m[8] * m[7] * m[13] +
                     m[12] * m[5] * m[11] - m[12] * m[7] * m[9];

        double c12 = -m[4] * m[9] * m[14] + m[4] * m[10] * m[13] +
                       m[8] * m[5] * m[14] - m[8] * m[6] * m[13] -
                       m[12] * m[5] * m[10] + m[12] * m[6] * m[9];

        return m[0] * c0 + m[1] * c4 + m[2] * c8 + m[3] * c12;
    }

    /// <summary>
    /// Transposes a 4×4 matrix.
    /// </summary>
    /// <param name="matrix">The 4×4 matrix to transpose.</param>
    /// <returns>The transposed 4×4 matrix.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="matrix"/> is <c>null</c>.
    /// </exception>
    public static double[,] Transpose(double[,] matrix)
    {
        ArgumentNullException.ThrowIfNull(matrix);
        var result = new double[4, 4];
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                result[i, j] = matrix[j, i];
            }
        }

        return result;
    }

    /// <summary>
    /// Builds a 4×4 translation matrix.
    /// </summary>
    /// <param name="x">Translation along the X axis.</param>
    /// <param name="y">Translation along the Y axis.</param>
    /// <param name="z">Translation along the Z axis.</param>
    /// <returns>A 4×4 translation matrix.</returns>
    public static double[,] TranslationMatrix(double x, double y, double z)
    {
        var m = Identity4x4();
        m[0, 3] = x;
        m[1, 3] = y;
        m[2, 3] = z;
        return m;
    }

    /// <summary>
    /// Builds a 4×4 uniform or non-uniform scaling matrix.
    /// </summary>
    /// <param name="sx">Scale factor along the X axis.</param>
    /// <param name="sy">Scale factor along the Y axis.</param>
    /// <param name="sz">Scale factor along the Z axis.</param>
    /// <returns>A 4×4 scaling matrix.</returns>
    public static double[,] ScalingMatrix(double sx, double sy, double sz)
    {
        var m = Identity4x4();
        m[0, 0] = sx;
        m[1, 1] = sy;
        m[2, 2] = sz;
        return m;
    }

    /// <summary>
    /// Checks whether two 4×4 matrices are equal within the specified tolerance.
    /// </summary>
    /// <param name="a">First matrix.</param>
    /// <param name="b">Second matrix.</param>
    /// <param name="tolerance">
    /// Maximum absolute difference allowed per element. Defaults to <c>1e-10</c>.
    /// </param>
    /// <returns>
    /// <c>true</c> if all corresponding elements differ by less than
    /// <paramref name="tolerance"/>; otherwise, <c>false</c>.
    /// </returns>
    public static bool MatricesAreClose(double[,] a, double[,] b, double tolerance = 1e-10)
    {
        ArgumentNullException.ThrowIfNull(a);
        ArgumentNullException.ThrowIfNull(b);

        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                if (Math.Abs(a[i, j] - b[i, j]) > tolerance)
                {
                    return false;
                }
            }
        }

        return true;
    }
}
