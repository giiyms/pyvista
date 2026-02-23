using System;
using System.Collections.Generic;
using System.Linq;

namespace PyVista.Core.Utilities;

/// <summary>
/// Static utility methods for point, line, and vector operations.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.core.utilities.points</c> module.
/// It provides helpers for constructing line segments, triangle meshes, fitting planes
/// and lines to point sets, and computing principal axes — all without a VTK dependency.
/// </para>
/// </summary>
public static class PointUtils
{
    /// <summary>
    /// Generates triangle-strip connectivity from a flat point-index array.
    /// <para>
    /// Triangle strips encode a sequence of connected triangles where each successive
    /// triangle shares an edge with the previous one. The output is a VTK-style
    /// connectivity array prefixed with the strip length.
    /// </para>
    /// </summary>
    /// <param name="pointIndices">
    /// An ordered sequence of point indices forming the triangle strip.
    /// Must have at least 3 elements.
    /// </param>
    /// <returns>
    /// A VTK-style strip connectivity array: <c>[n, idx0, idx1, ..., idxN-1]</c>.
    /// </returns>
    /// <exception cref="ArgumentException">
    /// Thrown when fewer than 3 point indices are provided.
    /// </exception>
    public static int[] MakeTristripConnectivity(int[] pointIndices)
    {
        ArgumentNullException.ThrowIfNull(pointIndices);

        if (pointIndices.Length < 3)
        {
            throw new ArgumentException(
                "At least 3 point indices are required to form a triangle strip.");
        }

        var result = new int[pointIndices.Length + 1];
        result[0] = pointIndices.Length;
        Array.Copy(pointIndices, 0, result, 1, pointIndices.Length);
        return result;
    }

    /// <summary>
    /// Generates non-connected line segments from an array of point coordinates.
    /// <para>
    /// Points are assumed to come in pairs: every two consecutive points form one
    /// line segment. An even number of points is required.
    /// </para>
    /// </summary>
    /// <param name="points">
    /// Flat array of point coordinates (length must be a multiple of 6, i.e., pairs of 3D points).
    /// </param>
    /// <returns>
    /// A tuple of (<c>points</c>, <c>lines</c>) where <c>lines</c> is the VTK-style line
    /// connectivity array.
    /// </returns>
    /// <exception cref="ArgumentException">
    /// Thrown when an odd number of points is provided or the array length is not a multiple of 3.
    /// </exception>
    public static (double[] Points, int[] Lines) LineSegmentsFromPoints(double[] points)
    {
        ArgumentNullException.ThrowIfNull(points);

        if (points.Length % 3 != 0)
        {
            throw new ArgumentException("Points array length must be a multiple of 3.");
        }

        int nPoints = points.Length / 3;
        if (nPoints % 2 != 0)
        {
            throw new ArgumentException(
                "An even number of points must be given to define each segment.");
        }

        int nLines = nPoints / 2;
        var lines = new int[nLines * 3];
        for (int i = 0; i < nLines; i++)
        {
            lines[i * 3]     = 2;
            lines[i * 3 + 1] = i * 2;
            lines[i * 3 + 2] = i * 2 + 1;
        }

        return ((double[])points.Clone(), lines);
    }

    /// <summary>
    /// Makes a connected polyline from an ordered array of 3D points.
    /// <para>
    /// Each consecutive pair of points defines one line segment, producing
    /// <c>N-1</c> segments for <c>N</c> points.
    /// </para>
    /// </summary>
    /// <param name="points">
    /// Flat array of point coordinates (length must be a multiple of 3).
    /// </param>
    /// <param name="close">
    /// When <c>true</c>, an additional segment connecting the last point back to the
    /// first is appended, creating a closed loop.
    /// </param>
    /// <returns>
    /// A tuple of (<c>points</c>, <c>lines</c>) where <c>lines</c> is the VTK-style
    /// connectivity array for the polyline segments.
    /// </returns>
    /// <exception cref="ArgumentException">
    /// Thrown when the points array is too short or its length is not a multiple of 3.
    /// </exception>
    public static (double[] Points, int[] Lines) LinesFromPoints(double[] points, bool close = false)
    {
        ArgumentNullException.ThrowIfNull(points);

        if (points.Length % 3 != 0)
        {
            throw new ArgumentException("Points array length must be a multiple of 3.");
        }

        int nPoints = points.Length / 3;
        if (nPoints < 2)
        {
            throw new ArgumentException("At least 2 points are required to form a line.");
        }

        int segmentCount = nPoints - 1 + (close ? 1 : 0);
        var lines = new int[segmentCount * 3];

        for (int i = 0; i < nPoints - 1; i++)
        {
            lines[i * 3]     = 2;
            lines[i * 3 + 1] = i;
            lines[i * 3 + 2] = i + 1;
        }

        if (close)
        {
            int last = (nPoints - 1) * 3;
            lines[last]     = 2;
            lines[last + 1] = nPoints - 1;
            lines[last + 2] = 0;
        }

        return ((double[])points.Clone(), lines);
    }

    /// <summary>
    /// Constructs triangle cells from an Nx3 array of triangle face indices
    /// and returns VTK-style connectivity.
    /// <para>
    /// Each face array element contains 3 point indices defining a triangle.
    /// </para>
    /// </summary>
    /// <param name="points">Flat array of point coordinates (Nx3 packed).</param>
    /// <param name="faces">
    /// Flat array of triangle face indices. Length must be a multiple of 3.
    /// Each group of 3 values represents one triangle.
    /// </param>
    /// <returns>
    /// A tuple of (<c>points</c>, <c>cells</c>) where <c>cells</c> is the VTK-style
    /// cell connectivity array with padding.
    /// </returns>
    /// <exception cref="ArgumentException">
    /// Thrown when points or faces have invalid lengths.
    /// </exception>
    public static (double[] Points, int[] Cells) MakeTriMesh(double[] points, int[] faces)
    {
        ArgumentNullException.ThrowIfNull(points);
        ArgumentNullException.ThrowIfNull(faces);

        if (points.Length % 3 != 0)
        {
            throw new ArgumentException("Points array length must be a multiple of 3.");
        }

        if (faces.Length % 3 != 0)
        {
            throw new ArgumentException("Faces array length must be a multiple of 3.");
        }

        int nFaces = faces.Length / 3;
        var cells = new int[nFaces * 4];
        for (int i = 0; i < nFaces; i++)
        {
            cells[i * 4]     = 3;
            cells[i * 4 + 1] = faces[i * 3];
            cells[i * 4 + 2] = faces[i * 3 + 1];
            cells[i * 4 + 3] = faces[i * 3 + 2];
        }

        return ((double[])points.Clone(), cells);
    }

    /// <summary>
    /// Creates vector field data from origin points and their associated vectors.
    /// <para>
    /// Returns parallel arrays of origins, vectors, and magnitudes suitable for
    /// glyph-based visualization.
    /// </para>
    /// </summary>
    /// <param name="origins">Flat array of origin coordinates (Nx3 packed).</param>
    /// <param name="vectors">Flat array of vector components (Nx3 packed).</param>
    /// <returns>
    /// A tuple of (<c>origins</c>, <c>vectors</c>, <c>magnitudes</c>).
    /// </returns>
    /// <exception cref="ArgumentException">
    /// Thrown when origins and vectors have different lengths or are not multiples of 3.
    /// </exception>
    public static (double[] Origins, double[] Vectors, double[] Magnitudes) VectorPolyData(
        double[] origins, double[] vectors)
    {
        ArgumentNullException.ThrowIfNull(origins);
        ArgumentNullException.ThrowIfNull(vectors);

        if (origins.Length != vectors.Length)
        {
            throw new ArgumentException("Origins and vectors must have the same length.");
        }

        if (origins.Length % 3 != 0)
        {
            throw new ArgumentException("Array lengths must be multiples of 3.");
        }

        int nPoints = origins.Length / 3;
        var magnitudes = new double[nPoints];
        for (int i = 0; i < nPoints; i++)
        {
            int off = i * 3;
            double vx = vectors[off], vy = vectors[off + 1], vz = vectors[off + 2];
            magnitudes[i] = Math.Sqrt(vx * vx + vy * vy + vz * vz);
        }

        return ((double[])origins.Clone(), (double[])vectors.Clone(), magnitudes);
    }

    /// <summary>
    /// Computes the principal axes of a set of 3D points using eigendecomposition
    /// of the covariance matrix.
    /// <para>
    /// The returned axes are orthonormal row vectors ordered by decreasing variance.
    /// The first axis explains the most variance. If <paramref name="returnStd"/> is
    /// <c>true</c>, the standard deviations along each axis are also returned.
    /// </para>
    /// </summary>
    /// <param name="points">Flat array of point coordinates (Nx3 packed).</param>
    /// <param name="returnStd">
    /// When <c>true</c>, also computes the standard deviation along each principal axis.
    /// </param>
    /// <returns>
    /// A tuple of (<c>axes</c>, <c>standardDeviations</c>). The <c>axes</c> array is
    /// a 3×3 row-major matrix (9 elements). <c>standardDeviations</c> is <c>null</c>
    /// when <paramref name="returnStd"/> is <c>false</c>.
    /// </returns>
    /// <exception cref="ArgumentException">
    /// Thrown when the points array has fewer than 3 points or is not a multiple of 3.
    /// </exception>
    public static (double[] Axes, double[]? StandardDeviations) PrincipalAxes(
        double[] points, bool returnStd = false)
    {
        ArgumentNullException.ThrowIfNull(points);

        if (points.Length % 3 != 0)
        {
            throw new ArgumentException("Points array length must be a multiple of 3.");
        }

        int n = points.Length / 3;
        if (n < 3)
        {
            throw new ArgumentException("At least 3 points are required to compute principal axes.");
        }

        // Compute mean
        double mx = 0, my = 0, mz = 0;
        for (int i = 0; i < n; i++)
        {
            mx += points[i * 3];
            my += points[i * 3 + 1];
            mz += points[i * 3 + 2];
        }
        mx /= n; my /= n; mz /= n;

        // Compute 3×3 covariance matrix (upper triangle, symmetric)
        double c00 = 0, c01 = 0, c02 = 0, c11 = 0, c12 = 0, c22 = 0;
        for (int i = 0; i < n; i++)
        {
            double dx = points[i * 3]     - mx;
            double dy = points[i * 3 + 1] - my;
            double dz = points[i * 3 + 2] - mz;
            c00 += dx * dx; c01 += dx * dy; c02 += dx * dz;
            c11 += dy * dy; c12 += dy * dz; c22 += dz * dz;
        }

        // Eigendecomposition via Jacobi iteration for 3×3 symmetric matrix
        var cov = new double[9]
        {
            c00, c01, c02,
            c01, c11, c12,
            c02, c12, c22,
        };

        Jacobi3x3(cov, out var eigenvalues, out var eigenvectors);

        // Sort by descending eigenvalue and return as row vectors
        var indices = new int[] { 0, 1, 2 };
        Array.Sort(indices, (a, b) => eigenvalues[b].CompareTo(eigenvalues[a]));

        var axes = new double[9];
        double[]? std = returnStd ? new double[3] : null;

        for (int r = 0; r < 3; r++)
        {
            int src = indices[r];
            axes[r * 3]     = eigenvectors[src * 3];
            axes[r * 3 + 1] = eigenvectors[src * 3 + 1];
            axes[r * 3 + 2] = eigenvectors[src * 3 + 2];

            if (std is not null)
            {
                std[r] = Math.Sqrt(Math.Abs(eigenvalues[src]) / n);
            }
        }

        // Ensure right-handed coordinate frame
        double det = Determinant3x3(axes);
        if (det < 0)
        {
            axes[6] = -axes[6];
            axes[7] = -axes[7];
            axes[8] = -axes[8];
        }

        return (axes, std);
    }

    /// <summary>
    /// Computes the centroid (arithmetic mean) of a set of 3D points.
    /// </summary>
    /// <param name="points">Flat array of point coordinates (Nx3 packed).</param>
    /// <returns>A 3-element array representing the centroid.</returns>
    /// <exception cref="ArgumentException">
    /// Thrown when the points array is empty or not a multiple of 3.
    /// </exception>
    public static double[] Centroid(double[] points)
    {
        ArgumentNullException.ThrowIfNull(points);

        if (points.Length == 0 || points.Length % 3 != 0)
        {
            throw new ArgumentException("Points array must be non-empty with length a multiple of 3.");
        }

        int n = points.Length / 3;
        double cx = 0, cy = 0, cz = 0;
        for (int i = 0; i < n; i++)
        {
            cx += points[i * 3];
            cy += points[i * 3 + 1];
            cz += points[i * 3 + 2];
        }

        return new[] { cx / n, cy / n, cz / n };
    }

    /// <summary>
    /// Computes the axis-aligned bounding box of a set of 3D points.
    /// </summary>
    /// <param name="points">Flat array of point coordinates (Nx3 packed).</param>
    /// <returns>
    /// A 6-element array <c>(xMin, xMax, yMin, yMax, zMin, zMax)</c>.
    /// </returns>
    /// <exception cref="ArgumentException">
    /// Thrown when the points array is empty or not a multiple of 3.
    /// </exception>
    public static double[] ComputeBounds(double[] points)
    {
        ArgumentNullException.ThrowIfNull(points);

        if (points.Length == 0 || points.Length % 3 != 0)
        {
            throw new ArgumentException("Points array must be non-empty with length a multiple of 3.");
        }

        double xMin = double.MaxValue, xMax = double.MinValue;
        double yMin = double.MaxValue, yMax = double.MinValue;
        double zMin = double.MaxValue, zMax = double.MinValue;

        int n = points.Length / 3;
        for (int i = 0; i < n; i++)
        {
            double x = points[i * 3], y = points[i * 3 + 1], z = points[i * 3 + 2];
            if (x < xMin) xMin = x; if (x > xMax) xMax = x;
            if (y < yMin) yMin = y; if (y > yMax) yMax = y;
            if (z < zMin) zMin = z; if (z > zMax) zMax = z;
        }

        return new[] { xMin, xMax, yMin, yMax, zMin, zMax };
    }

    /// <summary>
    /// Performs Jacobi eigenvalue iteration on a 3×3 symmetric matrix.
    /// </summary>
    /// <param name="matrix">Row-major 3×3 symmetric matrix (9 elements).</param>
    /// <param name="eigenvalues">Output: 3 eigenvalues.</param>
    /// <param name="eigenvectors">Output: 3 eigenvectors stored as rows in a 9-element array.</param>
    private static void Jacobi3x3(double[] matrix, out double[] eigenvalues, out double[] eigenvectors)
    {
        // Work on a copy
        var a = (double[])matrix.Clone();

        // Initialize eigenvectors to identity
        eigenvectors = new double[] { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

        const int maxIter = 50;
        for (int iter = 0; iter < maxIter; iter++)
        {
            // Find largest off-diagonal element
            double maxVal = 0;
            int p = 0, q = 1;
            for (int i = 0; i < 3; i++)
            {
                for (int j = i + 1; j < 3; j++)
                {
                    double v = Math.Abs(a[i * 3 + j]);
                    if (v > maxVal)
                    {
                        maxVal = v;
                        p = i;
                        q = j;
                    }
                }
            }

            if (maxVal < 1e-15) break;

            double app = a[p * 3 + p];
            double aqq = a[q * 3 + q];
            double apq = a[p * 3 + q];

            double theta = 0.5 * Math.Atan2(2.0 * apq, app - aqq);
            double c = Math.Cos(theta);
            double s = Math.Sin(theta);

            // Update matrix
            var aNew = (double[])a.Clone();
            aNew[p * 3 + p] = c * c * app + 2 * s * c * apq + s * s * aqq;
            aNew[q * 3 + q] = s * s * app - 2 * s * c * apq + c * c * aqq;
            aNew[p * 3 + q] = 0;
            aNew[q * 3 + p] = 0;

            for (int i = 0; i < 3; i++)
            {
                if (i == p || i == q) continue;
                double aip = a[i * 3 + p];
                double aiq = a[i * 3 + q];
                aNew[i * 3 + p] = c * aip + s * aiq;
                aNew[p * 3 + i] = aNew[i * 3 + p];
                aNew[i * 3 + q] = -s * aip + c * aiq;
                aNew[q * 3 + i] = aNew[i * 3 + q];
            }
            a = aNew;

            // Update eigenvectors
            var vNew = (double[])eigenvectors.Clone();
            for (int i = 0; i < 3; i++)
            {
                double vip = eigenvectors[p * 3 + i];
                double viq = eigenvectors[q * 3 + i];
                vNew[p * 3 + i] = c * vip + s * viq;
                vNew[q * 3 + i] = -s * vip + c * viq;
            }
            eigenvectors = vNew;
        }

        eigenvalues = new[] { a[0], a[4], a[8] };
    }

    /// <summary>
    /// Computes the determinant of a 3×3 row-major matrix.
    /// </summary>
    private static double Determinant3x3(double[] m)
    {
        return m[0] * (m[4] * m[8] - m[5] * m[7])
             - m[1] * (m[3] * m[8] - m[5] * m[6])
             + m[2] * (m[3] * m[7] - m[4] * m[6]);
    }
}
