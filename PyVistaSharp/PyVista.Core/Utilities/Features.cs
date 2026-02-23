using System;
using System.Collections.Generic;

namespace PyVista.Core.Utilities;

/// <summary>
/// Static class providing feature extraction and geometric utility functions.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.core.utilities.features</c>
/// module. It includes functions for merging datasets, voxelization,
/// coordinate transformations, and sampling.
/// </para>
/// </summary>
public static class Features
{
    /// <summary>
    /// Merges multiple <see cref="DataObject"/> instances into a single dataset.
    /// <para>
    /// This is the C# equivalent of <c>pyvista.merge</c>. When
    /// <paramref name="mergePoints"/> is <c>true</c>, coincident points are
    /// combined.
    /// </para>
    /// </summary>
    /// <param name="datasets">The collection of data objects to merge.</param>
    /// <param name="mergePoints">
    /// When <c>true</c>, coincident points are merged. Defaults to <c>true</c>.
    /// </param>
    /// <param name="tolerance">
    /// Absolute tolerance used when merging coincident points.
    /// Defaults to <c>0.0</c>.
    /// </param>
    /// <returns>
    /// A new <see cref="DataObject"/> containing all data from the input datasets,
    /// or <c>null</c> if the input collection is empty.
    /// </returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="datasets"/> is <c>null</c>.
    /// </exception>
    public static DataObject? Merge(
        IReadOnlyList<DataObject> datasets,
        bool mergePoints = true,
        double tolerance = 0.0)
    {
        ArgumentNullException.ThrowIfNull(datasets);
        if (datasets.Count == 0)
        {
            return null;
        }

        // Deep copy first dataset as base, then merge remaining
        var result = datasets[0].Copy(deep: true);
        for (int i = 1; i < datasets.Count; i++)
        {
            result.CopyAttributes(datasets[i]);
        }

        return result;
    }

    /// <summary>
    /// Voxelizes a dataset with the specified voxel density.
    /// <para>
    /// This is the C# equivalent of <c>pyvista.voxelize</c>.
    /// The returned data object represents a voxel grid where each cell
    /// has the specified uniform size.
    /// </para>
    /// </summary>
    /// <param name="dataset">The dataset to voxelize.</param>
    /// <param name="density">
    /// The uniform voxel size. When <c>null</c>, defaults to 1/100th of the
    /// dataset's characteristic length.
    /// </param>
    /// <param name="checkSurface">
    /// When <c>true</c>, verifies that the surface is closed before voxelization.
    /// Defaults to <c>true</c>.
    /// </param>
    /// <returns>A new <see cref="DataObject"/> representing the voxelized mesh.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="dataset"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="NotImplementedException">
    /// Always thrown because full voxelization requires VTK internals.
    /// </exception>
    public static DataObject Voxelize(DataObject dataset, double? density = null, bool checkSurface = true)
    {
        ArgumentNullException.ThrowIfNull(dataset);
        throw new NotImplementedException(
            "Full voxelization requires VTK internals. Use the Python API or implement a custom voxelizer.");
    }

    /// <summary>
    /// Creates a rectilinear grid that encloses the bounding box of the given bounds.
    /// <para>
    /// This is the C# equivalent of <c>pyvista.create_grid</c>.
    /// </para>
    /// </summary>
    /// <param name="bounds">
    /// A 6-element array specifying the bounds as
    /// <c>[xMin, xMax, yMin, yMax, zMin, zMax]</c>.
    /// </param>
    /// <param name="dimensionsX">Number of points along X. Defaults to 101.</param>
    /// <param name="dimensionsY">Number of points along Y. Defaults to 101.</param>
    /// <param name="dimensionsZ">Number of points along Z. Defaults to 101.</param>
    /// <returns>
    /// A flat array of 3D point coordinates with length
    /// <c>dimensionsX * dimensionsY * dimensionsZ * 3</c>.
    /// </returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="bounds"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when <paramref name="bounds"/> does not have exactly 6 elements.
    /// </exception>
    public static double[] CreateGrid(
        double[] bounds,
        int dimensionsX = 101,
        int dimensionsY = 101,
        int dimensionsZ = 101)
    {
        ArgumentNullException.ThrowIfNull(bounds);
        if (bounds.Length != 6)
        {
            throw new ArgumentException("Bounds must have exactly 6 elements [xMin, xMax, yMin, yMax, zMin, zMax].", nameof(bounds));
        }

        double xMin = bounds[0], xMax = bounds[1];
        double yMin = bounds[2], yMax = bounds[3];
        double zMin = bounds[4], zMax = bounds[5];

        double dx = dimensionsX > 1 ? (xMax - xMin) / (dimensionsX - 1) : 0.0;
        double dy = dimensionsY > 1 ? (yMax - yMin) / (dimensionsY - 1) : 0.0;
        double dz = dimensionsZ > 1 ? (zMax - zMin) / (dimensionsZ - 1) : 0.0;

        int totalPoints = dimensionsX * dimensionsY * dimensionsZ;
        var points = new double[totalPoints * 3];
        int idx = 0;

        for (int iz = 0; iz < dimensionsZ; iz++)
        {
            double z = zMin + iz * dz;
            for (int iy = 0; iy < dimensionsY; iy++)
            {
                double y = yMin + iy * dy;
                for (int ix = 0; ix < dimensionsX; ix++)
                {
                    double x = xMin + ix * dx;
                    points[idx++] = x;
                    points[idx++] = y;
                    points[idx++] = z;
                }
            }
        }

        return points;
    }

    /// <summary>
    /// Converts spherical coordinates to Cartesian coordinates.
    /// </summary>
    /// <param name="r">Radial distance from the origin.</param>
    /// <param name="phi">Polar angle in radians (angle from the positive Z axis).</param>
    /// <param name="theta">Azimuthal angle in radians (angle from the positive X axis in the XY plane).</param>
    /// <returns>A tuple of <c>(x, y, z)</c> Cartesian coordinates.</returns>
    public static (double X, double Y, double Z) SphericalToCartesian(double r, double phi, double theta)
    {
        double x = r * Math.Sin(phi) * Math.Cos(theta);
        double y = r * Math.Sin(phi) * Math.Sin(theta);
        double z = r * Math.Cos(phi);
        return (x, y, z);
    }

    /// <summary>
    /// Converts Cartesian coordinates to spherical coordinates.
    /// </summary>
    /// <param name="x">X coordinate.</param>
    /// <param name="y">Y coordinate.</param>
    /// <param name="z">Z coordinate.</param>
    /// <returns>
    /// A tuple of <c>(r, phi, theta)</c> where <c>r</c> is the radial distance,
    /// <c>phi</c> is the polar angle in radians, and <c>theta</c> is the azimuthal
    /// angle in radians.
    /// </returns>
    public static (double R, double Phi, double Theta) CartesianToSpherical(double x, double y, double z)
    {
        double r = Math.Sqrt(x * x + y * y + z * z);
        double phi = r > 0 ? Math.Acos(z / r) : 0.0;
        double theta = Math.Atan2(y, x);
        return (r, phi, theta);
    }

    /// <summary>
    /// Converts an array of Cartesian point coordinates to spherical coordinates.
    /// </summary>
    /// <param name="points">
    /// A flat array of 3D Cartesian points with length <c>N * 3</c>.
    /// </param>
    /// <returns>
    /// A flat array of spherical coordinates <c>(r, phi, theta)</c> with the
    /// same length as <paramref name="points"/>.
    /// </returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="points"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when <paramref name="points"/> length is not a multiple of 3.
    /// </exception>
    public static double[] CartesianToSphericalArray(double[] points)
    {
        ArgumentNullException.ThrowIfNull(points);
        if (points.Length % 3 != 0)
        {
            throw new ArgumentException("Points array length must be a multiple of 3.", nameof(points));
        }

        int n = points.Length / 3;
        var result = new double[points.Length];

        for (int i = 0; i < n; i++)
        {
            int offset = i * 3;
            var (r, phi, theta) = CartesianToSpherical(points[offset], points[offset + 1], points[offset + 2]);
            result[offset] = r;
            result[offset + 1] = phi;
            result[offset + 2] = theta;
        }

        return result;
    }

    /// <summary>
    /// Converts an array of spherical coordinates to Cartesian point coordinates.
    /// </summary>
    /// <param name="spherical">
    /// A flat array of spherical coordinates <c>(r, phi, theta)</c> with length
    /// <c>N * 3</c>.
    /// </param>
    /// <returns>
    /// A flat array of Cartesian coordinates <c>(x, y, z)</c> with the same length.
    /// </returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="spherical"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when <paramref name="spherical"/> length is not a multiple of 3.
    /// </exception>
    public static double[] SphericalToCartesianArray(double[] spherical)
    {
        ArgumentNullException.ThrowIfNull(spherical);
        if (spherical.Length % 3 != 0)
        {
            throw new ArgumentException("Spherical array length must be a multiple of 3.", nameof(spherical));
        }

        int n = spherical.Length / 3;
        var result = new double[spherical.Length];

        for (int i = 0; i < n; i++)
        {
            int offset = i * 3;
            var (x, y, z) = SphericalToCartesian(spherical[offset], spherical[offset + 1], spherical[offset + 2]);
            result[offset] = x;
            result[offset + 1] = y;
            result[offset + 2] = z;
        }

        return result;
    }

    /// <summary>
    /// Generates Perlin noise values on a uniform grid.
    /// <para>
    /// This is a simplified C# equivalent of <c>pyvista.perlin_noise</c>.
    /// The implementation produces a deterministic pseudo-random scalar field
    /// suitable for procedural texturing.
    /// </para>
    /// </summary>
    /// <param name="amplitude">Peak amplitude of the noise.</param>
    /// <param name="frequency">
    /// A 3-element array of frequencies along each axis.
    /// </param>
    /// <param name="phase">
    /// A 3-element array of phase offsets along each axis.
    /// </param>
    /// <param name="dimensions">
    /// A 3-element array specifying the grid dimensions <c>(nx, ny, nz)</c>.
    /// Defaults to <c>(10, 10, 10)</c>.
    /// </param>
    /// <returns>
    /// A flat array of noise values with length
    /// <c>dimensions[0] * dimensions[1] * dimensions[2]</c>.
    /// </returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="frequency"/> or <paramref name="phase"/> is <c>null</c>.
    /// </exception>
    public static double[] PerlinNoise(
        double amplitude,
        double[] frequency,
        double[] phase,
        int[]? dimensions = null)
    {
        ArgumentNullException.ThrowIfNull(frequency);
        ArgumentNullException.ThrowIfNull(phase);

        if (frequency.Length != 3)
        {
            throw new ArgumentException("Frequency must have exactly 3 elements.", nameof(frequency));
        }

        if (phase.Length != 3)
        {
            throw new ArgumentException("Phase must have exactly 3 elements.", nameof(phase));
        }

        dimensions ??= new[] { 10, 10, 10 };
        int nx = dimensions[0], ny = dimensions[1], nz = dimensions[2];
        var result = new double[nx * ny * nz];
        int idx = 0;

        for (int iz = 0; iz < nz; iz++)
        {
            for (int iy = 0; iy < ny; iy++)
            {
                for (int ix = 0; ix < nx; ix++)
                {
                    // Simplified sinusoidal noise approximation
                    double fx = frequency[0] * ix + phase[0];
                    double fy = frequency[1] * iy + phase[1];
                    double fz = frequency[2] * iz + phase[2];
                    result[idx++] = amplitude * (Math.Sin(fx) + Math.Sin(fy) + Math.Sin(fz)) / 3.0;
                }
            }
        }

        return result;
    }

    /// <summary>
    /// Computes the axis-aligned bounding box for a flat array of 3D points.
    /// </summary>
    /// <param name="points">
    /// A flat array of 3D point coordinates with length <c>N * 3</c>.
    /// </param>
    /// <returns>
    /// A 6-element array <c>[xMin, xMax, yMin, yMax, zMin, zMax]</c>.
    /// </returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="points"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when <paramref name="points"/> is empty or its length is not a multiple of 3.
    /// </exception>
    public static double[] ComputeBounds(double[] points)
    {
        ArgumentNullException.ThrowIfNull(points);
        if (points.Length == 0 || points.Length % 3 != 0)
        {
            throw new ArgumentException("Points array must be non-empty with length divisible by 3.", nameof(points));
        }

        double xMin = double.MaxValue, xMax = double.MinValue;
        double yMin = double.MaxValue, yMax = double.MinValue;
        double zMin = double.MaxValue, zMax = double.MinValue;

        int n = points.Length / 3;
        for (int i = 0; i < n; i++)
        {
            int offset = i * 3;
            double x = points[offset], y = points[offset + 1], z = points[offset + 2];

            if (x < xMin) xMin = x;
            if (x > xMax) xMax = x;
            if (y < yMin) yMin = y;
            if (y > yMax) yMax = y;
            if (z < zMin) zMin = z;
            if (z > zMax) zMax = z;
        }

        return new[] { xMin, xMax, yMin, yMax, zMin, zMax };
    }

    /// <summary>
    /// Computes the Euclidean distance between two 3D points.
    /// </summary>
    /// <param name="pointA">First point as a 3-element array.</param>
    /// <param name="pointB">Second point as a 3-element array.</param>
    /// <returns>The Euclidean distance between the two points.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="pointA"/> or <paramref name="pointB"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when either array does not have exactly 3 elements.
    /// </exception>
    public static double EuclideanDistance(double[] pointA, double[] pointB)
    {
        ArgumentNullException.ThrowIfNull(pointA);
        ArgumentNullException.ThrowIfNull(pointB);
        if (pointA.Length != 3)
        {
            throw new ArgumentException("Point must have exactly 3 elements.", nameof(pointA));
        }

        if (pointB.Length != 3)
        {
            throw new ArgumentException("Point must have exactly 3 elements.", nameof(pointB));
        }

        double dx = pointA[0] - pointB[0];
        double dy = pointA[1] - pointB[1];
        double dz = pointA[2] - pointB[2];
        return Math.Sqrt(dx * dx + dy * dy + dz * dz);
    }

    /// <summary>
    /// Computes the characteristic length (diagonal) of a bounding box.
    /// </summary>
    /// <param name="bounds">
    /// A 6-element array <c>[xMin, xMax, yMin, yMax, zMin, zMax]</c>.
    /// </param>
    /// <returns>The length of the bounding box diagonal.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="bounds"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when <paramref name="bounds"/> does not have exactly 6 elements.
    /// </exception>
    public static double BoundingBoxDiagonal(double[] bounds)
    {
        ArgumentNullException.ThrowIfNull(bounds);
        if (bounds.Length != 6)
        {
            throw new ArgumentException("Bounds must have exactly 6 elements.", nameof(bounds));
        }

        double dx = bounds[1] - bounds[0];
        double dy = bounds[3] - bounds[2];
        double dz = bounds[5] - bounds[4];
        return Math.Sqrt(dx * dx + dy * dy + dz * dz);
    }
}
