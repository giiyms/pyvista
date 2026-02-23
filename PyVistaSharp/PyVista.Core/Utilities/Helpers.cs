using System;
using System.Collections.Generic;
using System.Linq;

namespace PyVista.Core.Utilities;

/// <summary>
/// Static helper methods for general PyVista operations.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.core.utilities.helpers</c> module.
/// It provides wrapping, type-checking, point-cloud generation, and geometric helpers
/// such as plane/line fitting and axis rotation.
/// </para>
/// </summary>
public static class Helpers
{
    /// <summary>
    /// Predefined normal vectors indexed by axis name.
    /// </summary>
    private static readonly Dictionary<string, double[]> Normals = new(StringComparer.OrdinalIgnoreCase)
    {
        ["x"]  = new[] { 1.0,  0.0,  0.0 },
        ["y"]  = new[] { 0.0,  1.0,  0.0 },
        ["z"]  = new[] { 0.0,  0.0,  1.0 },
        ["-x"] = new[] { -1.0, 0.0,  0.0 },
        ["-y"] = new[] { 0.0, -1.0,  0.0 },
        ["-z"] = new[] { 0.0,  0.0, -1.0 },
    };

    /// <summary>
    /// Wraps a raw data object into its corresponding PyVista type.
    /// <para>
    /// If the input is already a <see cref="DataObject"/>, it is returned unchanged.
    /// If the input is a 2D point array (Nx3), it is wrapped as a point cloud.
    /// Otherwise, an exception is thrown.
    /// </para>
    /// </summary>
    /// <param name="dataset">
    /// The object to wrap. May be a <see cref="DataObject"/>, a jagged <c>double[][]</c>
    /// representing Nx3 points, or <c>null</c>.
    /// </param>
    /// <returns>
    /// The wrapped <see cref="DataObject"/>, or <c>null</c> when <paramref name="dataset"/> is <c>null</c>.
    /// </returns>
    /// <exception cref="NotSupportedException">
    /// Thrown when the input type cannot be wrapped.
    /// </exception>
    public static DataObject? WrapDataObject(object? dataset)
    {
        if (dataset is null)
        {
            return null;
        }

        if (dataset is DataObject existing)
        {
            return existing;
        }

        if (dataset is double[][] points)
        {
            return GeneratePointCloud(points);
        }

        throw new NotSupportedException(
            $"Unable to wrap object of type {dataset.GetType()} into a PyVista type.");
    }

    /// <summary>
    /// Returns <c>true</c> if the given object is a PyVista <see cref="DataObject"/>.
    /// </summary>
    /// <param name="obj">The object to test.</param>
    /// <returns><c>true</c> when <paramref name="obj"/> is a <see cref="DataObject"/>.</returns>
    public static bool IsPyVistaDataSet(object? obj)
    {
        return obj is DataObject;
    }

    /// <summary>
    /// Creates a simple point cloud <see cref="DataObject"/> from an Nx3 jagged array of points.
    /// <para>
    /// Each inner array must have exactly 3 elements representing (x, y, z).
    /// The returned object stores the points as flat field data under the key <c>"Points"</c>.
    /// </para>
    /// </summary>
    /// <param name="points">Nx3 jagged array of point coordinates.</param>
    /// <returns>A <see cref="PointCloudDataObject"/> containing the point cloud.</returns>
    /// <exception cref="ArgumentException">
    /// Thrown when any inner array does not have exactly 3 elements.
    /// </exception>
    public static PointCloudDataObject GeneratePointCloud(double[][] points)
    {
        ArgumentNullException.ThrowIfNull(points);

        var flat = new double[points.Length * 3];
        for (int i = 0; i < points.Length; i++)
        {
            if (points[i].Length != 3)
            {
                throw new ArgumentException(
                    $"Point at index {i} has {points[i].Length} components; expected 3.");
            }

            flat[i * 3]     = points[i][0];
            flat[i * 3 + 1] = points[i][1];
            flat[i * 3 + 2] = points[i][2];
        }

        var cloud = new PointCloudDataObject();
        cloud.AddFieldData(flat, "Points");
        return cloud;
    }

    /// <summary>
    /// Resolves a normal vector from either a 3-element array or an axis name string.
    /// <para>
    /// Accepted axis names are <c>"x"</c>, <c>"y"</c>, <c>"z"</c>, <c>"-x"</c>,
    /// <c>"-y"</c>, and <c>"-z"</c> (case-insensitive).
    /// </para>
    /// </summary>
    /// <param name="normalOrAxis">
    /// Either a <c>double[3]</c> normal vector or a <see cref="string"/> axis name.
    /// </param>
    /// <returns>A 3-component unit normal vector.</returns>
    /// <exception cref="ArgumentException">
    /// Thrown when the input is neither a valid axis name nor a 3-element array.
    /// </exception>
    public static double[] ResolveNormal(object normalOrAxis)
    {
        ArgumentNullException.ThrowIfNull(normalOrAxis);

        if (normalOrAxis is string axisName)
        {
            if (!Normals.TryGetValue(axisName.Trim(), out var n))
            {
                throw new ArgumentException(
                    $"Unknown axis name '{axisName}'. Must be one of: x, y, z, -x, -y, -z.");
            }

            return (double[])n.Clone();
        }

        if (normalOrAxis is double[] vec)
        {
            MiscUtils.CheckValidVector(vec, "normal");
            double mag = Math.Sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
            if (mag < 1e-12)
            {
                throw new ArgumentException("Normal vector must have non-zero magnitude.");
            }

            return new[] { vec[0] / mag, vec[1] / mag, vec[2] / mag };
        }

        throw new ArgumentException(
            $"Expected a double[3] vector or axis name string, got {normalOrAxis.GetType()}.");
    }

    /// <summary>
    /// Rotates an array of 3D points about a specified axis by the given angle.
    /// </summary>
    /// <param name="points">
    /// Flat array of point coordinates with length that is a multiple of 3.
    /// </param>
    /// <param name="angleDegrees">Rotation angle in degrees.</param>
    /// <param name="axis">
    /// Axis of rotation. Must be <c>"x"</c>, <c>"y"</c>, or <c>"z"</c>.
    /// </param>
    /// <returns>A new flat array of rotated point coordinates.</returns>
    /// <exception cref="ArgumentException">
    /// Thrown when the axis is invalid or the points array length is not a multiple of 3.
    /// </exception>
    public static double[] AxisRotation(double[] points, double angleDegrees, string axis = "z")
    {
        ArgumentNullException.ThrowIfNull(points);
        ArgumentNullException.ThrowIfNull(axis);

        if (points.Length % 3 != 0)
        {
            throw new ArgumentException("Points array length must be a multiple of 3.");
        }

        double rad = angleDegrees * Math.PI / 180.0;
        double cos = Math.Cos(rad);
        double sin = Math.Sin(rad);

        // Build 3×3 rotation matrix
        double[] m = axis.ToLowerInvariant() switch
        {
            "x" => new[] { 1, 0, 0, 0, cos, -sin, 0, sin, cos },
            "y" => new[] { cos, 0, sin, 0, 1, 0, -sin, 0, cos },
            "z" => new[] { cos, -sin, 0, sin, cos, 0, 0, 0, 1 },
            _ => throw new ArgumentException(
                $"Invalid axis '{axis}'. Must be \"x\", \"y\", or \"z\"."),
        };

        int nPoints = points.Length / 3;
        var result = new double[points.Length];

        for (int i = 0; i < nPoints; i++)
        {
            int off = i * 3;
            double x = points[off], y = points[off + 1], z = points[off + 2];
            result[off]     = m[0] * x + m[1] * y + m[2] * z;
            result[off + 1] = m[3] * x + m[4] * y + m[5] * z;
            result[off + 2] = m[6] * x + m[7] * y + m[8] * z;
        }

        return result;
    }

    /// <summary>
    /// Checks whether a 3D point lies within the specified axis-aligned bounds.
    /// <para>
    /// This is a convenience wrapper around <see cref="MiscUtils.IsInsideBounds"/>.
    /// </para>
    /// </summary>
    /// <param name="point">A 3-component point <c>(x, y, z)</c>.</param>
    /// <param name="bounds">
    /// A 6-element array <c>(xMin, xMax, yMin, yMax, zMin, zMax)</c>.
    /// </param>
    /// <returns><c>true</c> if the point is inside the bounds (inclusive).</returns>
    public static bool IsInsideBounds(double[] point, double[] bounds)
    {
        return MiscUtils.IsInsideBounds(point, bounds);
    }

    /// <summary>
    /// Normalizes a 3D vector to unit length.
    /// </summary>
    /// <param name="vector">A 3-component vector.</param>
    /// <returns>A new array containing the normalized vector.</returns>
    /// <exception cref="ArgumentException">
    /// Thrown when the vector has zero magnitude.
    /// </exception>
    public static double[] NormalizeVector(double[] vector)
    {
        MiscUtils.CheckValidVector(vector, "vector");

        double mag = Math.Sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
        if (mag < 1e-12)
        {
            throw new ArgumentException("Cannot normalize a zero-length vector.");
        }

        return new[] { vector[0] / mag, vector[1] / mag, vector[2] / mag };
    }

    /// <summary>
    /// Computes the Euclidean distance between two 3D points.
    /// </summary>
    /// <param name="a">First point as <c>(x, y, z)</c>.</param>
    /// <param name="b">Second point as <c>(x, y, z)</c>.</param>
    /// <returns>The Euclidean distance.</returns>
    public static double Distance(double[] a, double[] b)
    {
        MiscUtils.CheckValidVector(a, "a");
        MiscUtils.CheckValidVector(b, "b");

        double dx = a[0] - b[0];
        double dy = a[1] - b[1];
        double dz = a[2] - b[2];
        return Math.Sqrt(dx * dx + dy * dy + dz * dz);
    }
}

/// <summary>
/// A minimal concrete <see cref="DataObject"/> used to represent a point cloud.
/// <para>
/// Points are stored as a flat <c>double[]</c> in the field data under the key <c>"Points"</c>.
/// </para>
/// </summary>
public sealed class PointCloudDataObject : DataObject
{
    /// <summary>
    /// Gets a value indicating whether this point cloud contains no points.
    /// </summary>
    public override bool IsEmpty
    {
        get
        {
            if (!FieldData.ContainsKey("Points")) return true;
            return FieldData.GetArray("Points").Length == 0;
        }
    }

    /// <summary>
    /// Gets the number of points in this point cloud.
    /// </summary>
    public int NumberOfPoints
    {
        get
        {
            if (!FieldData.ContainsKey("Points")) return 0;
            return FieldData.GetArray("Points").Length / 3;
        }
    }

    /// <inheritdoc />
    public override (double Min, double Max) GetDataRange(
        string? name = null,
        FieldAssociation preference = FieldAssociation.Point)
    {
        string key = name ?? "Points";
        if (!FieldData.ContainsKey(key))
        {
            return (double.NaN, double.NaN);
        }

        var arr = FieldData.GetArray(key);
        if (arr.Length == 0)
        {
            return (double.NaN, double.NaN);
        }

        double min = double.MaxValue;
        double max = double.MinValue;
        foreach (double v in arr)
        {
            if (double.IsNaN(v)) continue;
            if (v < min) min = v;
            if (v > max) max = v;
        }

        return (min, max);
    }
}
