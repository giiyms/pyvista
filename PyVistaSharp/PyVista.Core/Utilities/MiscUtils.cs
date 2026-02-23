using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;

namespace PyVista.Core.Utilities;

/// <summary>
/// Miscellaneous static utility methods used across PyVista.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.core.utilities.misc</c> module.
/// It provides validation helpers, bounds checking, vector checks, and other small
/// utilities that do not belong to a more specific category.
/// </para>
/// </summary>
public static class MiscUtils
{
    /// <summary>
    /// Returns <c>true</c> if the specified array contains any duplicate values.
    /// </summary>
    /// <typeparam name="T">The element type. Must implement <see cref="IEquatable{T}"/>.</typeparam>
    /// <param name="values">The collection to check.</param>
    /// <returns><c>true</c> when at least one duplicate exists; otherwise <c>false</c>.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="values"/> is <c>null</c>.
    /// </exception>
    public static bool HasDuplicates<T>(IEnumerable<T> values) where T : IEquatable<T>
    {
        ArgumentNullException.ThrowIfNull(values);

        var seen = new HashSet<T>();
        foreach (var item in values)
        {
            if (!seen.Add(item))
            {
                return true;
            }
        }

        return false;
    }

    /// <summary>
    /// Throws an <see cref="InvalidOperationException"/> if any value in the array
    /// is <see cref="double.NaN"/> or <see cref="double.PositiveInfinity"/>
    /// or <see cref="double.NegativeInfinity"/>.
    /// </summary>
    /// <param name="values">The array to validate.</param>
    /// <param name="name">Optional name of the array, used in the exception message.</param>
    /// <exception cref="InvalidOperationException">
    /// Thrown when a non-finite value is detected.
    /// </exception>
    public static void RaiseHasNonFinite(double[] values, string? name = null)
    {
        ArgumentNullException.ThrowIfNull(values);

        for (int i = 0; i < values.Length; i++)
        {
            if (!double.IsFinite(values[i]))
            {
                string arrayName = name ?? "Array";
                throw new InvalidOperationException(
                    $"{arrayName} contains non-finite values at index {i} (value: {values[i]}).");
            }
        }
    }

    /// <summary>
    /// Checks whether a 3D point lies within the given axis-aligned bounds.
    /// <para>
    /// Bounds are specified as <c>(xMin, xMax, yMin, yMax, zMin, zMax)</c>.
    /// This mirrors the Python <c>is_inside_bounds</c> utility.
    /// </para>
    /// </summary>
    /// <param name="point">A 3-component point <c>(x, y, z)</c>.</param>
    /// <param name="bounds">
    /// A 6-element array specifying <c>(xMin, xMax, yMin, yMax, zMin, zMax)</c>.
    /// </param>
    /// <returns><c>true</c> if the point lies inside the bounds (inclusive); otherwise <c>false</c>.</returns>
    /// <exception cref="ArgumentException">
    /// Thrown when the point or bounds array has an incorrect length.
    /// </exception>
    public static bool IsInsideBounds(double[] point, double[] bounds)
    {
        ArgumentNullException.ThrowIfNull(point);
        ArgumentNullException.ThrowIfNull(bounds);

        if (point.Length != 3)
        {
            throw new ArgumentException("Point must have exactly 3 components.", nameof(point));
        }

        if (bounds.Length != 6)
        {
            throw new ArgumentException(
                "Bounds must have exactly 6 elements (xMin, xMax, yMin, yMax, zMin, zMax).",
                nameof(bounds));
        }

        return point[0] >= bounds[0] && point[0] <= bounds[1]
            && point[1] >= bounds[2] && point[1] <= bounds[3]
            && point[2] >= bounds[4] && point[2] <= bounds[5];
    }

    /// <summary>
    /// N-dimensional bounds check. Tests whether each coordinate of the point falls
    /// within the corresponding pair in the bounds array.
    /// <para>
    /// The bounds array must have <c>2 × point.Length</c> elements, arranged as
    /// sequential <c>(min, max)</c> pairs.
    /// </para>
    /// </summary>
    /// <param name="point">The N-dimensional point.</param>
    /// <param name="bounds">
    /// A flat array of <c>(min, max)</c> pairs, with length <c>2 × point.Length</c>.
    /// </param>
    /// <returns><c>true</c> if all coordinates lie within their respective bounds.</returns>
    /// <exception cref="ArgumentException">
    /// Thrown when the bounds length does not match <c>2 × point.Length</c>.
    /// </exception>
    public static bool IsInsideBoundsNd(double[] point, double[] bounds)
    {
        ArgumentNullException.ThrowIfNull(point);
        ArgumentNullException.ThrowIfNull(bounds);

        if (bounds.Length != 2 * point.Length)
        {
            throw new ArgumentException("Bounds length must be 2 × point dimensions.");
        }

        for (int i = 0; i < point.Length; i++)
        {
            double lower = bounds[i * 2];
            double upper = bounds[i * 2 + 1];
            if (point[i] < lower || point[i] > upper)
            {
                return false;
            }
        }

        return true;
    }

    /// <summary>
    /// Asserts that a condition is <c>true</c>; otherwise, throws an exception.
    /// <para>
    /// This is a convenience wrapper analogous to Python's <c>assert</c> statement,
    /// intended for internal validation checks.
    /// </para>
    /// </summary>
    /// <param name="condition">The condition to test.</param>
    /// <param name="message">
    /// Optional message for the exception. When <c>null</c>, a default message is used.
    /// </param>
    /// <param name="callerName">
    /// Automatically populated with the name of the calling method.
    /// </param>
    /// <exception cref="InvalidOperationException">
    /// Thrown when <paramref name="condition"/> is <c>false</c>.
    /// </exception>
    public static void Assert(
        bool condition,
        string? message = null,
        [CallerMemberName] string callerName = "")
    {
        if (!condition)
        {
            throw new InvalidOperationException(
                message ?? $"Assertion failed in {callerName}.");
        }
    }

    /// <summary>
    /// Validates that a vector has exactly three components.
    /// </summary>
    /// <param name="vector">The vector to validate.</param>
    /// <param name="name">
    /// Optional name for the vector, used in exception messages. Defaults to <c>"Vector"</c>.
    /// </param>
    /// <exception cref="ArgumentException">
    /// Thrown when the vector does not have exactly 3 elements.
    /// </exception>
    public static void CheckValidVector(double[] vector, string? name = null)
    {
        ArgumentNullException.ThrowIfNull(vector);

        if (vector.Length != 3)
        {
            string label = string.IsNullOrEmpty(name) ? "Vector" : name;
            throw new ArgumentException($"{label} must be a length three array of floats.");
        }
    }

    /// <summary>
    /// Asserts that no unexpected keyword arguments remain.
    /// <para>
    /// In Python, this is used to verify that all keyword arguments have been consumed.
    /// In C#, this can be used to validate option dictionaries passed to flexible APIs.
    /// </para>
    /// </summary>
    /// <param name="kwargs">Dictionary of remaining keyword arguments.</param>
    /// <param name="callerName">
    /// Automatically populated with the calling method name.
    /// </param>
    /// <exception cref="ArgumentException">
    /// Thrown when the dictionary is not empty.
    /// </exception>
    public static void AssertEmptyKwargs(
        IDictionary<string, object>? kwargs,
        [CallerMemberName] string callerName = "")
    {
        if (kwargs is null || kwargs.Count == 0)
        {
            return;
        }

        var keys = string.Join(", ", kwargs.Keys.Select(k => $"\"{k}\""));
        string grammar = kwargs.Count == 1
            ? "is an invalid keyword argument"
            : "are invalid keyword arguments";
        throw new ArgumentException($"{keys} {grammar} for `{callerName}`.");
    }

    /// <summary>
    /// Checks that a numeric value falls within the specified inclusive range.
    /// </summary>
    /// <param name="value">The value to check.</param>
    /// <param name="min">The minimum acceptable value (inclusive).</param>
    /// <param name="max">The maximum acceptable value (inclusive).</param>
    /// <param name="parameterName">The parameter name, used in the exception message.</param>
    /// <exception cref="ArgumentOutOfRangeException">
    /// Thrown when <paramref name="value"/> is outside <c>[min, max]</c>.
    /// </exception>
    public static void CheckRange(double value, double min, double max, string parameterName)
    {
        if (value < min || value > max)
        {
            throw new ArgumentOutOfRangeException(
                parameterName,
                value,
                $"The value {value} for `{parameterName}` is outside the acceptable range [{min}, {max}].");
        }
    }

    /// <summary>
    /// Computes the element-wise reciprocal of an array, avoiding division by zero.
    /// <para>
    /// Elements whose absolute value is less than <paramref name="tolerance"/> are
    /// assigned <paramref name="valueIfZero"/> instead of their reciprocal.
    /// </para>
    /// </summary>
    /// <param name="values">Input array.</param>
    /// <param name="tolerance">
    /// Absolute value threshold below which an element is treated as zero. Defaults to <c>1e-8</c>.
    /// </param>
    /// <param name="valueIfZero">
    /// Value assigned to near-zero elements. Defaults to <c>0.0</c>.
    /// </param>
    /// <returns>A new array containing the element-wise reciprocals.</returns>
    public static double[] Reciprocal(double[] values, double tolerance = 1e-8, double valueIfZero = 0.0)
    {
        ArgumentNullException.ThrowIfNull(values);

        var result = new double[values.Length];
        for (int i = 0; i < values.Length; i++)
        {
            result[i] = Math.Abs(values[i]) < tolerance
                ? valueIfZero
                : 1.0 / values[i];
        }

        return result;
    }
}
