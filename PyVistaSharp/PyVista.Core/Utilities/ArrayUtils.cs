using System;
using System.Collections.Generic;
using System.Linq;

namespace PyVista.Core.Utilities;

/// <summary>
/// Static utility methods for array operations.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.core.utilities.arrays</c> module.
/// It provides helper methods for converting, searching, and validating arrays used
/// throughout PyVista data structures.
/// </para>
/// </summary>
public static class ArrayUtils
{
    /// <summary>
    /// Converts a numeric array between common representations.
    /// <para>
    /// When the input is a <see cref="float"/> array it is returned as-is.
    /// When the input is a <see cref="double"/> array it is narrowed to <see cref="float"/>.
    /// Integer arrays are widened to <see cref="double"/>.
    /// </para>
    /// </summary>
    /// <param name="array">The source array to convert.</param>
    /// <param name="name">Optional name to associate with the converted array.</param>
    /// <param name="deep">
    /// When <c>true</c>, the returned array is always a new copy.
    /// When <c>false</c>, the original reference may be returned if no conversion is needed.
    /// </param>
    /// <returns>The converted array as <see cref="double"/>[], or <c>null</c> if the input is <c>null</c>.</returns>
    public static double[]? ConvertArray(Array? array, string? name = null, bool deep = false)
    {
        if (array is null)
        {
            return null;
        }

        double[] result;

        if (array is double[] doubleArr)
        {
            result = deep ? (double[])doubleArr.Clone() : doubleArr;
        }
        else if (array is float[] floatArr)
        {
            result = new double[floatArr.Length];
            for (int i = 0; i < floatArr.Length; i++)
            {
                result[i] = floatArr[i];
            }
        }
        else if (array is int[] intArr)
        {
            result = new double[intArr.Length];
            for (int i = 0; i < intArr.Length; i++)
            {
                result[i] = intArr[i];
            }
        }
        else if (array is long[] longArr)
        {
            result = new double[longArr.Length];
            for (int i = 0; i < longArr.Length; i++)
            {
                result[i] = longArr[i];
            }
        }
        else
        {
            throw new ArgumentException($"Unsupported array element type: {array.GetType().GetElementType()}.");
        }

        return result;
    }

    /// <summary>
    /// Searches point data, cell data, and field data for an array with the given name.
    /// <para>
    /// This mirrors the Python <c>get_array</c> helper which searches across all
    /// association types and returns the first match, preferring the specified association.
    /// </para>
    /// </summary>
    /// <param name="arrays">
    /// A dictionary mapping array names to their values, representing all arrays in a dataset.
    /// </param>
    /// <param name="name">The name of the array to retrieve.</param>
    /// <param name="throwOnMissing">
    /// When <c>true</c>, a <see cref="KeyNotFoundException"/> is thrown if the array is not found.
    /// </param>
    /// <returns>The matching array, or <c>null</c> when not found and <paramref name="throwOnMissing"/> is <c>false</c>.</returns>
    /// <exception cref="KeyNotFoundException">
    /// Thrown when <paramref name="throwOnMissing"/> is <c>true</c> and no matching array exists.
    /// </exception>
    public static double[]? GetArray(
        IDictionary<string, double[]> arrays,
        string name,
        bool throwOnMissing = false)
    {
        ArgumentNullException.ThrowIfNull(arrays);
        ArgumentNullException.ThrowIfNull(name);

        if (arrays.TryGetValue(name, out var result))
        {
            return result;
        }

        if (throwOnMissing)
        {
            throw new KeyNotFoundException($"Data array ({name}) not present in this dataset.");
        }

        return null;
    }

    /// <summary>
    /// Sets default active scalar, vector, and normal arrays on a dataset's data dictionaries.
    /// <para>
    /// When a dataset is first loaded or constructed, this method ensures that the active
    /// arrays are initialized to the first available array in each category.
    /// </para>
    /// </summary>
    /// <param name="pointData">Dictionary of point data arrays.</param>
    /// <param name="cellData">Dictionary of cell data arrays.</param>
    /// <param name="activeScalars">
    /// Receives the name of the active scalars array, or <c>null</c> if none are available.
    /// </param>
    public static void SetDefaultActiveArrays(
        IDictionary<string, double[]> pointData,
        IDictionary<string, double[]> cellData,
        out string? activeScalars)
    {
        ArgumentNullException.ThrowIfNull(pointData);
        ArgumentNullException.ThrowIfNull(cellData);

        activeScalars = null;

        if (pointData.Count > 0)
        {
            activeScalars = pointData.Keys.First();
        }
        else if (cellData.Count > 0)
        {
            activeScalars = cellData.Keys.First();
        }
    }

    /// <summary>
    /// Raises a <see cref="InvalidOperationException"/> when the number of scalars does not
    /// match the number of points or cells in a dataset.
    /// </summary>
    /// <param name="scalarsLength">The length of the scalars array.</param>
    /// <param name="numberOfPoints">The number of points in the dataset.</param>
    /// <param name="numberOfCells">The number of cells in the dataset.</param>
    /// <exception cref="InvalidOperationException">
    /// Always thrown, describing the size mismatch.
    /// </exception>
    public static void RaiseNotMatching(int scalarsLength, int numberOfPoints, int numberOfCells)
    {
        throw new InvalidOperationException(
            $"Number of scalars ({scalarsLength}) must match either the number of points " +
            $"({numberOfPoints}) or the number of cells ({numberOfCells}).");
    }

    /// <summary>
    /// Returns <c>true</c> if the given array contains any duplicate values.
    /// </summary>
    /// <param name="array">The array to check.</param>
    /// <returns><c>true</c> if duplicates exist; otherwise <c>false</c>.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="array"/> is <c>null</c>.
    /// </exception>
    public static bool HasDuplicates(double[] array)
    {
        ArgumentNullException.ThrowIfNull(array);
        var set = new HashSet<double>(array.Length);
        foreach (var value in array)
        {
            if (!set.Add(value))
            {
                return true;
            }
        }

        return false;
    }

    /// <summary>
    /// Returns <c>true</c> if the given integer array contains any duplicate values.
    /// </summary>
    /// <param name="array">The array to check.</param>
    /// <returns><c>true</c> if duplicates exist; otherwise <c>false</c>.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="array"/> is <c>null</c>.
    /// </exception>
    public static bool HasDuplicates(int[] array)
    {
        ArgumentNullException.ThrowIfNull(array);
        var set = new HashSet<int>(array.Length);
        foreach (var value in array)
        {
            if (!set.Add(value))
            {
                return true;
            }
        }

        return false;
    }

    /// <summary>
    /// Throws a <see cref="InvalidOperationException"/> if the array contains duplicate values.
    /// </summary>
    /// <param name="array">The array to validate.</param>
    /// <exception cref="InvalidOperationException">
    /// Thrown when the array contains one or more duplicate values.
    /// </exception>
    public static void RaiseHasDuplicates(double[] array)
    {
        if (HasDuplicates(array))
        {
            throw new InvalidOperationException("Array contains duplicate values.");
        }
    }

    /// <summary>
    /// Parses a field association string into a <see cref="FieldAssociation"/> value.
    /// <para>
    /// Accepts values such as <c>"point"</c>, <c>"cell"</c>, <c>"field"</c>, and <c>"row"</c>
    /// (case-insensitive), as well as their single-letter abbreviations.
    /// </para>
    /// </summary>
    /// <param name="field">The string representation of the field association.</param>
    /// <returns>The corresponding <see cref="FieldAssociation"/> enumeration value.</returns>
    /// <exception cref="ArgumentException">
    /// Thrown when <paramref name="field"/> does not map to a known association.
    /// </exception>
    public static FieldAssociation ParseFieldChoice(string field)
    {
        ArgumentNullException.ThrowIfNull(field);

        return field.Trim().ToLowerInvariant() switch
        {
            "cell" or "c" or "cells" => FieldAssociation.Cell,
            "point" or "p" or "points" => FieldAssociation.Point,
            "field" or "f" or "fields" => FieldAssociation.None,
            "row" or "r" => FieldAssociation.Row,
            _ => throw new ArgumentException($"Data field ({field}) not supported."),
        };
    }

    /// <summary>
    /// Coerces a point-like input into an Nx3 array representation.
    /// <para>
    /// If the input is a single 3-element point, it is wrapped into a 1×3 array.
    /// </para>
    /// </summary>
    /// <param name="points">
    /// A flat array of coordinates. Must have a length that is a multiple of 3.
    /// </param>
    /// <param name="isSingular">
    /// Set to <c>true</c> when the input represents a single point (length 3).
    /// </param>
    /// <returns>
    /// A jagged array where each element is a 3-component point.
    /// </returns>
    /// <exception cref="ArgumentException">
    /// Thrown when the array length is not a multiple of 3.
    /// </exception>
    public static double[][] CoercePointsLike(double[] points, out bool isSingular)
    {
        ArgumentNullException.ThrowIfNull(points);

        if (points.Length % 3 != 0)
        {
            throw new ArgumentException("Array of points must have a length that is a multiple of 3.");
        }

        int count = points.Length / 3;
        isSingular = count == 1;

        var result = new double[count][];
        for (int i = 0; i < count; i++)
        {
            result[i] = new double[] { points[i * 3], points[i * 3 + 1], points[i * 3 + 2] };
        }

        return result;
    }

    /// <summary>
    /// Returns the association of a named array within the given data dictionaries.
    /// </summary>
    /// <param name="pointData">Dictionary of point data arrays.</param>
    /// <param name="cellData">Dictionary of cell data arrays.</param>
    /// <param name="fieldData">Dictionary of field data arrays.</param>
    /// <param name="name">Name of the array to look up.</param>
    /// <param name="preference">Preferred association when the name exists in multiple dictionaries.</param>
    /// <returns>The <see cref="FieldAssociation"/> where the array was found.</returns>
    /// <exception cref="KeyNotFoundException">
    /// Thrown when no array with the given name exists in any of the dictionaries.
    /// </exception>
    public static FieldAssociation GetArrayAssociation(
        IDictionary<string, double[]> pointData,
        IDictionary<string, double[]> cellData,
        IDictionary<string, double[]> fieldData,
        string name,
        FieldAssociation preference = FieldAssociation.Cell)
    {
        ArgumentNullException.ThrowIfNull(name);

        bool inPoint = pointData.ContainsKey(name);
        bool inCell = cellData.ContainsKey(name);
        bool inField = fieldData.ContainsKey(name);

        if (!inPoint && !inCell && !inField)
        {
            throw new KeyNotFoundException($"Data array ({name}) not present in this dataset.");
        }

        if (preference == FieldAssociation.Point && inPoint) return FieldAssociation.Point;
        if (preference == FieldAssociation.Cell && inCell) return FieldAssociation.Cell;
        if (preference == FieldAssociation.None && inField) return FieldAssociation.None;

        if (inPoint) return FieldAssociation.Point;
        if (inCell) return FieldAssociation.Cell;
        return FieldAssociation.None;
    }

    /// <summary>
    /// Converts a string array into a <see cref="double"/> array by parsing each element.
    /// </summary>
    /// <param name="values">Array of string values to convert.</param>
    /// <returns>An array of parsed <see cref="double"/> values.</returns>
    /// <exception cref="FormatException">
    /// Thrown when any element cannot be parsed as a <see cref="double"/>.
    /// </exception>
    public static double[] ConvertStringArray(string[] values)
    {
        ArgumentNullException.ThrowIfNull(values);

        var result = new double[values.Length];
        for (int i = 0; i < values.Length; i++)
        {
            result[i] = double.Parse(values[i], System.Globalization.CultureInfo.InvariantCulture);
        }

        return result;
    }
}
