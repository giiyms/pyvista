using System;
using System.Collections.Generic;
using System.Linq;

namespace PyVista.Core.Utilities;

/// <summary>
/// Static utility methods for cell array operations.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.core.utilities.cells</c> module.
/// It provides helpers for creating and interrogating VTK-style cell connectivity arrays
/// without requiring a direct VTK dependency.
/// </para>
/// </summary>
public static class CellUtils
{
    /// <summary>
    /// Gets the number of cells encoded in a VTK cell connectivity array.
    /// <para>
    /// A VTK cell connectivity array stores cells sequentially where each cell
    /// is preceded by a count of the number of points that define it.
    /// For example, two triangles: <c>[3, 0, 1, 2, 3, 3, 4, 5]</c>.
    /// </para>
    /// </summary>
    /// <param name="cells">
    /// A flat array of VTK cell connectivity data.
    /// </param>
    /// <returns>The total number of cells encoded in the array.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="cells"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when the connectivity array is malformed.
    /// </exception>
    public static int NumCellsFromCells(int[] cells)
    {
        ArgumentNullException.ThrowIfNull(cells);

        int index = 0;
        int count = 0;

        while (index < cells.Length)
        {
            int nPoints = cells[index];
            if (nPoints < 0 || index + nPoints >= cells.Length)
            {
                throw new ArgumentException(
                    $"Malformed cell connectivity array at index {index}: cell size {nPoints} is invalid.");
            }

            index += nPoints + 1;
            count++;
        }

        return count;
    }

    /// <summary>
    /// Gets the number of cells from a VTK offset array.
    /// <para>
    /// The offset array has <c>n_cells + 1</c> entries, so the cell count
    /// is simply one less than its length.
    /// </para>
    /// </summary>
    /// <param name="offsets">VTK offset array.</param>
    /// <returns>The number of cells.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="offsets"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when the offset array is empty.
    /// </exception>
    public static int NumCellsFromOffset(int[] offsets)
    {
        ArgumentNullException.ThrowIfNull(offsets);

        if (offsets.Length == 0)
        {
            throw new ArgumentException("Offset array must not be empty.");
        }

        return offsets.Length - 1;
    }

    /// <summary>
    /// Generates cell type and cell connectivity arrays for the creation of an unstructured
    /// grid from a dictionary mapping cell type identifiers to point index arrays.
    /// <para>
    /// Each key is a cell type code (e.g., 5 for triangles, 10 for tetrahedra) and each
    /// value is a flat or Nx-D array of point indices. Only fixed-size cell types are supported.
    /// </para>
    /// </summary>
    /// <param name="mixedCellDict">
    /// Dictionary mapping cell type codes to arrays of connectivity indices.
    /// Each array must have a length that is a multiple of the cell type's point count.
    /// </param>
    /// <param name="pointsPerCellType">
    /// Dictionary mapping cell type codes to the number of points per cell for that type.
    /// </param>
    /// <param name="totalPoints">
    /// Optional total number of points in the grid. When provided, indices are validated
    /// against this bound.
    /// </param>
    /// <returns>
    /// A tuple of (<c>cellTypes</c>, <c>cellArray</c>) where <c>cellTypes</c> contains one
    /// entry per cell and <c>cellArray</c> is the VTK-style connectivity array.
    /// </returns>
    /// <exception cref="ArgumentException">
    /// Thrown when an unsupported cell type is encountered, array sizes are inconsistent,
    /// or indices are out of range.
    /// </exception>
    public static (byte[] CellTypes, int[] CellArray) CreateMixedCells(
        IDictionary<byte, int[]> mixedCellDict,
        IDictionary<byte, int> pointsPerCellType,
        int? totalPoints = null)
    {
        ArgumentNullException.ThrowIfNull(mixedCellDict);
        ArgumentNullException.ThrowIfNull(pointsPerCellType);

        var allCellTypes = new List<byte>();
        var allCellArr = new List<int>();

        foreach (var (cellType, indices) in mixedCellDict)
        {
            if (!pointsPerCellType.TryGetValue(cellType, out int pointsPerCell))
            {
                throw new ArgumentException(
                    $"Unknown or unsupported cell type: {cellType}.");
            }

            if (pointsPerCell <= 0)
            {
                throw new ArgumentException(
                    $"Cell type {cellType} has variable length and cannot be used with this method.");
            }

            if (indices.Length == 0 || indices.Length % pointsPerCell != 0)
            {
                throw new ArgumentException(
                    $"Expected connectivity array length to be a multiple of {pointsPerCell} " +
                    $"for cell type {cellType}, but got {indices.Length}.");
            }

            int numCells = indices.Length / pointsPerCell;

            for (int c = 0; c < numCells; c++)
            {
                int offset = c * pointsPerCell;

                // Validate indices
                for (int p = 0; p < pointsPerCell; p++)
                {
                    int idx = indices[offset + p];
                    if (idx < 0)
                    {
                        throw new ArgumentException(
                            $"Negative index ({idx}) found for cell type {cellType}.");
                    }

                    if (totalPoints.HasValue && idx >= totalPoints.Value)
                    {
                        throw new ArgumentException(
                            $"Index ({idx}) is out of range (>= {totalPoints.Value}) for cell type {cellType}.");
                    }
                }

                allCellTypes.Add(cellType);
                allCellArr.Add(pointsPerCell);
                for (int p = 0; p < pointsPerCell; p++)
                {
                    allCellArr.Add(indices[offset + p]);
                }
            }
        }

        return (allCellTypes.ToArray(), allCellArr.ToArray());
    }

    /// <summary>
    /// Extracts cells grouped by type from a flat VTK cell connectivity array and
    /// an associated cell types array.
    /// <para>
    /// This is the inverse of <see cref="CreateMixedCells"/>. Only fixed-size cell
    /// types are supported.
    /// </para>
    /// </summary>
    /// <param name="cellTypes">Array of cell type codes, one per cell.</param>
    /// <param name="cells">Flat VTK cell connectivity array.</param>
    /// <param name="pointsPerCellType">
    /// Dictionary mapping cell type codes to the number of points per cell for that type.
    /// </param>
    /// <returns>
    /// A dictionary mapping each unique cell type to a flat array of point indices
    /// for all cells of that type.
    /// </returns>
    /// <exception cref="ArgumentException">
    /// Thrown when an unsupported or variable-size cell type is encountered.
    /// </exception>
    public static Dictionary<byte, int[]> GetCellsByType(
        byte[] cellTypes,
        int[] cells,
        IDictionary<byte, int> pointsPerCellType)
    {
        ArgumentNullException.ThrowIfNull(cellTypes);
        ArgumentNullException.ThrowIfNull(cells);
        ArgumentNullException.ThrowIfNull(pointsPerCellType);

        var result = new Dictionary<byte, List<int>>();

        int index = 0;
        for (int i = 0; i < cellTypes.Length; i++)
        {
            byte cellType = cellTypes[i];

            if (!pointsPerCellType.TryGetValue(cellType, out int pointsPerCell))
            {
                throw new ArgumentException(
                    $"Unknown or unsupported cell type: {cellType}.");
            }

            if (pointsPerCell <= 0)
            {
                throw new ArgumentException(
                    $"Cell type {cellType} has variable length and is not supported.");
            }

            // Skip the count prefix in the connectivity array
            index++;

            if (!result.TryGetValue(cellType, out var list))
            {
                list = new List<int>();
                result[cellType] = list;
            }

            for (int p = 0; p < pointsPerCell; p++)
            {
                list.Add(cells[index + p]);
            }

            index += pointsPerCell;
        }

        return result.ToDictionary(kvp => kvp.Key, kvp => kvp.Value.ToArray());
    }

    /// <summary>
    /// Builds a flat VTK cell connectivity array from a jagged array of cells.
    /// <para>
    /// Each inner array represents a single cell's point indices. The method prepends
    /// the cell size before each cell's indices, producing the standard VTK format.
    /// </para>
    /// </summary>
    /// <param name="cells">
    /// Jagged array where each element contains the point indices for one cell.
    /// </param>
    /// <returns>A flat VTK-style cell connectivity array.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="cells"/> is <c>null</c>.
    /// </exception>
    public static int[] BuildCellArray(int[][] cells)
    {
        ArgumentNullException.ThrowIfNull(cells);

        int totalLength = 0;
        for (int i = 0; i < cells.Length; i++)
        {
            totalLength += 1 + cells[i].Length;
        }

        var result = new int[totalLength];
        int pos = 0;
        for (int i = 0; i < cells.Length; i++)
        {
            result[pos++] = cells[i].Length;
            Array.Copy(cells[i], 0, result, pos, cells[i].Length);
            pos += cells[i].Length;
        }

        return result;
    }

    /// <summary>
    /// Builds an offset array from a flat VTK cell connectivity array.
    /// <para>
    /// The offset array has <c>n_cells + 1</c> entries and records the starting
    /// position of each cell's point indices (after the count prefix) in the
    /// connectivity array.
    /// </para>
    /// </summary>
    /// <param name="cells">A flat VTK cell connectivity array.</param>
    /// <returns>An offset array suitable for modern VTK cell array representations.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="cells"/> is <c>null</c>.
    /// </exception>
    public static int[] BuildOffsetArray(int[] cells)
    {
        ArgumentNullException.ThrowIfNull(cells);

        var offsets = new List<int> { 0 };
        int index = 0;

        while (index < cells.Length)
        {
            int nPoints = cells[index];
            index += nPoints + 1;
            offsets.Add(index);
        }

        return offsets.ToArray();
    }
}
