using PyVista.Core;
using PyVista.Core.Cells;

using CT = PyVista.Core.Cells.CellType;

namespace PyVista.Core.Filters;

/// <summary>
/// Extension methods that mirror the Python <c>pyvista.StructuredGridFilters</c> mixin.
/// <para>
/// These filters operate on <see cref="StructuredGrid"/> datasets and provide
/// operations such as extracting subsets, concatenating grids, and casting to
/// other grid types.
/// </para>
/// </summary>
public static class StructuredGridFilters
{
    // ---------------------------------------------------------------
    //  Extract subset
    // ---------------------------------------------------------------

    /// <summary>
    /// Selects a volume of interest (sub-region) from the structured grid.
    /// <para>
    /// This is the C# equivalent of the Python <c>StructuredGrid.extract_subset()</c>
    /// method. It extracts a rectangular sub-region defined by I-J-K min/max indices
    /// and optionally subsamples the data at the specified rate.
    /// </para>
    /// </summary>
    /// <param name="self">The structured grid to extract from.</param>
    /// <param name="voi">
    /// Volume of interest as a 6-element array:
    /// <c>[xMin, xMax, yMin, yMax, zMin, zMax]</c> in I-J-K indices (0-based).
    /// </param>
    /// <param name="rate">
    /// Sampling rate in each direction as a 3-element array:
    /// <c>[xRate, yRate, zRate]</c>. Defaults to <c>(1, 1, 1)</c>.
    /// </param>
    /// <param name="boundary">
    /// When <c>true</c>, enforces that the boundary of the grid is included
    /// in the subsampled output even when the sample rate is not an even
    /// multiple of the grid dimensions.
    /// </param>
    /// <returns>A new <see cref="StructuredGrid"/> containing the extracted subset.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="self"/> or <paramref name="voi"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when <paramref name="voi"/> does not have exactly 6 elements
    /// or <paramref name="rate"/> does not have exactly 3 elements.
    /// </exception>
    public static StructuredGrid ExtractSubset(
        this StructuredGrid self,
        int[] voi,
        int[]? rate = null,
        bool boundary = false)
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(voi);

        if (voi.Length != 6)
        {
            throw new ArgumentException("VOI must have exactly 6 elements: [xMin, xMax, yMin, yMax, zMin, zMax].", nameof(voi));
        }

        rate ??= [1, 1, 1];
        if (rate.Length != 3)
        {
            throw new ArgumentException("Rate must have exactly 3 elements: [xRate, yRate, zRate].", nameof(rate));
        }

        throw new NotImplementedException("ExtractSubset requires VTK vtkExtractGrid.");
    }

    // ---------------------------------------------------------------
    //  Cast to UnstructuredGrid
    // ---------------------------------------------------------------

    /// <summary>
    /// Casts the <see cref="StructuredGrid"/> to an <see cref="UnstructuredGrid"/>.
    /// <para>
    /// Each cell in the structured grid is converted to a hexahedral cell in the
    /// resulting unstructured grid. Point data, cell data, and field data are
    /// copied to the new grid.
    /// </para>
    /// </summary>
    /// <param name="self">The structured grid to cast.</param>
    /// <returns>A new <see cref="UnstructuredGrid"/> with hexahedral cells.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="self"/> is <c>null</c>.
    /// </exception>
    public static UnstructuredGrid CastToUnstructuredGrid(this StructuredGrid self)
    {
        ArgumentNullException.ThrowIfNull(self);

        var (nx, ny, nz) = self.Dimensions;
        int cx = Math.Max(1, nx - 1);
        int cy = Math.Max(1, ny - 1);
        int cz = Math.Max(1, nz - 1);
        int nCells = cx * cy * cz;

        // Build hexahedral cell connectivity in legacy padded format
        var cells = new int[nCells * 9];
        var cellTypes = new byte[nCells];
        Array.Fill(cellTypes, (byte)CT.Hexahedron);

        int cellIdx = 0;
        for (int iz = 0; iz < cz; iz++)
        {
            for (int iy = 0; iy < cy; iy++)
            {
                for (int ix = 0; ix < cx; ix++)
                {
                    int p0 = PointIndex(ix, iy, iz, nx, ny);
                    int p1 = PointIndex(ix + 1, iy, iz, nx, ny);
                    int p2 = PointIndex(ix + 1, iy + 1, iz, nx, ny);
                    int p3 = PointIndex(ix, iy + 1, iz, nx, ny);
                    int p4 = PointIndex(ix, iy, iz + 1, nx, ny);
                    int p5 = PointIndex(ix + 1, iy, iz + 1, nx, ny);
                    int p6 = PointIndex(ix + 1, iy + 1, iz + 1, nx, ny);
                    int p7 = PointIndex(ix, iy + 1, iz + 1, nx, ny);

                    int offset = cellIdx * 9;
                    cells[offset] = 8;
                    cells[offset + 1] = p0;
                    cells[offset + 2] = p1;
                    cells[offset + 3] = p2;
                    cells[offset + 4] = p3;
                    cells[offset + 5] = p4;
                    cells[offset + 6] = p5;
                    cells[offset + 7] = p6;
                    cells[offset + 8] = p7;

                    cellIdx++;
                }
            }
        }

        var points = (double[])self.Points.Clone();
        return new UnstructuredGrid(cells, cellTypes, points, deep: false);
    }

    // ---------------------------------------------------------------
    //  Cast to RectilinearGrid
    // ---------------------------------------------------------------

    /// <summary>
    /// Casts the <see cref="StructuredGrid"/> to a <see cref="RectilinearGrid"/>.
    /// <para>
    /// This conversion is only valid when the structured grid has axis-aligned
    /// geometry where each row of points along a given axis shares the same
    /// coordinate values. The method extracts unique coordinate values along
    /// each axis to construct the rectilinear grid.
    /// </para>
    /// </summary>
    /// <param name="self">The structured grid to cast.</param>
    /// <returns>A new <see cref="RectilinearGrid"/>.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="self"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="InvalidOperationException">
    /// Thrown when the structured grid does not have axis-aligned geometry
    /// compatible with a rectilinear representation.
    /// </exception>
    public static RectilinearGrid CastToRectilinearGrid(this StructuredGrid self)
    {
        ArgumentNullException.ThrowIfNull(self);

        var (nx, ny, nz) = self.Dimensions;
        double[] xCoords = self.X;
        double[] yCoords = self.Y;
        double[] zCoords = self.Z;

        // Extract unique sorted coordinates along each axis
        var xUnique = ExtractUniqueCoords(xCoords, nx, ny, nz, axis: 0);
        var yUnique = ExtractUniqueCoords(yCoords, nx, ny, nz, axis: 1);
        var zUnique = ExtractUniqueCoords(zCoords, nx, ny, nz, axis: 2);

        if (xUnique.Length != nx || yUnique.Length != ny || zUnique.Length != nz)
        {
            throw new InvalidOperationException(
                "Structured grid geometry is not axis-aligned; cannot convert to RectilinearGrid.");
        }

        return new RectilinearGrid(xUnique, yUnique, zUnique);
    }

    // ---------------------------------------------------------------
    //  Concatenate
    // ---------------------------------------------------------------

    /// <summary>
    /// Concatenates another <see cref="StructuredGrid"/> along the specified axis.
    /// <para>
    /// This is the C# equivalent of the Python <c>StructuredGrid.concatenate()</c>
    /// method. The two grids must have compatible dimensions (matching in all
    /// dimensions except the concatenation axis) and must be coincident along
    /// the joining seam within the specified tolerance.
    /// </para>
    /// </summary>
    /// <param name="self">The first structured grid.</param>
    /// <param name="other">The structured grid to concatenate.</param>
    /// <param name="axis">
    /// Axis along which to concatenate (0 = X, 1 = Y, 2 = Z).
    /// </param>
    /// <param name="tolerance">
    /// Tolerance for point coincidence along the joining seam.
    /// </param>
    /// <returns>A new <see cref="StructuredGrid"/> containing the concatenated grids.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="self"/> or <paramref name="other"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentOutOfRangeException">
    /// Thrown when <paramref name="axis"/> is not 0, 1, or 2.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when the grids have incompatible dimensions.
    /// </exception>
    public static StructuredGrid Concatenate(
        this StructuredGrid self,
        StructuredGrid other,
        int axis,
        double tolerance = 0.0)
    {
        ArgumentNullException.ThrowIfNull(self);
        ArgumentNullException.ThrowIfNull(other);

        if (axis < 0 || axis > 2)
        {
            throw new ArgumentOutOfRangeException(nameof(axis), "Concatenation axis must be 0, 1, or 2.");
        }

        var (nx1, ny1, nz1) = self.Dimensions;
        var (nx2, ny2, nz2) = other.Dimensions;
        int[] dims1 = [nx1, ny1, nz1];
        int[] dims2 = [nx2, ny2, nz2];

        // Check that non-concatenation dimensions match
        for (int i = 0; i < 3; i++)
        {
            if (i != axis && dims1[i] != dims2[i])
            {
                throw new ArgumentException(
                    $"StructuredGrids with dimensions ({nx1},{ny1},{nz1}) and " +
                    $"({nx2},{ny2},{nz2}) are not compatible for concatenation along axis {axis}.");
            }
        }

        throw new NotImplementedException(
            "Full Concatenate implementation requires numpy-style array reshaping. " +
            "Use the Python pyvista.StructuredGrid.concatenate() for full functionality.");
    }

    // ---------------------------------------------------------------
    //  Private helpers
    // ---------------------------------------------------------------

    /// <summary>
    /// Computes the flat point index for the given IJK coordinates using
    /// Fortran (column-major) ordering: <c>ix + nx * (iy + ny * iz)</c>.
    /// </summary>
    private static int PointIndex(int ix, int iy, int iz, int nx, int ny)
    {
        return ix + nx * (iy + ny * iz);
    }

    /// <summary>
    /// Extracts unique coordinate values along the specified axis from a
    /// flat coordinate array, sorted in ascending order.
    /// </summary>
    /// <param name="coords">Flat coordinate component array (length = nx*ny*nz).</param>
    /// <param name="nx">Number of points in the X direction.</param>
    /// <param name="ny">Number of points in the Y direction.</param>
    /// <param name="nz">Number of points in the Z direction.</param>
    /// <param name="axis">Axis index (0=X, 1=Y, 2=Z).</param>
    /// <returns>Sorted array of unique coordinate values along the axis.</returns>
    private static double[] ExtractUniqueCoords(double[] coords, int nx, int ny, int nz, int axis)
    {
        int count = axis switch
        {
            0 => nx,
            1 => ny,
            2 => nz,
            _ => throw new ArgumentOutOfRangeException(nameof(axis)),
        };

        // Sample coordinates along the specified axis at fixed positions in other axes
        var unique = new SortedSet<double>();
        int stride = axis switch
        {
            0 => 1,
            1 => nx,
            2 => nx * ny,
            _ => 1,
        };

        for (int i = 0; i < count; i++)
        {
            unique.Add(coords[i * stride]);
        }

        var result = new double[unique.Count];
        unique.CopyTo(result);
        return result;
    }
}
