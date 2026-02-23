namespace PyVista.Core.Cells;

/// <summary>
/// Represents an array of cells, providing convenience methods for constructing and
/// querying cell connectivity data.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.CellArray</c> class. Internally,
/// cells are stored using offset and connectivity arrays following the VTK cell-array
/// convention.
/// </para>
/// <para>
/// The <em>legacy format</em> interleaves cell sizes with point IDs:
/// <c>{ n0, p0_0, p0_1, …, p0_n, n1, p1_0, p1_1, …, p1_n, … }</c>
/// where <c>n0</c> is the number of points in cell 0, and <c>pX_Y</c> is the Y-th
/// point in cell X.
/// </para>
/// </summary>
/// <example>
/// Create a cell array containing two triangles from the traditional interleaved format:
/// <code>
/// var cellArr = new CellArray(new[] { 3, 0, 1, 2, 3, 3, 4, 5 });
/// </code>
/// Create a cell array from separate offsets and connectivity arrays:
/// <code>
/// var cellArr = CellArray.FromArrays(new[] { 0, 3, 6 }, new[] { 0, 1, 2, 3, 4, 5 });
/// </code>
/// </example>
public sealed class CellArray
{
    private int[] _offsets;
    private int[] _connectivity;

    /// <summary>
    /// Initializes a new instance of the <see cref="CellArray"/> class from a legacy-format
    /// cell array.
    /// </summary>
    /// <param name="cells">
    /// A flat array in legacy VTK cell-array format:
    /// <c>{ n0, p0_0, …, p0_n0, n1, p1_0, …, p1_n1, … }</c>.
    /// Pass <c>null</c> or omit to create an empty cell array.
    /// </param>
    /// <exception cref="CellSizeError">
    /// Thrown when the array contains an invalid cell size (e.g. negative or exceeding the
    /// remaining array length).
    /// </exception>
    public CellArray(int[]? cells = null)
    {
        _offsets = new[] { 0 };
        _connectivity = Array.Empty<int>();

        if (cells is not null)
        {
            Cells = cells;
        }
    }

    /// <summary>
    /// Gets or sets the cells in legacy VTK interleaved format.
    /// <para>
    /// The format interleaves the number of points in each cell with the point IDs:
    /// <c>{ n0, p0_0, …, p0_n0, n1, p1_0, …, p1_n1, … }</c>.
    /// </para>
    /// </summary>
    /// <exception cref="CellSizeError">
    /// Thrown when the supplied array has an invalid layout.
    /// </exception>
    public int[] Cells
    {
        get
        {
            // Export to legacy format: interleave size + connectivity per cell.
            int nCells = NCells;
            int totalLength = nCells + _connectivity.Length;
            int[] result = new int[totalLength];
            int pos = 0;

            for (int c = 0; c < nCells; c++)
            {
                int start = _offsets[c];
                int end = _offsets[c + 1];
                int cellSize = end - start;
                result[pos++] = cellSize;
                for (int j = start; j < end; j++)
                {
                    result[pos++] = _connectivity[j];
                }
            }

            return result;
        }
        set
        {
            ArgumentNullException.ThrowIfNull(value);
            ImportLegacyFormat(value);
        }
    }

    /// <summary>
    /// Gets the number of cells in the array.
    /// </summary>
    public int NCells => _offsets.Length - 1;

    /// <summary>
    /// Gets the connectivity array containing the point IDs that define each cell.
    /// <para>
    /// The connectivity is a flat array of all point IDs concatenated. Use
    /// <see cref="OffsetArray"/> to determine where each cell's IDs begin and end.
    /// </para>
    /// </summary>
    /// <returns>A copy of the internal connectivity array.</returns>
    public int[] ConnectivityArray => (int[])_connectivity.Clone();

    /// <summary>
    /// Gets the offset array used to index into the <see cref="ConnectivityArray"/>.
    /// <para>
    /// The offset array has length <c><see cref="NCells"/> + 1</c>. Cell <c>i</c> uses
    /// connectivity indices from <c>offsets[i]</c> (inclusive) to <c>offsets[i+1]</c>
    /// (exclusive).
    /// </para>
    /// </summary>
    /// <returns>A copy of the internal offset array.</returns>
    public int[] OffsetArray => (int[])_offsets.Clone();

    /// <summary>
    /// Gets a two-dimensional array of shape <c>(NCells, cellSize)</c> for cell arrays
    /// where every cell has the same number of points.
    /// <para>
    /// This property does not validate that all cells are actually the same size. If they
    /// are not, the behaviour is undefined and may throw an exception or return incorrect
    /// data.
    /// </para>
    /// </summary>
    /// <returns>
    /// A two-dimensional <see cref="int"/> array of shape <c>(NCells, cellSize)</c>, or an
    /// empty array if the cell array contains no cells.
    /// </returns>
    /// <exception cref="InvalidOperationException">
    /// Thrown when the connectivity length is not evenly divisible by the cell size.
    /// </exception>
    public int[,] RegularCells
    {
        get
        {
            if (_connectivity.Length == 0)
            {
                return new int[0, 0];
            }

            int cellSize = _offsets[1] - _offsets[0];
            if (cellSize == 0)
            {
                return new int[NCells, 0];
            }

            int nCells = NCells;
            if (_connectivity.Length % cellSize != 0)
            {
                throw new InvalidOperationException(
                    $"Connectivity length ({_connectivity.Length}) is not divisible by cell size ({cellSize}).");
            }

            var result = new int[nCells, cellSize];
            for (int i = 0; i < nCells; i++)
            {
                int offset = _offsets[i];
                for (int j = 0; j < cellSize; j++)
                {
                    result[i, j] = _connectivity[offset + j];
                }
            }

            return result;
        }
    }

    /// <summary>
    /// Constructs a <see cref="CellArray"/> from separate offset and connectivity arrays.
    /// </summary>
    /// <param name="offsets">
    /// Offset array of length <c>nCells + 1</c>. Each pair of consecutive values defines the
    /// start and end indices into <paramref name="connectivity"/>.
    /// </param>
    /// <param name="connectivity">
    /// Flat array of point IDs.
    /// </param>
    /// <param name="deep">
    /// When <c>true</c>, the input arrays are deep-copied. When <c>false</c> (the default),
    /// the cell array takes ownership of the references.
    /// </param>
    /// <returns>A new <see cref="CellArray"/>.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="offsets"/> or <paramref name="connectivity"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when <paramref name="offsets"/> has fewer than 2 elements or does not start at 0.
    /// </exception>
    public static CellArray FromArrays(int[] offsets, int[] connectivity, bool deep = false)
    {
        ArgumentNullException.ThrowIfNull(offsets);
        ArgumentNullException.ThrowIfNull(connectivity);

        if (offsets.Length < 1)
        {
            throw new ArgumentException("Offsets array must have at least 1 element.", nameof(offsets));
        }

        var cellArr = new CellArray();
        cellArr.SetData(offsets, connectivity, deep);
        return cellArr;
    }

    /// <summary>
    /// Constructs a <see cref="CellArray"/> from a regular (equal-sized) cell array.
    /// </summary>
    /// <param name="cells">
    /// A two-dimensional array of shape <c>(nCells, cellSize)</c> where every cell has the
    /// same number of points.
    /// </param>
    /// <param name="deep">
    /// When <c>true</c>, the connectivity is deep-copied. When <c>false</c> (the default),
    /// the cell array takes ownership of the flattened array.
    /// </param>
    /// <returns>A new <see cref="CellArray"/>.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="cells"/> is <c>null</c>.
    /// </exception>
    public static CellArray FromRegularCells(int[,] cells, bool deep = false)
    {
        ArgumentNullException.ThrowIfNull(cells);

        int nCells = cells.GetLength(0);
        int cellSize = cells.GetLength(1);

        int[] offsets = new int[nCells + 1];
        for (int i = 0; i <= nCells; i++)
        {
            offsets[i] = cellSize * i;
        }

        int[] connectivity = new int[nCells * cellSize];
        int pos = 0;
        for (int i = 0; i < nCells; i++)
        {
            for (int j = 0; j < cellSize; j++)
            {
                connectivity[pos++] = cells[i, j];
            }
        }

        var cellArr = new CellArray();
        cellArr.SetData(offsets, connectivity, deep);
        return cellArr;
    }

    /// <summary>
    /// Constructs a <see cref="CellArray"/> from an array of variable-length cells.
    /// </summary>
    /// <param name="cells">
    /// A jagged array where each element contains the point IDs for one cell.
    /// Cells may have different sizes.
    /// </param>
    /// <returns>A new <see cref="CellArray"/>.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="cells"/> is <c>null</c>.
    /// </exception>
    public static CellArray FromIrregularCells(int[][] cells)
    {
        ArgumentNullException.ThrowIfNull(cells);

        int nCells = cells.Length;
        int[] offsets = new int[nCells + 1];
        int total = 0;
        for (int i = 0; i < nCells; i++)
        {
            offsets[i] = total;
            total += cells[i].Length;
        }

        offsets[nCells] = total;

        int[] connectivity = new int[total];
        int pos = 0;
        for (int i = 0; i < nCells; i++)
        {
            Array.Copy(cells[i], 0, connectivity, pos, cells[i].Length);
            pos += cells[i].Length;
        }

        return FromArrays(offsets, connectivity);
    }

    /// <summary>
    /// Returns a string representation of this cell array.
    /// </summary>
    /// <returns>A formatted string describing the cell array.</returns>
    public override string ToString()
    {
        return $"CellArray (NCells={NCells}, ConnectivityLength={_connectivity.Length})";
    }

    // ------------------------------------------------------------------
    //  Private helpers
    // ------------------------------------------------------------------

    /// <summary>
    /// Sets the internal offset and connectivity arrays.
    /// </summary>
    private void SetData(int[] offsets, int[] connectivity, bool deep)
    {
        _offsets = deep ? (int[])offsets.Clone() : offsets;
        _connectivity = deep ? (int[])connectivity.Clone() : connectivity;
    }

    /// <summary>
    /// Imports cell data from the legacy VTK interleaved format.
    /// </summary>
    private void ImportLegacyFormat(int[] cells)
    {
        // First pass: count cells and validate.
        var offsetList = new List<int> { 0 };
        int i = 0;
        while (i < cells.Length)
        {
            int nPts = cells[i];
            if (nPts < 0 || i + 1 + nPts > cells.Length)
            {
                throw new CellSizeError(
                    $"Cell array size is invalid. Encountered cell size {nPts} at position {i} " +
                    $"but only {cells.Length - i - 1} values remain. This is likely due to " +
                    "an invalid connectivity array.");
            }

            offsetList.Add(offsetList[^1] + nPts);
            i += 1 + nPts;
        }

        if (i != cells.Length)
        {
            throw new CellSizeError(
                $"Cell array size is invalid. Size ({cells.Length}) does not match expected " +
                $"size ({i}). This is likely due to an invalid connectivity array.");
        }

        // Second pass: extract connectivity.
        int nCells = offsetList.Count - 1;
        int totalConnectivity = offsetList[^1];
        int[] connectivity = new int[totalConnectivity];
        int pos = 0;
        int src = 0;

        for (int c = 0; c < nCells; c++)
        {
            int cellSize = cells[src];
            src++;
            for (int j = 0; j < cellSize; j++)
            {
                connectivity[pos++] = cells[src++];
            }
        }

        _offsets = offsetList.ToArray();
        _connectivity = connectivity;
    }
}
