using System.Globalization;

using PyVista.Core.Cells;

using CT = PyVista.Core.Cells.CellType;

namespace PyVista.Core;

/// <summary>
/// Dataset used for topologically regular arrays of data.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.StructuredGrid</c> class.
/// A structured grid has implicit topology defined by its <see cref="Dimensions"/>,
/// with explicit point coordinates stored in <see cref="DataSet.Points"/>.
/// </para>
/// </summary>
public class StructuredGrid : PointSet
{
    private int _dimX = 1;
    private int _dimY = 1;
    private int _dimZ = 1;

    /// <summary>
    /// Initializes a new instance of the <see cref="StructuredGrid"/> class (empty).
    /// </summary>
    public StructuredGrid()
    {
    }

    /// <summary>
    /// Initializes a new instance of the <see cref="StructuredGrid"/> class from
    /// separate X, Y, and Z coordinate arrays.
    /// <para>
    /// All three arrays must have the same shape. The shape defines the grid
    /// <see cref="Dimensions"/> and the points are assembled from the flattened arrays.
    /// </para>
    /// </summary>
    /// <param name="x">X-coordinates as a flat Fortran-order array.</param>
    /// <param name="y">Y-coordinates as a flat Fortran-order array.</param>
    /// <param name="z">Z-coordinates as a flat Fortran-order array.</param>
    /// <param name="dimX">Number of points in the X direction.</param>
    /// <param name="dimY">Number of points in the Y direction.</param>
    /// <param name="dimZ">Number of points in the Z direction.</param>
    /// <exception cref="ArgumentException">
    /// Thrown when the input arrays have different lengths or the dimensions do not match.
    /// </exception>
    public StructuredGrid(double[] x, double[] y, double[] z, int dimX, int dimY, int dimZ)
    {
        ArgumentNullException.ThrowIfNull(x);
        ArgumentNullException.ThrowIfNull(y);
        ArgumentNullException.ThrowIfNull(z);

        if (x.Length != y.Length || y.Length != z.Length)
        {
            throw new ArgumentException("Input point coordinate arrays must have the same length.");
        }

        int expectedSize = dimX * dimY * dimZ;
        if (x.Length != expectedSize)
        {
            throw new ArgumentException(
                $"Array length ({x.Length}) does not match dimensions ({dimX}×{dimY}×{dimZ}={expectedSize}).");
        }

        _dimX = dimX;
        _dimY = dimY;
        _dimZ = dimZ;

        var pts = new double[x.Length * 3];
        for (int i = 0; i < x.Length; i++)
        {
            int offset = i * 3;
            pts[offset] = x[i];
            pts[offset + 1] = y[i];
            pts[offset + 2] = z[i];
        }

        Points = pts;
    }

    /// <summary>
    /// Initializes a new instance of the <see cref="StructuredGrid"/> class from a flat
    /// points array and explicit dimensions.
    /// </summary>
    /// <param name="points">Flat row-major point coordinates (length = dimX × dimY × dimZ × 3).</param>
    /// <param name="dimX">Number of points in the X direction.</param>
    /// <param name="dimY">Number of points in the Y direction.</param>
    /// <param name="dimZ">Number of points in the Z direction.</param>
    /// <param name="deep">When <c>true</c>, copies the point array.</param>
    public StructuredGrid(double[] points, int dimX, int dimY, int dimZ, bool deep = false)
        : base(points, deep)
    {
        int expectedPoints = dimX * dimY * dimZ;
        if (NPoints != expectedPoints)
        {
            throw new ArgumentException(
                $"Number of points ({NPoints}) does not match dimensions ({dimX}×{dimY}×{dimZ}={expectedPoints}).");
        }

        _dimX = dimX;
        _dimY = dimY;
        _dimZ = dimZ;
    }

    // ---------------------------------------------------------------
    //  Dimensions
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets or sets the grid dimensions as a tuple (NX, NY, NZ).
    /// <para>
    /// Dimensions define the number of points in each direction. The number of cells
    /// in each direction is one less than the corresponding dimension.
    /// </para>
    /// </summary>
    /// <exception cref="ArgumentException">
    /// Thrown when setting dimensions whose product does not match the number of points.
    /// </exception>
    public (int NX, int NY, int NZ) Dimensions
    {
        get => (_dimX, _dimY, _dimZ);
        set
        {
            int expectedPoints = value.NX * value.NY * value.NZ;
            if (NPoints != 0 && NPoints != expectedPoints)
            {
                throw new ArgumentException(
                    $"Cannot set dimensions ({value.NX}×{value.NY}×{value.NZ}={expectedPoints}) " +
                    $"for grid with {NPoints} points.");
            }

            _dimX = value.NX;
            _dimY = value.NY;
            _dimZ = value.NZ;
        }
    }

    // ---------------------------------------------------------------
    //  X / Y / Z coordinate arrays
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets the X coordinates of all points, reshaped to the grid dimensions
    /// in Fortran (column-major) order.
    /// </summary>
    public double[] X => ExtractCoordinateComponent(0);

    /// <summary>
    /// Gets the Y coordinates of all points, reshaped to the grid dimensions
    /// in Fortran (column-major) order.
    /// </summary>
    public double[] Y => ExtractCoordinateComponent(1);

    /// <summary>
    /// Gets the Z coordinates of all points, reshaped to the grid dimensions
    /// in Fortran (column-major) order.
    /// </summary>
    public double[] Z => ExtractCoordinateComponent(2);

    // ---------------------------------------------------------------
    //  Points matrix
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets the points as a 4-D style array with x/y/z along the last dimension.
    /// <para>
    /// Returns a flat array where the first three logical dimensions correspond to
    /// <see cref="Dimensions"/> in Fortran order, and the last dimension is the
    /// (x, y, z) component.
    /// </para>
    /// </summary>
    public double[] PointsMatrix
    {
        get
        {
            int n = NPoints;
            var result = new double[n * 3];
            Array.Copy(Points, result, n * 3);
            return result;
        }
    }

    // ---------------------------------------------------------------
    //  Implicit cell count
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets the number of cells implied by the grid dimensions.
    /// <para>
    /// For a structured grid, the number of cells is
    /// <c>max(1, NX-1) × max(1, NY-1) × max(1, NZ-1)</c>.
    /// </para>
    /// </summary>
    public int NCellsFromDimensions
    {
        get
        {
            int cx = Math.Max(1, _dimX - 1);
            int cy = Math.Max(1, _dimY - 1);
            int cz = Math.Max(1, _dimZ - 1);
            return cx * cy * cz;
        }
    }

    // ---------------------------------------------------------------
    //  Hide cells / points
    // ---------------------------------------------------------------

    /// <summary>
    /// Hides cells without deleting them by setting ghost cell flags.
    /// <para>
    /// Creates a cell data array named <c>"vtkGhostType"</c> with the hidden cell flag
    /// set at the specified indices.
    /// </para>
    /// </summary>
    /// <param name="indices">Cell indices to hide.</param>
    /// <param name="inplace">When <c>true</c>, modifies this grid; otherwise returns a copy.</param>
    /// <returns>The grid with hidden cells.</returns>
    public StructuredGrid HideCells(int[] indices, bool inplace = false)
    {
        ArgumentNullException.ThrowIfNull(indices);
        StructuredGrid target = inplace ? this : (StructuredGrid)Copy(deep: true);
        int nCells = target.NCellsFromDimensions;
        var ghostCells = new double[nCells];

        foreach (int idx in indices)
        {
            if (idx >= 0 && idx < nCells)
            {
                ghostCells[idx] = 2; // HIDDENCELL flag value
            }
        }

        target.CellData.SetArray(ghostCells, "vtkGhostType");
        return target;
    }

    /// <summary>
    /// Hides points without deleting them by setting ghost point flags.
    /// <para>
    /// Creates a point data array named <c>"vtkGhostType"</c> with the hidden point flag
    /// set at the specified indices.
    /// </para>
    /// </summary>
    /// <param name="indices">Point indices to hide.</param>
    public void HidePoints(int[] indices)
    {
        ArgumentNullException.ThrowIfNull(indices);
        int n = NPoints;
        var ghostPoints = new double[n];

        foreach (int idx in indices)
        {
            if (idx >= 0 && idx < n)
            {
                ghostPoints[idx] = 1; // HIDDENPOINT flag value
            }
        }

        PointData.SetArray(ghostPoints, "vtkGhostType");
    }

    // ---------------------------------------------------------------
    //  Copy helpers
    // ---------------------------------------------------------------

    /// <inheritdoc />
    public override void DeepCopy(DataObject source)
    {
        base.DeepCopy(source);
        if (source is StructuredGrid sg)
        {
            _dimX = sg._dimX;
            _dimY = sg._dimY;
            _dimZ = sg._dimZ;
        }
    }

    /// <inheritdoc />
    public override void ShallowCopy(DataObject source)
    {
        base.ShallowCopy(source);
        if (source is StructuredGrid sg)
        {
            _dimX = sg._dimX;
            _dimY = sg._dimY;
            _dimZ = sg._dimZ;
        }
    }

    /// <inheritdoc />
    protected override List<(string Name, string Value)> GetAttributes()
    {
        var attrs = base.GetAttributes();
        attrs.Add(("Dimensions", $"{_dimX}, {_dimY}, {_dimZ}"));
        return attrs;
    }

    // ---------------------------------------------------------------
    //  Private helpers
    // ---------------------------------------------------------------

    /// <summary>
    /// Extracts a single coordinate component (0=X, 1=Y, 2=Z) from the points array.
    /// </summary>
    private double[] ExtractCoordinateComponent(int component)
    {
        int n = NPoints;
        var result = new double[n];
        var pts = Points;
        for (int i = 0; i < n; i++)
        {
            result[i] = pts[i * 3 + component];
        }

        return result;
    }
}

/// <summary>
/// Explicit structured grid dataset.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.ExplicitStructuredGrid</c> class.
/// An explicit structured grid extends the <see cref="PointSet"/> with explicit cell connectivity
/// and topological dimensions, supporting cells that may have degenerate or irregular geometry.
/// </para>
/// </summary>
public class ExplicitStructuredGrid : PointSet
{
    private CellArray _cellArray = new();
    private int _dimX = 1;
    private int _dimY = 1;
    private int _dimZ = 1;

    /// <summary>
    /// Initializes a new instance of the <see cref="ExplicitStructuredGrid"/> class (empty).
    /// </summary>
    public ExplicitStructuredGrid()
    {
    }

    /// <summary>
    /// Initializes a new instance of the <see cref="ExplicitStructuredGrid"/> class
    /// from dimensions, cell connectivity, and points.
    /// </summary>
    /// <param name="dimensions">A tuple of (NI, NJ, NK) topological dimensions.</param>
    /// <param name="cells">
    /// Padded cell connectivity array in legacy format. Each cell is expected to be a hexahedron
    /// (8 points per cell).
    /// </param>
    /// <param name="points">Flat row-major point coordinates.</param>
    public ExplicitStructuredGrid((int NI, int NJ, int NK) dimensions, int[] cells, double[] points)
        : base(points, deep: false)
    {
        ArgumentNullException.ThrowIfNull(cells);
        _dimX = dimensions.NI;
        _dimY = dimensions.NJ;
        _dimZ = dimensions.NK;
        _cellArray = new CellArray(cells);
    }

    /// <summary>
    /// Initializes a new instance of the <see cref="ExplicitStructuredGrid"/> class
    /// from dimensions and corner coordinates.
    /// </summary>
    /// <param name="dimensions">A tuple of (NI, NJ, NK) topological dimensions.</param>
    /// <param name="corners">
    /// Flat row-major corner coordinates. The number of corners must match
    /// <c>2*(NI-1) × 2*(NJ-1) × 2*(NK-1)</c> corner points with 3 components each.
    /// </param>
    public ExplicitStructuredGrid((int NI, int NJ, int NK) dimensions, double[] corners)
    {
        ArgumentNullException.ThrowIfNull(corners);
        _dimX = dimensions.NI;
        _dimY = dimensions.NJ;
        _dimZ = dimensions.NK;

        int ni = _dimX - 1;
        int nj = _dimY - 1;
        int nk = _dimZ - 1;
        int nCells = ni * nj * nk;

        // Build simple sequential cell connectivity (8 points per hexahedron)
        int[] cells = new int[nCells * 9];
        for (int c = 0; c < nCells; c++)
        {
            cells[c * 9] = 8;
            for (int p = 0; p < 8; p++)
            {
                cells[c * 9 + 1 + p] = c * 8 + p;
            }
        }

        _cellArray = new CellArray(cells);
        Points = corners;
    }

    // ---------------------------------------------------------------
    //  Dimensions
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets the topological dimensions of the grid.
    /// </summary>
    /// <returns>
    /// A tuple of the number of sampling points in the I, J, and K directions.
    /// </returns>
    public (int NI, int NJ, int NK) Dimensions => (_dimX, _dimY, _dimZ);

    // ---------------------------------------------------------------
    //  Cells
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets the cell connectivity in legacy padded format.
    /// </summary>
    public int[] Cells => _cellArray.Cells;

    /// <summary>
    /// Gets the number of cells.
    /// </summary>
    public int NCellsTotal => _cellArray.NCells;

    // ---------------------------------------------------------------
    //  Cast to UnstructuredGrid
    // ---------------------------------------------------------------

    /// <summary>
    /// Casts this explicit structured grid to an <see cref="UnstructuredGrid"/>.
    /// <para>
    /// All cells are treated as hexahedra. The BLOCK_I, BLOCK_J, and BLOCK_K
    /// cell arrays are added to enable restoring the explicit structured grid.
    /// </para>
    /// </summary>
    /// <returns>A new <see cref="UnstructuredGrid"/>.</returns>
    public new UnstructuredGrid CastToUnstructuredGrid()
    {
        int nCells = _cellArray.NCells;
        var cellTypesArr = new byte[nCells];
        Array.Fill(cellTypesArr, (byte)CT.Hexahedron);

        var ugrid = new UnstructuredGrid(_cellArray.Cells, cellTypesArr, (double[])Points.Clone(), deep: false);

        // Add BLOCK_I, BLOCK_J, BLOCK_K cell arrays
        int ni = Math.Max(1, _dimX - 1);
        int nj = Math.Max(1, _dimY - 1);
        int nk = Math.Max(1, _dimZ - 1);

        var blockI = new double[nCells];
        var blockJ = new double[nCells];
        var blockK = new double[nCells];

        for (int idx = 0; idx < nCells; idx++)
        {
            int remainder = idx;
            blockI[idx] = remainder % ni;
            remainder /= ni;
            blockJ[idx] = remainder % nj;
            remainder /= nj;
            blockK[idx] = remainder;
        }

        ugrid.CellData.SetArray(blockI, "BLOCK_I");
        ugrid.CellData.SetArray(blockJ, "BLOCK_J");
        ugrid.CellData.SetArray(blockK, "BLOCK_K");

        return ugrid;
    }

    // ---------------------------------------------------------------
    //  Hide / show cells
    // ---------------------------------------------------------------

    /// <summary>
    /// Hides specific cells by setting the ghost cell array to HIDDENCELL.
    /// </summary>
    /// <param name="indices">Cell indices to hide.</param>
    /// <param name="inplace">When <c>true</c>, modifies this grid; otherwise returns a copy.</param>
    /// <returns>The grid with hidden cells.</returns>
    public ExplicitStructuredGrid HideCells(int[] indices, bool inplace = false)
    {
        ArgumentNullException.ThrowIfNull(indices);
        ExplicitStructuredGrid target = inplace ? this : (ExplicitStructuredGrid)Copy(deep: true);
        int nCells = target._cellArray.NCells;
        var ghostCells = new double[nCells];

        foreach (int idx in indices)
        {
            if (idx >= 0 && idx < nCells)
            {
                ghostCells[idx] = 2; // HIDDENCELL flag value
            }
        }

        target.CellData.SetArray(ghostCells, "vtkGhostType");
        return target;
    }

    /// <summary>
    /// Shows hidden cells by clearing the HIDDENCELL ghost cell flag.
    /// </summary>
    /// <param name="inplace">When <c>true</c>, modifies this grid; otherwise returns a copy.</param>
    /// <returns>The grid with all cells visible.</returns>
    public ExplicitStructuredGrid ShowCells(bool inplace = false)
    {
        ExplicitStructuredGrid target = inplace ? this : (ExplicitStructuredGrid)Copy(deep: true);
        if (target.CellData.ContainsKey("vtkGhostType"))
        {
            var array = target.CellData["vtkGhostType"];
            for (int i = 0; i < array.Length; i++)
            {
                if ((int)array[i] == 2) // HIDDENCELL
                {
                    array[i] = 0;
                }
            }
        }

        return target;
    }

    // ---------------------------------------------------------------
    //  Visible bounds
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets the bounding box of the visible cells.
    /// <para>
    /// Different from <see cref="DataSet.Bounds"/> which returns the bounding box of the
    /// complete grid, this property returns the bounding box considering only cells
    /// whose ghost cell flag is not HIDDENCELL.
    /// </para>
    /// </summary>
    public BoundsTuple VisibleBounds
    {
        get
        {
            if (!CellData.ContainsKey("vtkGhostType"))
            {
                return Bounds;
            }

            var ghost = CellData["vtkGhostType"];
            var offsets = _cellArray.OffsetArray;
            var conn = _cellArray.ConnectivityArray;
            var pts = Points;

            double xMin = double.MaxValue, xMax = double.MinValue;
            double yMin = double.MaxValue, yMax = double.MinValue;
            double zMin = double.MaxValue, zMax = double.MinValue;
            bool anyVisible = false;

            int nCells = _cellArray.NCells;
            for (int c = 0; c < nCells; c++)
            {
                if ((int)ghost[c] == 2) continue; // HIDDENCELL

                anyVisible = true;
                int start = offsets[c];
                int end = offsets[c + 1];
                for (int j = start; j < end; j++)
                {
                    int pid = conn[j];
                    int pOff = pid * 3;
                    double x = pts[pOff], y = pts[pOff + 1], z = pts[pOff + 2];
                    if (x < xMin) xMin = x;
                    if (x > xMax) xMax = x;
                    if (y < yMin) yMin = y;
                    if (y > yMax) yMax = y;
                    if (z < zMin) zMin = z;
                    if (z > zMax) zMax = z;
                }
            }

            return anyVisible
                ? new BoundsTuple(xMin, xMax, yMin, yMax, zMin, zMax)
                : Bounds;
        }
    }

    // ---------------------------------------------------------------
    //  Cell ID / coords
    // ---------------------------------------------------------------

    /// <summary>
    /// Returns the cell ID for the given structured coordinates (i, j, k).
    /// </summary>
    /// <param name="i">I-index.</param>
    /// <param name="j">J-index.</param>
    /// <param name="k">K-index.</param>
    /// <returns>The cell ID, or <c>-1</c> if the coordinates are outside the grid extent.</returns>
    public int CellId(int i, int j, int k)
    {
        int ni = _dimX - 1;
        int nj = _dimY - 1;
        int nk = _dimZ - 1;

        if (i < 0 || i >= ni || j < 0 || j >= nj || k < 0 || k >= nk)
        {
            return -1;
        }

        return i + ni * (j + nj * k);
    }

    /// <summary>
    /// Returns the structured coordinates (i, j, k) for the given cell ID.
    /// </summary>
    /// <param name="cellId">The cell ID.</param>
    /// <returns>A tuple of (I, J, K) coordinates, or <c>null</c> if the ID is invalid.</returns>
    public (int I, int J, int K)? CellCoords(int cellId)
    {
        int ni = _dimX - 1;
        int nj = _dimY - 1;
        int nk = _dimZ - 1;
        int total = ni * nj * nk;

        if (cellId < 0 || cellId >= total)
        {
            return null;
        }

        int remainder = cellId;
        int ci = remainder % ni;
        remainder /= ni;
        int cj = remainder % nj;
        int ck = remainder / nj;

        return (ci, cj, ck);
    }

    // ---------------------------------------------------------------
    //  Compute connectivity / connections
    // ---------------------------------------------------------------

    /// <summary>
    /// Computes the faces connectivity flags array.
    /// <para>
    /// Each value is interpreted as a 6-bit binary number indicating face connectivity
    /// with topological neighbors -Z, +Z, -Y, +Y, -X, +X.
    /// The result is stored in the <c>"ConnectivityFlags"</c> cell data array.
    /// </para>
    /// </summary>
    /// <param name="inplace">When <c>true</c>, modifies this grid; otherwise returns a copy.</param>
    /// <returns>The grid with connectivity flags computed.</returns>
    public ExplicitStructuredGrid ComputeConnectivity(bool inplace = false)
    {
        ExplicitStructuredGrid target = inplace ? this : (ExplicitStructuredGrid)Copy(deep: true);
        int nCells = target._cellArray.NCells;
        var flags = new double[nCells];

        int ni = _dimX - 1;
        int nj = _dimY - 1;

        // Simplified connectivity: mark topological neighbors as connected
        for (int c = 0; c < nCells; c++)
        {
            var coords = target.CellCoords(c);
            if (coords is null) continue;

            var (ci, cj, ck) = coords.Value;
            int flag = 0;

            // -X neighbor
            if (ci > 0) flag |= 1;
            // +X neighbor
            if (ci < ni - 1) flag |= 2;
            // -Y neighbor
            if (cj > 0) flag |= 4;
            // +Y neighbor
            if (cj < nj - 1) flag |= 8;
            // -Z neighbor
            if (ck > 0) flag |= 16;
            // +Z neighbor
            int nk = _dimZ - 1;
            if (ck < nk - 1) flag |= 32;

            flags[c] = flag;
        }

        target.CellData.SetArray(flags, "ConnectivityFlags");
        return target;
    }

    /// <summary>
    /// Computes an array with the number of connected cell faces.
    /// <para>
    /// The results are stored in the <c>"number_of_connections"</c> cell data array.
    /// </para>
    /// </summary>
    /// <param name="inplace">When <c>true</c>, modifies this grid; otherwise returns a copy.</param>
    /// <returns>The grid with connection counts computed.</returns>
    public ExplicitStructuredGrid ComputeConnections(bool inplace = false)
    {
        ExplicitStructuredGrid target = inplace ? this : (ExplicitStructuredGrid)Copy(deep: true);

        double[]? flagsArr = null;
        if (target.CellData.ContainsKey("ConnectivityFlags"))
        {
            flagsArr = target.CellData["ConnectivityFlags"];
        }
        else
        {
            var temp = target.ComputeConnectivity(inplace: false);
            flagsArr = temp.CellData["ConnectivityFlags"];
        }

        int nCells = flagsArr.Length;
        var connections = new double[nCells];
        for (int i = 0; i < nCells; i++)
        {
            int val = (int)flagsArr[i];
            int count = 0;
            while (val > 0)
            {
                count += val & 1;
                val >>= 1;
            }

            connections[i] = count;
        }

        target.CellData.SetArray(connections, "number_of_connections");
        return target;
    }

    // ---------------------------------------------------------------
    //  Copy helpers
    // ---------------------------------------------------------------

    /// <inheritdoc />
    public override void DeepCopy(DataObject source)
    {
        base.DeepCopy(source);
        if (source is ExplicitStructuredGrid esg)
        {
            _cellArray = new CellArray(esg._cellArray.Cells);
            _dimX = esg._dimX;
            _dimY = esg._dimY;
            _dimZ = esg._dimZ;
        }
    }

    /// <inheritdoc />
    public override void ShallowCopy(DataObject source)
    {
        base.ShallowCopy(source);
        if (source is ExplicitStructuredGrid esg)
        {
            _cellArray = esg._cellArray;
            _dimX = esg._dimX;
            _dimY = esg._dimY;
            _dimZ = esg._dimZ;
        }
    }

    /// <inheritdoc />
    protected override List<(string Name, string Value)> GetAttributes()
    {
        var attrs = base.GetAttributes();
        attrs.Add(("Dimensions", $"{_dimX}, {_dimY}, {_dimZ}"));
        return attrs;
    }
}
