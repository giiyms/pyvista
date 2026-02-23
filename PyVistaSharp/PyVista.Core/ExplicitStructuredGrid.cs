using System.Globalization;

using PyVista.Core.Cells;

using CT = PyVista.Core.Cells.CellType;

namespace PyVista.Core;

/// <summary>
/// Explicit structured grid dataset.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.ExplicitStructuredGrid</c> class.
/// An explicit structured grid is a point-based dataset that combines structured grid
/// topology (I, J, K indexing) with explicit cell connectivity. Unlike
/// <see cref="StructuredGrid"/>, which has purely implicit connectivity, this class
/// stores the cell-to-point mapping explicitly, allowing cells with degenerate or
/// irregular geometry.
/// </para>
/// <para>
/// Cells are expected to be hexahedra (8 points per cell). The grid supports hiding
/// and showing individual cells via a ghost cell array, computing face connectivity
/// flags, and converting between linear cell IDs and structured (i, j, k) coordinates.
/// </para>
/// </summary>
public class ExplicitStructuredGrid : PointSet
{
    /// <summary>Ghost cell flag value indicating a hidden cell.</summary>
    private const int HiddenCellFlag = 2;

    /// <summary>Six-connected topological neighborhood offsets.</summary>
    private static readonly (int di, int dj, int dk)[] NeighborOffsets =
    [
        (-1, 0, 0), (1, 0, 0),
        (0, -1, 0), (0, 1, 0),
        (0, 0, -1), (0, 0, 1),
    ];

    /// <summary>Name of the VTK ghost cell array.</summary>
    private const string GhostCellArrayName = "vtkGhostType";

    /// <summary>Name of the connectivity flags cell data array.</summary>
    private const string ConnectivityFlagsArrayName = "ConnectivityFlags";

    /// <summary>Name of the number-of-connections cell data array.</summary>
    private const string NumberOfConnectionsArrayName = "number_of_connections";

    /// <summary>Name of the visibility cell data array.</summary>
    private const string VisibilityArrayName = "Visibility";

    private CellArray _cellArray = new();
    private int _dimX = 1;
    private int _dimY = 1;
    private int _dimZ = 1;

    // ---------------------------------------------------------------
    //  Constructors
    // ---------------------------------------------------------------

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
    /// <para>
    /// Dimensions define the number of sampling points in the I, J and K directions.
    /// The number of cells in each direction is one less than the corresponding
    /// dimension value.
    /// </para>
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
    /// Gets the total number of cells in the grid.
    /// </summary>
    public int NCellsTotal => _cellArray.NCells;

    // ---------------------------------------------------------------
    //  Visible bounds
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets the bounding box of the visible cells.
    /// <para>
    /// Different from <see cref="DataSet.Bounds"/> which returns the bounding box of the
    /// complete grid, this property returns the bounding box considering only cells
    /// whose ghost cell flag is not <c>HIDDENCELL</c>. When no ghost cell array exists
    /// or all cells are visible, the full <see cref="DataSet.Bounds"/> is returned.
    /// </para>
    /// </summary>
    public BoundsTuple VisibleBounds
    {
        get
        {
            if (!CellData.ContainsKey(GhostCellArrayName))
            {
                return Bounds;
            }

            var ghost = CellData[GhostCellArrayName];
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
                if ((int)ghost[c] == HiddenCellFlag) continue;

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
    //  Cast to UnstructuredGrid
    // ---------------------------------------------------------------

    /// <summary>
    /// Casts this explicit structured grid to an <see cref="UnstructuredGrid"/>.
    /// <para>
    /// All cells are treated as hexahedra. The <c>BLOCK_I</c>, <c>BLOCK_J</c>,
    /// and <c>BLOCK_K</c> cell arrays are added to enable restoring the explicit
    /// structured grid from the unstructured representation.
    /// </para>
    /// </summary>
    /// <returns>A new <see cref="UnstructuredGrid"/>.</returns>
    public new UnstructuredGrid CastToUnstructuredGrid()
    {
        int nCells = _cellArray.NCells;
        var cellTypesArr = new byte[nCells];
        Array.Fill(cellTypesArr, (byte)CT.Hexahedron);

        var ugrid = new UnstructuredGrid(_cellArray.Cells, cellTypesArr, (double[])Points.Clone(), deep: false);

        int ni = Math.Max(1, _dimX - 1);
        int nj = Math.Max(1, _dimY - 1);

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
    /// Hides specific cells by setting the ghost cell array to <c>HIDDENCELL</c>.
    /// <para>
    /// Hidden cells are excluded from rendering and from the
    /// <see cref="VisibleBounds"/> computation. The original cell data is preserved;
    /// hiding is reversible via <see cref="ShowCells"/>.
    /// </para>
    /// </summary>
    /// <param name="indices">
    /// Cell indices to hide. Indices outside the valid range are silently ignored.
    /// </param>
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
                ghostCells[idx] = HiddenCellFlag;
            }
        }

        target.CellData.SetArray(ghostCells, GhostCellArrayName);
        return target;
    }

    /// <summary>
    /// Shows hidden cells by clearing the <c>HIDDENCELL</c> ghost cell flag.
    /// <para>
    /// Cells whose ghost flag equals <c>HIDDENCELL</c> are reset to <c>0</c> (visible).
    /// Other ghost flag values are left untouched.
    /// </para>
    /// </summary>
    /// <param name="inplace">When <c>true</c>, modifies this grid; otherwise returns a copy.</param>
    /// <returns>The grid with all previously hidden cells made visible.</returns>
    public ExplicitStructuredGrid ShowCells(bool inplace = false)
    {
        ExplicitStructuredGrid target = inplace ? this : (ExplicitStructuredGrid)Copy(deep: true);
        if (target.CellData.ContainsKey(GhostCellArrayName))
        {
            var array = target.CellData[GhostCellArrayName];
            for (int i = 0; i < array.Length; i++)
            {
                if ((int)array[i] == HiddenCellFlag)
                {
                    array[i] = 0;
                }
            }
        }

        return target;
    }

    // ---------------------------------------------------------------
    //  Cell ID / coords conversion
    // ---------------------------------------------------------------

    /// <summary>
    /// Returns the linear cell ID for the given structured coordinates (i, j, k).
    /// <para>
    /// The mapping follows column-major (Fortran) ordering:
    /// <c>id = i + ni * (j + nj * k)</c>, where <c>ni = NI - 1</c> and
    /// <c>nj = NJ - 1</c>.
    /// </para>
    /// </summary>
    /// <param name="i">I-index (must be in <c>[0, NI-2]</c>).</param>
    /// <param name="j">J-index (must be in <c>[0, NJ-2]</c>).</param>
    /// <param name="k">K-index (must be in <c>[0, NK-2]</c>).</param>
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
    /// Returns the structured coordinates (i, j, k) for the given linear cell ID.
    /// </summary>
    /// <param name="cellId">The cell ID.</param>
    /// <returns>
    /// A tuple of (I, J, K) structured coordinates, or <c>null</c> if the ID is invalid.
    /// </returns>
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
    //  Neighbors
    // ---------------------------------------------------------------

    /// <summary>
    /// Returns the indices of topological neighboring cells.
    /// <para>
    /// For a given cell, the six topological neighbors are at
    /// <c>(i±1, j, k)</c>, <c>(i, j±1, k)</c>, and <c>(i, j, k±1)</c>.
    /// Only neighbors that lie within the grid extent are returned.
    /// </para>
    /// </summary>
    /// <param name="cellId">The cell ID whose neighbors are requested.</param>
    /// <returns>A sorted list of neighboring cell IDs.</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// Thrown when <paramref name="cellId"/> is outside the valid cell range.
    /// </exception>
    public List<int> Neighbors(int cellId)
    {
        var coords = CellCoords(cellId);
        if (coords is null)
        {
            throw new ArgumentOutOfRangeException(nameof(cellId), cellId,
                "Cell ID is outside the valid range.");
        }

        var (ci, cj, ck) = coords.Value;
        var indices = new List<int>(6);

        foreach (var (di, dj, dk) in NeighborOffsets)
        {
            int id = CellId(ci + di, cj + dj, ck + dk);
            if (id >= 0)
            {
                indices.Add(id);
            }
        }

        indices.Sort();
        return indices;
    }

    /// <summary>
    /// Returns the indices of neighboring cells for multiple cell IDs.
    /// <para>
    /// The returned list contains all unique neighbors of all supplied cells,
    /// sorted in ascending order. A cell that appears as a neighbor of more
    /// than one input cell is included only once.
    /// </para>
    /// </summary>
    /// <param name="cellIds">Cell IDs whose neighbors are requested.</param>
    /// <returns>A sorted list of unique neighboring cell IDs.</returns>
    public List<int> Neighbors(int[] cellIds)
    {
        ArgumentNullException.ThrowIfNull(cellIds);
        var unique = new SortedSet<int>();
        foreach (int id in cellIds)
        {
            foreach (int neighbor in Neighbors(id))
            {
                unique.Add(neighbor);
            }
        }

        return new List<int>(unique);
    }

    // ---------------------------------------------------------------
    //  Compute connectivity / connections / visibility
    // ---------------------------------------------------------------

    /// <summary>
    /// Computes the faces connectivity flags array.
    /// <para>
    /// Each value is interpreted as a 6-bit binary number indicating face connectivity
    /// with topological neighbors. The bits represent, from least-significant to
    /// most-significant: <c>-X</c>, <c>+X</c>, <c>-Y</c>, <c>+Y</c>, <c>-Z</c>,
    /// <c>+Z</c>. For example, a value of <c>27</c> (<c>011011</c>) means the cell
    /// shares faces with neighbors <c>(-1,0,0)</c>, <c>(+1,0,0)</c>, <c>(0,-1,0)</c>,
    /// and <c>(0,0,+1)</c>.
    /// </para>
    /// <para>
    /// The result is stored in the <c>"ConnectivityFlags"</c> cell data array.
    /// </para>
    /// </summary>
    /// <param name="inplace">When <c>true</c>, modifies this grid; otherwise returns a copy.</param>
    /// <returns>The grid with connectivity flags computed.</returns>
    /// <seealso cref="ComputeConnections"/>
    public ExplicitStructuredGrid ComputeConnectivity(bool inplace = false)
    {
        ExplicitStructuredGrid target = inplace ? this : (ExplicitStructuredGrid)Copy(deep: true);
        int nCells = target._cellArray.NCells;
        var flags = new double[nCells];

        int ni = _dimX - 1;
        int nj = _dimY - 1;
        int nk = _dimZ - 1;

        for (int c = 0; c < nCells; c++)
        {
            var coords = target.CellCoords(c);
            if (coords is null) continue;

            var (ci, cj, ck) = coords.Value;
            int flag = 0;

            if (ci > 0) flag |= 1;        // -X neighbor
            if (ci < ni - 1) flag |= 2;   // +X neighbor
            if (cj > 0) flag |= 4;        // -Y neighbor
            if (cj < nj - 1) flag |= 8;   // +Y neighbor
            if (ck > 0) flag |= 16;       // -Z neighbor
            if (ck < nk - 1) flag |= 32;  // +Z neighbor

            flags[c] = flag;
        }

        target.CellData.SetArray(flags, ConnectivityFlagsArrayName);
        return target;
    }

    /// <summary>
    /// Computes an array with the number of connected cell faces.
    /// <para>
    /// For each cell, counts how many of its six faces are shared with topological
    /// neighbors. If the <c>"ConnectivityFlags"</c> array does not already exist,
    /// <see cref="ComputeConnectivity"/> is called first. The result is stored in the
    /// <c>"number_of_connections"</c> cell data array.
    /// </para>
    /// </summary>
    /// <param name="inplace">When <c>true</c>, modifies this grid; otherwise returns a copy.</param>
    /// <returns>The grid with connection counts computed.</returns>
    /// <seealso cref="ComputeConnectivity"/>
    public ExplicitStructuredGrid ComputeConnections(bool inplace = false)
    {
        ExplicitStructuredGrid target = inplace ? this : (ExplicitStructuredGrid)Copy(deep: true);

        double[]? flagsArr = null;
        if (target.CellData.ContainsKey(ConnectivityFlagsArrayName))
        {
            flagsArr = target.CellData[ConnectivityFlagsArrayName];
        }
        else
        {
            var temp = target.ComputeConnectivity(inplace: false);
            flagsArr = temp.CellData[ConnectivityFlagsArrayName];
        }

        int nCells = flagsArr.Length;
        var connections = new double[nCells];
        for (int i = 0; i < nCells; i++)
        {
            connections[i] = PopCount((int)flagsArr[i]);
        }

        target.CellData.SetArray(connections, NumberOfConnectionsArrayName);
        return target;
    }

    /// <summary>
    /// Computes a visibility array for each cell.
    /// <para>
    /// Each cell is assigned a value of <c>1</c> (visible) or <c>0</c> (hidden)
    /// based on the ghost cell array. If no ghost cell array exists, all cells
    /// are marked as visible. The result is stored in the <c>"Visibility"</c>
    /// cell data array.
    /// </para>
    /// </summary>
    /// <param name="inplace">When <c>true</c>, modifies this grid; otherwise returns a copy.</param>
    /// <returns>The grid with the visibility array computed.</returns>
    public ExplicitStructuredGrid ComputeVisibility(bool inplace = false)
    {
        ExplicitStructuredGrid target = inplace ? this : (ExplicitStructuredGrid)Copy(deep: true);
        int nCells = target._cellArray.NCells;
        var visibility = new double[nCells];

        if (target.CellData.ContainsKey(GhostCellArrayName))
        {
            var ghost = target.CellData[GhostCellArrayName];
            for (int i = 0; i < nCells; i++)
            {
                visibility[i] = (int)ghost[i] == HiddenCellFlag ? 0.0 : 1.0;
            }
        }
        else
        {
            Array.Fill(visibility, 1.0);
        }

        target.CellData.SetArray(visibility, VisibilityArrayName);
        return target;
    }

    // ---------------------------------------------------------------
    //  Number of visible / hidden cells
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets the number of visible (non-hidden) cells.
    /// </summary>
    public int NVisibleCells
    {
        get
        {
            if (!CellData.ContainsKey(GhostCellArrayName))
            {
                return _cellArray.NCells;
            }

            var ghost = CellData[GhostCellArrayName];
            int count = 0;
            for (int i = 0; i < ghost.Length; i++)
            {
                if ((int)ghost[i] != HiddenCellFlag) count++;
            }

            return count;
        }
    }

    /// <summary>
    /// Gets the number of hidden cells.
    /// </summary>
    public int NHiddenCells => _cellArray.NCells - NVisibleCells;

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
        attrs.Add(("N Cells (total)", _cellArray.NCells.ToString(CultureInfo.InvariantCulture)));
        attrs.Add(("N Visible Cells", NVisibleCells.ToString(CultureInfo.InvariantCulture)));
        return attrs;
    }

    // ---------------------------------------------------------------
    //  Private helpers
    // ---------------------------------------------------------------

    /// <summary>
    /// Counts the number of set bits (population count) in an integer.
    /// </summary>
    private static int PopCount(int value)
    {
        int count = 0;
        while (value > 0)
        {
            count += value & 1;
            value >>= 1;
        }

        return count;
    }
}
