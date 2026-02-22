using System.Globalization;

using PyVista.Core.Cells;

using CT = PyVista.Core.Cells.CellType;

namespace PyVista.Core;

/// <summary>
/// Dataset used for arbitrary combinations of all possible cell types.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.UnstructuredGrid</c> class.
/// An unstructured grid stores cells of arbitrary type, defined by separate
/// cell connectivity (<see cref="Cells"/>), cell type (<see cref="CellTypes"/>),
/// and offset (<see cref="Offset"/>) arrays.
/// </para>
/// </summary>
public class UnstructuredGrid : PointSet
{
    private CellArray _cellArray = new();
    private byte[] _cellTypes = Array.Empty<byte>();

    /// <summary>
    /// Initializes a new instance of the <see cref="UnstructuredGrid"/> class (empty).
    /// </summary>
    public UnstructuredGrid()
    {
    }

    /// <summary>
    /// Initializes a new instance of the <see cref="UnstructuredGrid"/> class from
    /// cell connectivity, cell types, and point arrays.
    /// </summary>
    /// <param name="cells">
    /// Padded cell connectivity array in legacy format:
    /// <c>[n0, p0_0, …, p0_n, n1, p1_0, …, p1_n, …]</c>.
    /// </param>
    /// <param name="cellTypes">
    /// Array of cell type codes (one per cell). Values correspond to
    /// <see cref="CellType"/> enum values.
    /// </param>
    /// <param name="points">
    /// Flat row-major point coordinates (length must be divisible by 3).
    /// </param>
    /// <param name="deep">When <c>true</c>, copies all input arrays.</param>
    /// <exception cref="ArgumentNullException">
    /// Thrown when any parameter is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when the number of cell types does not match the number of cells.
    /// </exception>
    public UnstructuredGrid(int[] cells, byte[] cellTypes, double[] points, bool deep = false)
        : base(points, deep)
    {
        ArgumentNullException.ThrowIfNull(cells);
        ArgumentNullException.ThrowIfNull(cellTypes);

        _cellArray = new CellArray(deep ? (int[])cells.Clone() : cells);
        _cellTypes = deep ? (byte[])cellTypes.Clone() : cellTypes;

        CheckForConsistency();
    }

    /// <summary>
    /// Initializes a new instance of the <see cref="UnstructuredGrid"/> class from
    /// cell connectivity, cell types (as <see cref="CellType"/> array), and points.
    /// </summary>
    /// <param name="cells">Padded cell connectivity array in legacy format.</param>
    /// <param name="cellTypes">Array of <see cref="CT"/> values (one per cell).</param>
    /// <param name="points">Flat row-major point coordinates.</param>
    /// <param name="deep">When <c>true</c>, copies all input arrays.</param>
    public UnstructuredGrid(int[] cells, CT[] cellTypes, double[] points, bool deep = false)
        : base(points, deep)
    {
        ArgumentNullException.ThrowIfNull(cells);
        ArgumentNullException.ThrowIfNull(cellTypes);

        _cellArray = new CellArray(deep ? (int[])cells.Clone() : cells);
        _cellTypes = new byte[cellTypes.Length];
        for (int i = 0; i < cellTypes.Length; i++)
        {
            _cellTypes[i] = (byte)cellTypes[i];
        }

        CheckForConsistency();
    }

    // ---------------------------------------------------------------
    //  Cells property
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets or sets the cell data in legacy padded format.
    /// <para>
    /// The format is:
    /// <c>[n0, p0_0, …, p0_n, n1, p1_0, …, p1_n, …]</c>
    /// where <c>n0</c> is the number of points in cell 0.
    /// </para>
    /// </summary>
    public int[] Cells
    {
        get => _cellArray.Cells;
        set
        {
            ArgumentNullException.ThrowIfNull(value);
            _cellArray = new CellArray(value);
        }
    }

    // ---------------------------------------------------------------
    //  Cell types
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets or sets the cell types array.
    /// <para>
    /// Each element is a <see cref="byte"/> corresponding to a <see cref="CellType"/> value.
    /// </para>
    /// </summary>
    public byte[] CellTypes
    {
        get => (byte[])_cellTypes.Clone();
        set
        {
            ArgumentNullException.ThrowIfNull(value);
            _cellTypes = (byte[])value.Clone();
        }
    }

    /// <summary>
    /// Gets the cell types as a <see cref="CT"/> array.
    /// </summary>
    public CT[] CellTypesEnum
    {
        get
        {
            var result = new CT[_cellTypes.Length];
            for (int i = 0; i < _cellTypes.Length; i++)
            {
                result[i] = (CT)_cellTypes[i];
            }

            return result;
        }
    }

    // ---------------------------------------------------------------
    //  Offset
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets the cell locations (offset) array.
    /// <para>
    /// This array has length <c>NCells + 1</c> and indicates the start of each cell
    /// in the <see cref="CellConnectivity"/> array.
    /// </para>
    /// </summary>
    public int[] Offset => _cellArray.OffsetArray;

    // ---------------------------------------------------------------
    //  Cell connectivity
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets the cell connectivity array (cells without padding).
    /// <para>
    /// This is effectively <see cref="Cells"/> without the interleaved cell-size counts.
    /// </para>
    /// </summary>
    public int[] CellConnectivity => _cellArray.ConnectivityArray;

    // ---------------------------------------------------------------
    //  Cell count
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets the number of cells.
    /// </summary>
    public int NCellsTotal => _cellArray.NCells;

    // ---------------------------------------------------------------
    //  Cells dict
    // ---------------------------------------------------------------

    /// <summary>
    /// Returns a dictionary that maps each distinct cell type to a 2-D array
    /// of cell connectivity (point indices). Only works with fixed-size cell types.
    /// </summary>
    /// <returns>
    /// A dictionary mapping <see cref="byte"/> cell type codes to 2-D arrays of
    /// shape <c>(nCells, nPointsPerCell)</c>.
    /// </returns>
    public Dictionary<byte, List<int[]>> CellsDict
    {
        get
        {
            var result = new Dictionary<byte, List<int[]>>();
            int nCells = _cellArray.NCells;
            var offsets = _cellArray.OffsetArray;
            var conn = _cellArray.ConnectivityArray;

            for (int i = 0; i < nCells; i++)
            {
                byte ct = _cellTypes[i];
                int start = offsets[i];
                int end = offsets[i + 1];
                int len = end - start;
                var cellConn = new int[len];
                Array.Copy(conn, start, cellConn, 0, len);

                if (!result.TryGetValue(ct, out var list))
                {
                    list = new List<int[]>();
                    result[ct] = list;
                }

                list.Add(cellConn);
            }

            return result;
        }
    }

    // ---------------------------------------------------------------
    //  Linear copy
    // ---------------------------------------------------------------

    /// <summary>
    /// Returns a copy of the unstructured grid containing only linear cells.
    /// <para>
    /// Converts quadratic cell types to their linear equivalents by updating the
    /// cell type codes. Point coordinates and cell connectivity are preserved.
    /// </para>
    /// </summary>
    /// <param name="deep">When <c>true</c>, deep copies the point array.</param>
    /// <returns>A new <see cref="UnstructuredGrid"/> with linear cell types.</returns>
    public UnstructuredGrid LinearCopy(bool deep = false)
    {
        var copy = (UnstructuredGrid)Copy(deep: deep);

        for (int i = 0; i < copy._cellTypes.Length; i++)
        {
            copy._cellTypes[i] = (CT)copy._cellTypes[i] switch
            {
                CT.QuadraticTriangle => (byte)CT.Triangle,
                CT.QuadraticQuad => (byte)CT.Quad,
                CT.QuadraticTetra => (byte)CT.Tetra,
                CT.QuadraticPyramid => (byte)CT.Pyramid,
                CT.QuadraticWedge => (byte)CT.Wedge,
                CT.QuadraticHexahedron => (byte)CT.Hexahedron,
                _ => copy._cellTypes[i],
            };
        }

        return copy;
    }

    // ---------------------------------------------------------------
    //  Copy helpers
    // ---------------------------------------------------------------

    /// <inheritdoc />
    public override void DeepCopy(DataObject source)
    {
        base.DeepCopy(source);
        if (source is UnstructuredGrid ug)
        {
            _cellArray = new CellArray(ug._cellArray.Cells);
            _cellTypes = (byte[])ug._cellTypes.Clone();
        }
    }

    /// <inheritdoc />
    public override void ShallowCopy(DataObject source)
    {
        base.ShallowCopy(source);
        if (source is UnstructuredGrid ug)
        {
            _cellArray = ug._cellArray;
            _cellTypes = ug._cellTypes;
        }
    }

    /// <inheritdoc />
    protected override List<(string Name, string Value)> GetAttributes()
    {
        var attrs = base.GetAttributes();
        return attrs;
    }

    // ---------------------------------------------------------------
    //  Validation
    // ---------------------------------------------------------------

    /// <summary>
    /// Checks that the cell types and offset arrays are consistent with the
    /// number of cells.
    /// </summary>
    /// <exception cref="ArgumentException">
    /// Thrown when the cell types array size does not match the cell count,
    /// or when the offset array size is inconsistent.
    /// </exception>
    private void CheckForConsistency()
    {
        int nCells = _cellArray.NCells;

        if (nCells != _cellTypes.Length)
        {
            throw new ArgumentException(
                $"Number of cell types ({_cellTypes.Length}) must match the number of cells ({nCells}).");
        }

        var offsets = _cellArray.OffsetArray;
        if (nCells != offsets.Length - 1)
        {
            throw new ArgumentException(
                $"Size of the offset ({offsets.Length}) must be one greater than the number of cells ({nCells}).");
        }
    }
}
