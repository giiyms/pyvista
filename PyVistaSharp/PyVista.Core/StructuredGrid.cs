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
