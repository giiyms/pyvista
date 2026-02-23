using PyVista.Core.Cells;

using CT = PyVista.Core.Cells.CellType;

namespace PyVista.Core;

/// <summary>
/// Dataset with variable spacing in the three coordinate directions.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.RectilinearGrid</c> class.
/// A rectilinear grid has axis-aligned cells whose spacing can vary along each
/// coordinate direction. The grid geometry is defined by three arrays of coordinates
/// (one per axis) rather than explicit point locations.
/// </para>
/// </summary>
public class RectilinearGrid : Grid
{
    private double[] _xCoords = [0.0];
    private double[] _yCoords = [0.0];
    private double[] _zCoords = [0.0];

    // ---------------------------------------------------------------
    //  Constructors
    // ---------------------------------------------------------------

    /// <summary>
    /// Initializes a new instance of the <see cref="RectilinearGrid"/> class (empty).
    /// </summary>
    public RectilinearGrid()
    {
        UpdateDimensions();
    }

    /// <summary>
    /// Initializes a new instance of the <see cref="RectilinearGrid"/> class
    /// from X coordinate arrays only.
    /// </summary>
    /// <param name="x">Coordinates of the points in the X direction.</param>
    /// <param name="checkDuplicates">
    /// When <c>true</c>, an error is raised if there are duplicate values.
    /// </param>
    public RectilinearGrid(double[] x, bool checkDuplicates = false)
    {
        ArgumentNullException.ThrowIfNull(x);
        FromArrays(x, null, null, checkDuplicates);
    }

    /// <summary>
    /// Initializes a new instance of the <see cref="RectilinearGrid"/> class
    /// from X and Y coordinate arrays.
    /// </summary>
    /// <param name="x">Coordinates of the points in the X direction.</param>
    /// <param name="y">Coordinates of the points in the Y direction.</param>
    /// <param name="checkDuplicates">
    /// When <c>true</c>, an error is raised if there are duplicate values.
    /// </param>
    public RectilinearGrid(double[] x, double[] y, bool checkDuplicates = false)
    {
        ArgumentNullException.ThrowIfNull(x);
        ArgumentNullException.ThrowIfNull(y);
        FromArrays(x, y, null, checkDuplicates);
    }

    /// <summary>
    /// Initializes a new instance of the <see cref="RectilinearGrid"/> class
    /// from X, Y, and Z coordinate arrays.
    /// </summary>
    /// <param name="x">Coordinates of the points in the X direction.</param>
    /// <param name="y">Coordinates of the points in the Y direction.</param>
    /// <param name="z">Coordinates of the points in the Z direction.</param>
    /// <param name="checkDuplicates">
    /// When <c>true</c>, an error is raised if there are duplicate values.
    /// </param>
    public RectilinearGrid(double[] x, double[] y, double[] z, bool checkDuplicates = false)
    {
        ArgumentNullException.ThrowIfNull(x);
        ArgumentNullException.ThrowIfNull(y);
        ArgumentNullException.ThrowIfNull(z);
        FromArrays(x, y, z, checkDuplicates);
    }

    // ---------------------------------------------------------------
    //  Dimensions (read-only; derived from coordinate arrays)
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets the grid dimensions.
    /// <para>
    /// The dimensions of a <see cref="RectilinearGrid"/> are implicitly defined
    /// by the lengths of the <see cref="X"/>, <see cref="Y"/>, and <see cref="Z"/>
    /// coordinate arrays and cannot be set directly.
    /// </para>
    /// </summary>
    /// <exception cref="InvalidOperationException">
    /// Thrown when attempting to set dimensions directly.
    /// </exception>
    public override (int NX, int NY, int NZ) Dimensions
    {
        get => base.Dimensions;
        set => throw new InvalidOperationException(
            "The dimensions of a RectilinearGrid are implicitly defined and cannot be set.");
    }

    // ---------------------------------------------------------------
    //  X / Y / Z coordinate arrays
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets or sets the coordinates along the X direction.
    /// <para>
    /// Setting this property updates the grid <see cref="Dimensions"/> automatically.
    /// </para>
    /// </summary>
    public double[] X
    {
        get => (double[])_xCoords.Clone();
        set
        {
            ArgumentNullException.ThrowIfNull(value);
            _xCoords = (double[])value.Clone();
            UpdateDimensions();
        }
    }

    /// <summary>
    /// Gets or sets the coordinates along the Y direction.
    /// <para>
    /// Setting this property updates the grid <see cref="Dimensions"/> automatically.
    /// </para>
    /// </summary>
    public double[] Y
    {
        get => (double[])_yCoords.Clone();
        set
        {
            ArgumentNullException.ThrowIfNull(value);
            _yCoords = (double[])value.Clone();
            UpdateDimensions();
        }
    }

    /// <summary>
    /// Gets or sets the coordinates along the Z direction.
    /// <para>
    /// Setting this property updates the grid <see cref="Dimensions"/> automatically.
    /// </para>
    /// </summary>
    public double[] Z
    {
        get => (double[])_zCoords.Clone();
        set
        {
            ArgumentNullException.ThrowIfNull(value);
            _zCoords = (double[])value.Clone();
            UpdateDimensions();
        }
    }

    // ---------------------------------------------------------------
    //  Points (computed from coordinate arrays)
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets a copy of the points as a flat row-major (N×3) array built from the
    /// meshgrid of the X, Y, and Z coordinates.
    /// <para>
    /// Points of a <see cref="RectilinearGrid"/> cannot be set directly.
    /// Set point coordinates with <see cref="X"/>, <see cref="Y"/>, or <see cref="Z"/>.
    /// </para>
    /// </summary>
    /// <exception cref="InvalidOperationException">
    /// Thrown when attempting to set points directly.
    /// </exception>
    public new double[] Points
    {
        get
        {
            int nx = _xCoords.Length;
            int ny = _yCoords.Length;
            int nz = _zCoords.Length;
            int nPoints = nx * ny * nz;
            var pts = new double[nPoints * 3];

            // Iterate in Fortran (column-major) order: x varies slowest, z varies fastest
            int idx = 0;
            for (int iz = 0; iz < nz; iz++)
            {
                for (int iy = 0; iy < ny; iy++)
                {
                    for (int ix = 0; ix < nx; ix++)
                    {
                        pts[idx++] = _xCoords[ix];
                        pts[idx++] = _yCoords[iy];
                        pts[idx++] = _zCoords[iz];
                    }
                }
            }

            return pts;
        }
    }

    /// <summary>
    /// Gets the number of points in this rectilinear grid.
    /// <para>
    /// The number of points is the product of the lengths of the X, Y, and Z
    /// coordinate arrays.
    /// </para>
    /// </summary>
    public new int NPoints => _xCoords.Length * _yCoords.Length * _zCoords.Length;

    // ---------------------------------------------------------------
    //  Meshgrid
    // ---------------------------------------------------------------

    /// <summary>
    /// Returns a meshgrid of the X, Y, and Z coordinates.
    /// <para>
    /// Produces three flat arrays corresponding to the full meshgrid expansion
    /// of the coordinate arrays in IJK (Fortran column-major) order.
    /// This is the C# equivalent of calling <c>numpy.meshgrid(x, y, z, indexing='ij')</c>
    /// and then flattening with Fortran order.
    /// </para>
    /// </summary>
    /// <returns>
    /// A tuple of three arrays (XX, YY, ZZ) where each has length NX×NY×NZ.
    /// </returns>
    public (double[] XX, double[] YY, double[] ZZ) Meshgrid
    {
        get
        {
            int nx = _xCoords.Length;
            int ny = _yCoords.Length;
            int nz = _zCoords.Length;
            int total = nx * ny * nz;
            var xx = new double[total];
            var yy = new double[total];
            var zz = new double[total];

            int idx = 0;
            for (int iz = 0; iz < nz; iz++)
            {
                for (int iy = 0; iy < ny; iy++)
                {
                    for (int ix = 0; ix < nx; ix++)
                    {
                        xx[idx] = _xCoords[ix];
                        yy[idx] = _yCoords[iy];
                        zz[idx] = _zCoords[iz];
                        idx++;
                    }
                }
            }

            return (xx, yy, zz);
        }
    }

    // ---------------------------------------------------------------
    //  Bounds
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets the axis-aligned bounding box of this rectilinear grid.
    /// <para>
    /// Bounds are computed from the min/max values of each coordinate array.
    /// </para>
    /// </summary>
    public new BoundsTuple Bounds
    {
        get
        {
            double xMin = MinOf(_xCoords);
            double xMax = MaxOf(_xCoords);
            double yMin = MinOf(_yCoords);
            double yMax = MaxOf(_yCoords);
            double zMin = MinOf(_zCoords);
            double zMax = MaxOf(_zCoords);
            return new BoundsTuple(xMin, xMax, yMin, yMax, zMin, zMax);
        }
    }

    // ---------------------------------------------------------------
    //  Implicit cell count
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets the number of cells implied by the grid dimensions.
    /// <para>
    /// For a rectilinear grid, the number of cells is
    /// <c>max(1, NX-1) × max(1, NY-1) × max(1, NZ-1)</c>.
    /// </para>
    /// </summary>
    public int NCellsFromDimensions
    {
        get
        {
            int cx = Math.Max(1, _xCoords.Length - 1);
            int cy = Math.Max(1, _yCoords.Length - 1);
            int cz = Math.Max(1, _zCoords.Length - 1);
            return cx * cy * cz;
        }
    }

    // ---------------------------------------------------------------
    //  Cast to StructuredGrid
    // ---------------------------------------------------------------

    /// <summary>
    /// Casts this rectilinear grid to a <see cref="StructuredGrid"/>.
    /// <para>
    /// The resulting structured grid has explicit point coordinates computed from
    /// the meshgrid of the X, Y, and Z coordinate arrays. Point data, cell data,
    /// and field data are copied to the new grid.
    /// </para>
    /// </summary>
    /// <returns>A new <see cref="StructuredGrid"/> with the same geometry.</returns>
    public StructuredGrid CastToStructuredGrid()
    {
        var (nx, ny, nz) = base.Dimensions;
        var pts = Points;
        var sg = new StructuredGrid(pts, nx, ny, nz, deep: false);

        CopyDataArraysTo(sg);
        return sg;
    }

    // ---------------------------------------------------------------
    //  Cast to UnstructuredGrid
    // ---------------------------------------------------------------

    /// <summary>
    /// Casts this rectilinear grid to an <see cref="UnstructuredGrid"/>.
    /// <para>
    /// All cells of the rectilinear grid are converted to hexahedral cells.
    /// Each cell is defined by the eight corner points of its axis-aligned
    /// bounding box. Point data, cell data, and field data are copied to
    /// the new grid.
    /// </para>
    /// <para>
    /// This is the C# equivalent of the Python
    /// <c>RectilinearGrid.cast_to_unstructured_grid()</c> method.
    /// </para>
    /// </summary>
    /// <returns>A new <see cref="UnstructuredGrid"/> with hexahedral cells.</returns>
    public new UnstructuredGrid CastToUnstructuredGrid()
    {
        int nx = _xCoords.Length;
        int ny = _yCoords.Length;
        int nz = _zCoords.Length;

        // Number of cells in each direction
        int cx = Math.Max(1, nx - 1);
        int cy = Math.Max(1, ny - 1);
        int cz = Math.Max(1, nz - 1);
        int nCells = cx * cy * cz;

        // Build flat points array from the meshgrid
        var pts = Points;

        // Build hexahedral cell connectivity
        // Each hex cell has 8 corner points; legacy format: [8, p0, p1, ..., p7]
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
                    // Point indices for the 8 corners of the hexahedron
                    // using the same IJK ordering as the points array
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

        var ugrid = new UnstructuredGrid(cells, cellTypes, pts, deep: false);

        CopyDataArraysTo(ugrid);
        return ugrid;
    }

    // ---------------------------------------------------------------
    //  Copy helpers
    // ---------------------------------------------------------------

    /// <inheritdoc />
    public override void DeepCopy(DataObject source)
    {
        base.DeepCopy(source);
        if (source is RectilinearGrid rg)
        {
            _xCoords = (double[])rg._xCoords.Clone();
            _yCoords = (double[])rg._yCoords.Clone();
            _zCoords = (double[])rg._zCoords.Clone();
        }
    }

    /// <inheritdoc />
    public override void ShallowCopy(DataObject source)
    {
        base.ShallowCopy(source);
        if (source is RectilinearGrid rg)
        {
            _xCoords = rg._xCoords;
            _yCoords = rg._yCoords;
            _zCoords = rg._zCoords;
        }
    }

    /// <inheritdoc />
    protected override List<(string Name, string Value)> GetAttributes()
    {
        var attrs = base.GetAttributes();
        var b = Bounds;
        attrs.Add(("Bounds", $"{b.XMin:E1}, {b.XMax:E1}, {b.YMin:E1}, {b.YMax:E1}, {b.ZMin:E1}, {b.ZMax:E1}"));
        return attrs;
    }

    // ---------------------------------------------------------------
    //  Private helpers
    // ---------------------------------------------------------------

    /// <summary>
    /// Initializes the coordinate arrays from the given inputs.
    /// </summary>
    private void FromArrays(double[] x, double[]? y, double[]? z, bool checkDuplicates)
    {
        if (checkDuplicates)
        {
            ThrowIfHasDuplicates(x, nameof(x));
        }

        _xCoords = (double[])x.Clone();

        if (y is not null)
        {
            if (checkDuplicates)
            {
                ThrowIfHasDuplicates(y, nameof(y));
            }

            _yCoords = (double[])y.Clone();
        }

        if (z is not null)
        {
            if (checkDuplicates)
            {
                ThrowIfHasDuplicates(z, nameof(z));
            }

            _zCoords = (double[])z.Clone();
        }

        UpdateDimensions();
    }

    /// <summary>
    /// Updates the base dimensions from the current coordinate array lengths.
    /// </summary>
    private void UpdateDimensions()
    {
        base.Dimensions = (_xCoords.Length, _yCoords.Length, _zCoords.Length);
    }

    /// <summary>
    /// Copies point data, cell data, and field data to a target dataset.
    /// </summary>
    private void CopyDataArraysTo(DataSet target)
    {
        foreach (var name in PointData.Keys)
        {
            target.PointData.SetArray((double[])PointData[name].Clone(), name);
        }

        foreach (var name in CellData.Keys)
        {
            target.CellData.SetArray((double[])CellData[name].Clone(), name);
        }

        foreach (var name in FieldData.Keys)
        {
            target.FieldData.SetArray((double[])FieldData[name].Clone(), name);
        }
    }

    /// <summary>
    /// Computes the flat point index for the given IJK coordinates using the
    /// Fortran (column-major) ordering: <c>ix + nx * (iy + ny * iz)</c>.
    /// </summary>
    private static int PointIndex(int ix, int iy, int iz, int nx, int ny)
    {
        return ix + nx * (iy + ny * iz);
    }

    /// <summary>
    /// Throws <see cref="ArgumentException"/> if the array contains duplicate values.
    /// </summary>
    private static void ThrowIfHasDuplicates(double[] array, string paramName)
    {
        var seen = new HashSet<double>();
        foreach (double v in array)
        {
            if (!seen.Add(v))
            {
                throw new ArgumentException($"Array '{paramName}' contains duplicate values.", paramName);
            }
        }
    }

    /// <summary>Returns the minimum value in the array.</summary>
    private static double MinOf(double[] arr)
    {
        double min = arr[0];
        for (int i = 1; i < arr.Length; i++)
        {
            if (arr[i] < min) min = arr[i];
        }

        return min;
    }

    /// <summary>Returns the maximum value in the array.</summary>
    private static double MaxOf(double[] arr)
    {
        double max = arr[0];
        for (int i = 1; i < arr.Length; i++)
        {
            if (arr[i] > max) max = arr[i];
        }

        return max;
    }
}
