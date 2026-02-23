namespace PyVista.Core;

/// <summary>
/// Abstract base class for non-pointset grids (rectilinear grids, image data).
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.Grid</c> class.
/// It provides common functionality shared by grid types that define their
/// geometry implicitly through dimensions and coordinates rather than
/// explicit point arrays.
/// </para>
/// </summary>
public abstract class Grid : DataSet
{
    private int _dimX = 1;
    private int _dimY = 1;
    private int _dimZ = 1;

    /// <summary>
    /// Initializes a new instance of the <see cref="Grid"/> class.
    /// </summary>
    protected Grid()
    {
    }

    /// <summary>
    /// Gets or sets the grid dimensions as a tuple (NX, NY, NZ).
    /// <para>
    /// These are effectively the number of points along each of the
    /// three dataset axes.
    /// </para>
    /// </summary>
    public virtual (int NX, int NY, int NZ) Dimensions
    {
        get => (_dimX, _dimY, _dimZ);
        set
        {
            _dimX = value.NX;
            _dimY = value.NY;
            _dimZ = value.NZ;
        }
    }

    /// <inheritdoc />
    protected override List<(string Name, string Value)> GetAttributes()
    {
        var attrs = base.GetAttributes();
        attrs.Add(("Dimensions", $"{_dimX}, {_dimY}, {_dimZ}"));
        return attrs;
    }

    /// <inheritdoc />
    public override void DeepCopy(DataObject source)
    {
        base.DeepCopy(source);
        if (source is Grid g)
        {
            _dimX = g._dimX;
            _dimY = g._dimY;
            _dimZ = g._dimZ;
        }
    }

    /// <inheritdoc />
    public override void ShallowCopy(DataObject source)
    {
        base.ShallowCopy(source);
        if (source is Grid g)
        {
            _dimX = g._dimX;
            _dimY = g._dimY;
            _dimZ = g._dimZ;
        }
    }
}

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

    // ---------------------------------------------------------------
    //  Meshgrid
    // ---------------------------------------------------------------

    /// <summary>
    /// Returns a meshgrid of the X, Y, and Z coordinates.
    /// <para>
    /// Produces three flat arrays corresponding to the full meshgrid expansion
    /// of the coordinate arrays in IJK (Fortran column-major) order.
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
        var (nx, ny, nz) = Dimensions;
        var pts = Points;
        var sg = new StructuredGrid(pts, nx, ny, nz, deep: false);

        CopyDataArraysTo(sg);
        return sg;
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
