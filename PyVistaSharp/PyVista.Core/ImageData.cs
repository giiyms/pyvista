namespace PyVista.Core;

/// <summary>
/// Models datasets with uniform spacing in the three coordinate directions.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.ImageData</c> class.
/// An ImageData defines a regular, axis-aligned grid whose point locations are
/// implicitly determined by its <see cref="Grid.Dimensions"/>, <see cref="Spacing"/>,
/// <see cref="Origin"/>, and <see cref="DirectionMatrix"/>.
/// </para>
/// </summary>
public class ImageData : Grid
{
    private double _originX;
    private double _originY;
    private double _originZ;
    private double _spacingX = 1.0;
    private double _spacingY = 1.0;
    private double _spacingZ = 1.0;
    private int _offsetX;
    private int _offsetY;
    private int _offsetZ;
    private double[] _directionMatrix = [1, 0, 0, 0, 1, 0, 0, 0, 1];

    /// <summary>
    /// Initializes a new instance of the <see cref="ImageData"/> class (empty).
    /// </summary>
    public ImageData()
    {
    }

    /// <summary>
    /// Initializes a new instance of the <see cref="ImageData"/> class with the
    /// given dimensions, spacing, origin, and optional direction matrix and offset.
    /// </summary>
    /// <param name="dimensions">Number of points along each axis (NX, NY, NZ).</param>
    /// <param name="spacing">Spacing between points along each axis. Must be non-negative.</param>
    /// <param name="origin">Origin of the grid (bottom south-west corner).</param>
    /// <param name="directionMatrix">
    /// Optional 3×3 direction matrix stored as a 9-element row-major array.
    /// Defaults to the identity matrix.
    /// </param>
    /// <param name="offset">
    /// Optional index offset for each axis. Defines the minimum extent.
    /// </param>
    public ImageData(
        (int NX, int NY, int NZ) dimensions,
        (double X, double Y, double Z)? spacing = null,
        (double X, double Y, double Z)? origin = null,
        double[]? directionMatrix = null,
        (int X, int Y, int Z)? offset = null)
    {
        Dimensions = dimensions;

        if (origin.HasValue)
        {
            _originX = origin.Value.X;
            _originY = origin.Value.Y;
            _originZ = origin.Value.Z;
        }

        if (spacing.HasValue)
        {
            Spacing = spacing.Value;
        }

        if (directionMatrix is not null)
        {
            DirectionMatrix = directionMatrix;
        }

        if (offset.HasValue)
        {
            Offset = offset.Value;
        }
    }

    // ---------------------------------------------------------------
    //  Origin
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets or sets the origin of the grid (bottom south-west corner).
    /// </summary>
    public (double X, double Y, double Z) Origin
    {
        get => (_originX, _originY, _originZ);
        set
        {
            _originX = value.X;
            _originY = value.Y;
            _originZ = value.Z;
        }
    }

    // ---------------------------------------------------------------
    //  Spacing
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets or sets the spacing between points along each axis.
    /// <para>
    /// Spacing must be non-negative. Negative spacing results in
    /// unexpected behavior and is not allowed.
    /// </para>
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException">
    /// Thrown when any spacing component is negative.
    /// </exception>
    public (double X, double Y, double Z) Spacing
    {
        get => (_spacingX, _spacingY, _spacingZ);
        set
        {
            if (value.X < 0 || value.Y < 0 || value.Z < 0)
            {
                throw new ArgumentOutOfRangeException(
                    nameof(value), "Spacing must be non-negative.");
            }

            _spacingX = value.X;
            _spacingY = value.Y;
            _spacingZ = value.Z;
        }
    }

    // ---------------------------------------------------------------
    //  Extent
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets or sets the extent of the ImageData.
    /// <para>
    /// The extent is six values (XMin, XMax, YMin, YMax, ZMin, ZMax) representing
    /// the first and last point indices for each axis. Setting the extent updates
    /// both the <see cref="Grid.Dimensions"/> and <see cref="Offset"/>.
    /// </para>
    /// </summary>
    /// <exception cref="ArgumentException">
    /// Thrown when the extent array does not have exactly 6 elements.
    /// </exception>
    public (int XMin, int XMax, int YMin, int YMax, int ZMin, int ZMax) Extent
    {
        get
        {
            var (nx, ny, nz) = Dimensions;
            return (
                _offsetX, _offsetX + nx - 1,
                _offsetY, _offsetY + ny - 1,
                _offsetZ, _offsetZ + nz - 1);
        }
        set
        {
            int nx = value.XMax - value.XMin + 1;
            int ny = value.YMax - value.YMin + 1;
            int nz = value.ZMax - value.ZMin + 1;
            Dimensions = (nx, ny, nz);
            _offsetX = value.XMin;
            _offsetY = value.YMin;
            _offsetZ = value.ZMin;
        }
    }

    // ---------------------------------------------------------------
    //  Offset
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets or sets the index offset of the ImageData.
    /// <para>
    /// The offset defines the minimum extent for each axis and can be positive
    /// or negative. In physical space the offset is relative to the image's
    /// <see cref="Origin"/>.
    /// </para>
    /// </summary>
    public (int X, int Y, int Z) Offset
    {
        get => (_offsetX, _offsetY, _offsetZ);
        set
        {
            var (nx, ny, nz) = Dimensions;
            Extent = (
                value.X, value.X + nx - 1,
                value.Y, value.Y + ny - 1,
                value.Z, value.Z + nz - 1);
        }
    }

    // ---------------------------------------------------------------
    //  Direction matrix
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets or sets the 3×3 direction matrix stored as a 9-element row-major array.
    /// <para>
    /// The direction matrix controls the orientation of the image data in physical space.
    /// </para>
    /// </summary>
    /// <exception cref="ArgumentNullException">
    /// Thrown when the value is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when the array does not have exactly 9 elements.
    /// </exception>
    public double[] DirectionMatrix
    {
        get => (double[])_directionMatrix.Clone();
        set
        {
            ArgumentNullException.ThrowIfNull(value);
            if (value.Length != 9)
            {
                throw new ArgumentException(
                    "Direction matrix must have exactly 9 elements (3×3 row-major).", nameof(value));
            }

            _directionMatrix = (double[])value.Clone();
        }
    }

    // ---------------------------------------------------------------
    //  Points (computed from origin, spacing, dimensions, offset, direction)
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets a copy of the implicitly defined points as a flat row-major (N×3) array.
    /// <para>
    /// Points of an <see cref="ImageData"/> cannot be set directly. They are
    /// implicitly defined by the <see cref="Origin"/>, <see cref="Spacing"/>,
    /// <see cref="Grid.Dimensions"/>, <see cref="Offset"/>, and <see cref="DirectionMatrix"/>.
    /// </para>
    /// </summary>
    public new double[] Points
    {
        get
        {
            var (nx, ny, nz) = Dimensions;
            if (nx == 0 || ny == 0 || nz == 0)
            {
                return [];
            }

            double dx = _spacingX;
            double dy = _spacingY;
            double dz = _spacingZ;

            double ox = _originX + _offsetX * dx;
            double oy = _originY + _offsetY * dy;
            double oz = _originZ + _offsetZ * dz;

            // Build 1-D coordinate arrays
            var xArr = new double[nx];
            for (int i = 0; i < nx; i++) xArr[i] = ox + i * dx;
            var yArr = new double[ny];
            for (int i = 0; i < ny; i++) yArr[i] = oy + i * dy;
            var zArr = new double[nz];
            for (int i = 0; i < nz; i++) zArr[i] = oz + i * dz;

            int nPoints = nx * ny * nz;
            var pts = new double[nPoints * 3];

            bool hasDirection = !IsIdentity(_directionMatrix);

            // Iterate in Fortran (column-major) order: x varies slowest, z varies fastest
            int idx = 0;
            for (int iz = 0; iz < nz; iz++)
            {
                for (int iy = 0; iy < ny; iy++)
                {
                    for (int ix = 0; ix < nx; ix++)
                    {
                        double px = xArr[ix];
                        double py = yArr[iy];
                        double pz = zArr[iz];

                        if (hasDirection)
                        {
                            // Apply rotation around origin
                            double rx = px - _originX;
                            double ry = py - _originY;
                            double rz = pz - _originZ;
                            pts[idx++] = _directionMatrix[0] * rx + _directionMatrix[1] * ry + _directionMatrix[2] * rz + _originX;
                            pts[idx++] = _directionMatrix[3] * rx + _directionMatrix[4] * ry + _directionMatrix[5] * rz + _originY;
                            pts[idx++] = _directionMatrix[6] * rx + _directionMatrix[7] * ry + _directionMatrix[8] * rz + _originZ;
                        }
                        else
                        {
                            pts[idx++] = px;
                            pts[idx++] = py;
                            pts[idx++] = pz;
                        }
                    }
                }
            }

            return pts;
        }
    }

    // ---------------------------------------------------------------
    //  X / Y / Z point arrays
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets the X coordinates of all points.
    /// </summary>
    public double[] X => ExtractComponent(0);

    /// <summary>
    /// Gets the Y coordinates of all points.
    /// </summary>
    public double[] Y => ExtractComponent(1);

    /// <summary>
    /// Gets the Z coordinates of all points.
    /// </summary>
    public double[] Z => ExtractComponent(2);

    // ---------------------------------------------------------------
    //  Bounds
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets the axis-aligned bounding box of this image data.
    /// </summary>
    public new BoundsTuple Bounds
    {
        get
        {
            var pts = Points;
            int nPoints = pts.Length / 3;
            if (nPoints == 0)
            {
                return new BoundsTuple(0, 0, 0, 0, 0, 0);
            }

            double xMin = double.MaxValue, xMax = double.MinValue;
            double yMin = double.MaxValue, yMax = double.MinValue;
            double zMin = double.MaxValue, zMax = double.MinValue;

            for (int i = 0; i < nPoints; i++)
            {
                int off = i * 3;
                double x = pts[off], y = pts[off + 1], z = pts[off + 2];
                if (x < xMin) xMin = x;
                if (x > xMax) xMax = x;
                if (y < yMin) yMin = y;
                if (y > yMax) yMax = y;
                if (z < zMin) zMin = z;
                if (z > zMax) zMax = z;
            }

            return new BoundsTuple(xMin, xMax, yMin, yMax, zMin, zMax);
        }
    }

    // ---------------------------------------------------------------
    //  Index-to-physical and physical-to-index matrices
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets the 4×4 index-to-physical transformation matrix as a 16-element row-major array.
    /// <para>
    /// Transforms index space (IJK) to physical space (XYZ).
    /// </para>
    /// </summary>
    public double[] IndexToPhysicalMatrix
    {
        get
        {
            double[] m = new double[16];
            // Column 0: direction[col0] * spacingX
            m[0] = _directionMatrix[0] * _spacingX;
            m[4] = _directionMatrix[3] * _spacingX;
            m[8] = _directionMatrix[6] * _spacingX;
            // Column 1: direction[col1] * spacingY
            m[1] = _directionMatrix[1] * _spacingY;
            m[5] = _directionMatrix[4] * _spacingY;
            m[9] = _directionMatrix[7] * _spacingY;
            // Column 2: direction[col2] * spacingZ
            m[2] = _directionMatrix[2] * _spacingZ;
            m[6] = _directionMatrix[5] * _spacingZ;
            m[10] = _directionMatrix[8] * _spacingZ;
            // Column 3: origin
            m[3] = _originX;
            m[7] = _originY;
            m[11] = _originZ;
            // Last row
            m[12] = 0;
            m[13] = 0;
            m[14] = 0;
            m[15] = 1;
            return m;
        }
    }

    /// <summary>
    /// Gets the 4×4 physical-to-index transformation matrix as a 16-element row-major array.
    /// <para>
    /// Transforms physical space (XYZ) to index space (IJK). This is the inverse
    /// of <see cref="IndexToPhysicalMatrix"/>.
    /// </para>
    /// </summary>
    public double[] PhysicalToIndexMatrix
    {
        get
        {
            var fwd = IndexToPhysicalMatrix;
            return Invert4x4(fwd);
        }
    }

    // ---------------------------------------------------------------
    //  Cast to StructuredGrid
    // ---------------------------------------------------------------

    /// <summary>
    /// Casts this image data to a <see cref="StructuredGrid"/>.
    /// <para>
    /// The resulting structured grid has explicit point coordinates computed from
    /// the implicit geometry. Point data, cell data, and field data are copied.
    /// </para>
    /// </summary>
    /// <returns>A new <see cref="StructuredGrid"/>.</returns>
    public StructuredGrid CastToStructuredGrid()
    {
        var (nx, ny, nz) = Dimensions;
        var pts = Points;
        var sg = new StructuredGrid(pts, nx, ny, nz, deep: false);

        CopyDataArraysTo(sg);
        return sg;
    }

    // ---------------------------------------------------------------
    //  Cast to RectilinearGrid
    // ---------------------------------------------------------------

    /// <summary>
    /// Casts this image data to a <see cref="RectilinearGrid"/>.
    /// <para>
    /// Off-axis rotations (direction matrices that are not axis-aligned sign
    /// permutations) are not supported. An <see cref="InvalidOperationException"/>
    /// is thrown if the direction matrix contains off-axis rotation.
    /// Point data, cell data, and field data are copied.
    /// </para>
    /// </summary>
    /// <returns>A new <see cref="RectilinearGrid"/>.</returns>
    /// <exception cref="InvalidOperationException">
    /// Thrown when the direction matrix contains off-axis rotation.
    /// </exception>
    public RectilinearGrid CastToRectilinearGrid()
    {
        var (nx, ny, nz) = Dimensions;

        // Determine axis signs from the direction matrix
        double[] sign = new double[3];
        for (int row = 0; row < 3; row++)
        {
            for (int col = 0; col < 3; col++)
            {
                double val = _directionMatrix[row * 3 + col];
                if (row == col)
                {
                    // Diagonal: must be ±1 (or close)
                    if (Math.Abs(Math.Abs(val) - 1.0) > 1e-10 && Math.Abs(val) > 1e-10)
                    {
                        throw new InvalidOperationException(
                            "Rectilinear grid does not support off-axis rotations. " +
                            "Consider removing off-axis rotations from the DirectionMatrix, " +
                            "or casting to StructuredGrid instead.");
                    }

                    sign[row] = Math.Abs(val) < 1e-10 ? 1.0 : Math.Sign(val);
                }
                else if (Math.Abs(val) > 1e-10)
                {
                    throw new InvalidOperationException(
                        "Rectilinear grid does not support off-axis rotations. " +
                        "Consider removing off-axis rotations from the DirectionMatrix, " +
                        "or casting to StructuredGrid instead.");
                }
            }
        }

        // Build coordinate arrays using linspace-style generation
        var xCoords = GenerateCoords(nx, _offsetX, _spacingX, sign[0], _originX);
        var yCoords = GenerateCoords(ny, _offsetY, _spacingY, sign[1], _originY);
        var zCoords = GenerateCoords(nz, _offsetZ, _spacingZ, sign[2], _originZ);

        var rg = new RectilinearGrid(xCoords, yCoords, zCoords);

        CopyDataArraysTo(rg);
        return rg;
    }

    // ---------------------------------------------------------------
    //  GetAttributes
    // ---------------------------------------------------------------

    /// <inheritdoc />
    protected override List<(string Name, string Value)> GetAttributes()
    {
        var attrs = base.GetAttributes();
        attrs.Add(("Spacing", $"{_spacingX:E1}, {_spacingY:E1}, {_spacingZ:E1}"));
        return attrs;
    }

    // ---------------------------------------------------------------
    //  Copy helpers
    // ---------------------------------------------------------------

    /// <inheritdoc />
    public override void DeepCopy(DataObject source)
    {
        base.DeepCopy(source);
        if (source is ImageData img)
        {
            CopyFieldsFrom(img);
            _directionMatrix = (double[])img._directionMatrix.Clone();
        }
    }

    /// <inheritdoc />
    public override void ShallowCopy(DataObject source)
    {
        base.ShallowCopy(source);
        if (source is ImageData img)
        {
            CopyFieldsFrom(img);
            _directionMatrix = img._directionMatrix;
        }
    }

    // ---------------------------------------------------------------
    //  Private helpers
    // ---------------------------------------------------------------

    /// <summary>
    /// Copies scalar fields from another <see cref="ImageData"/> instance.
    /// </summary>
    private void CopyFieldsFrom(ImageData img)
    {
        _originX = img._originX;
        _originY = img._originY;
        _originZ = img._originZ;
        _spacingX = img._spacingX;
        _spacingY = img._spacingY;
        _spacingZ = img._spacingZ;
        _offsetX = img._offsetX;
        _offsetY = img._offsetY;
        _offsetZ = img._offsetZ;
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
    /// Extracts a single coordinate component (0=X, 1=Y, 2=Z) from the points array.
    /// </summary>
    private double[] ExtractComponent(int component)
    {
        var pts = Points;
        int n = pts.Length / 3;
        var result = new double[n];
        for (int i = 0; i < n; i++)
        {
            result[i] = pts[i * 3 + component];
        }

        return result;
    }

    /// <summary>
    /// Checks whether the given 9-element matrix is the 3×3 identity.
    /// </summary>
    private static bool IsIdentity(double[] m)
    {
        const double tol = 1e-15;
        return Math.Abs(m[0] - 1) < tol && Math.Abs(m[1]) < tol && Math.Abs(m[2]) < tol
            && Math.Abs(m[3]) < tol && Math.Abs(m[4] - 1) < tol && Math.Abs(m[5]) < tol
            && Math.Abs(m[6]) < tol && Math.Abs(m[7]) < tol && Math.Abs(m[8] - 1) < tol;
    }

    /// <summary>
    /// Generates a 1-D coordinate array for one axis using linspace-style spacing.
    /// </summary>
    private static double[] GenerateCoords(int dim, int offset, double spacing, double sign, double origin)
    {
        var coords = new double[dim];
        for (int i = 0; i < dim; i++)
        {
            coords[i] = (offset + i) * spacing * sign + origin;
        }

        return coords;
    }

    /// <summary>
    /// Inverts a 4×4 matrix stored as a 16-element row-major array.
    /// Returns the identity if the matrix is singular.
    /// </summary>
    private static double[] Invert4x4(double[] m)
    {
        // Compute cofactor matrix and determinant via Laplace expansion
        double a00 = m[0], a01 = m[1], a02 = m[2], a03 = m[3];
        double a10 = m[4], a11 = m[5], a12 = m[6], a13 = m[7];
        double a20 = m[8], a21 = m[9], a22 = m[10], a23 = m[11];
        double a30 = m[12], a31 = m[13], a32 = m[14], a33 = m[15];

        double b00 = a00 * a11 - a01 * a10;
        double b01 = a00 * a12 - a02 * a10;
        double b02 = a00 * a13 - a03 * a10;
        double b03 = a01 * a12 - a02 * a11;
        double b04 = a01 * a13 - a03 * a11;
        double b05 = a02 * a13 - a03 * a12;
        double b06 = a20 * a31 - a21 * a30;
        double b07 = a20 * a32 - a22 * a30;
        double b08 = a20 * a33 - a23 * a30;
        double b09 = a21 * a32 - a22 * a31;
        double b10 = a21 * a33 - a23 * a31;
        double b11 = a22 * a33 - a23 * a32;

        double det = b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;
        if (Math.Abs(det) < 1e-30)
        {
            return [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
        }

        double invDet = 1.0 / det;
        return
        [
            (a11 * b11 - a12 * b10 + a13 * b09) * invDet,
            (a02 * b10 - a01 * b11 - a03 * b09) * invDet,
            (a31 * b05 - a32 * b04 + a33 * b03) * invDet,
            (a22 * b04 - a21 * b05 - a23 * b03) * invDet,
            (a12 * b08 - a10 * b11 - a13 * b07) * invDet,
            (a00 * b11 - a02 * b08 + a03 * b07) * invDet,
            (a32 * b02 - a30 * b05 - a33 * b01) * invDet,
            (a20 * b05 - a22 * b02 + a23 * b01) * invDet,
            (a10 * b10 - a11 * b08 + a13 * b06) * invDet,
            (a01 * b08 - a00 * b10 - a03 * b06) * invDet,
            (a30 * b04 - a31 * b02 + a33 * b00) * invDet,
            (a21 * b02 - a20 * b04 - a23 * b00) * invDet,
            (a11 * b07 - a10 * b09 - a12 * b06) * invDet,
            (a00 * b09 - a01 * b07 + a02 * b06) * invDet,
            (a31 * b01 - a30 * b03 - a32 * b00) * invDet,
            (a20 * b03 - a21 * b01 + a22 * b00) * invDet,
        ];
    }
}
