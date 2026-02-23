namespace PyVista.Core.Utilities;

/// <summary>
/// Enumerates the supported cell types for quality evaluation.
/// <para>
/// These correspond to the VTK cell types commonly used with mesh quality filters.
/// </para>
/// </summary>
public enum CellType
{
    /// <summary>A triangle cell (VTK type 5).</summary>
    Triangle = 5,

    /// <summary>A quadrilateral cell (VTK type 9).</summary>
    Quad = 9,

    /// <summary>A tetrahedron cell (VTK type 10).</summary>
    Tetra = 10,

    /// <summary>A hexahedron cell (VTK type 12).</summary>
    Hexahedron = 12,

    /// <summary>A pyramid cell (VTK type 14).</summary>
    Pyramid = 14,

    /// <summary>A wedge cell (VTK type 13).</summary>
    Wedge = 13,
}

/// <summary>
/// Enumerates the available quality measures that can be computed for mesh cells.
/// <para>
/// Each measure captures a different geometric property of the cell. Not all
/// measures are applicable to every <see cref="CellType"/>.
/// </para>
/// </summary>
public enum QualityMeasure
{
    /// <summary>Cell area (2D cells).</summary>
    Area,

    /// <summary>Aspect ratio of the cell.</summary>
    AspectRatio,

    /// <summary>Frobenius-norm aspect ratio.</summary>
    AspectFrobenius,

    /// <summary>Gamma aspect metric (tetrahedra only).</summary>
    AspectGamma,

    /// <summary>Collapse ratio of the cell.</summary>
    CollapseRatio,

    /// <summary>Condition number of the Jacobian matrix.</summary>
    Condition,

    /// <summary>Diagonal ratio metric (hexahedra).</summary>
    Diagonal,

    /// <summary>Characteristic dimension of the cell.</summary>
    Dimension,

    /// <summary>Distortion of the cell shape.</summary>
    Distortion,

    /// <summary>Jacobian determinant at cell vertices.</summary>
    Jacobian,

    /// <summary>Maximum interior angle.</summary>
    MaxAngle,

    /// <summary>Maximum Frobenius-norm aspect ratio.</summary>
    MaxAspectFrobenius,

    /// <summary>Maximum edge length ratio.</summary>
    MaxEdgeRatio,

    /// <summary>Median Frobenius-norm aspect ratio.</summary>
    MedAspectFrobenius,

    /// <summary>Minimum interior angle.</summary>
    MinAngle,

    /// <summary>Oddy metric measuring combined deviations from ideal shape.</summary>
    Oddy,

    /// <summary>Ratio of inscribed to circumscribed sphere radii.</summary>
    RadiusRatio,

    /// <summary>Relative size squared (compared to average cell).</summary>
    RelativeSizeSquared,

    /// <summary>Scaled Jacobian determinant normalized by edge lengths.</summary>
    ScaledJacobian,

    /// <summary>Shape metric based on the Jacobian matrix.</summary>
    Shape,

    /// <summary>Combined shape and size metric.</summary>
    ShapeAndSize,

    /// <summary>Shear metric of the cell.</summary>
    Shear,

    /// <summary>Combined shear and size metric.</summary>
    ShearAndSize,

    /// <summary>Skew of the cell from ideal orientation.</summary>
    Skew,

    /// <summary>Stretch ratio measuring deviation from ideal edge lengths.</summary>
    Stretch,

    /// <summary>Taper metric for quad and hex cells.</summary>
    Taper,

    /// <summary>Cell volume (3D cells).</summary>
    Volume,

    /// <summary>Warpage of the cell surface.</summary>
    Warpage,
}

/// <summary>
/// Holds information about a specific quality measure for a given cell type.
/// <para>
/// This is the C# equivalent of the Python <c>CellQualityInfo</c> dataclass.
/// It describes the acceptable, normal, and full value ranges as well as the
/// reference unit-cell value for the measure.
/// </para>
/// </summary>
/// <remarks>
/// Range information is based on the
/// <see href="https://github.com/sandialabs/verdict/raw/master/SAND2007-2853p.pdf">
/// Verdict Library Reference Manual</see>.
/// </remarks>
public sealed class CellQualityInfo
{
    /// <summary>
    /// Initializes a new instance of the <see cref="CellQualityInfo"/> class.
    /// </summary>
    /// <param name="cellType">The cell type.</param>
    /// <param name="qualityMeasure">The quality measure.</param>
    /// <param name="acceptableRange">Range for well-behaved cells.</param>
    /// <param name="normalRange">Range for all non-degenerate cells.</param>
    /// <param name="fullRange">Range for all cells including degenerate ones.</param>
    /// <param name="unitCellValue">The measure value for a reference unit cell.</param>
    public CellQualityInfo(
        CellType cellType,
        QualityMeasure qualityMeasure,
        (double Min, double Max) acceptableRange,
        (double Min, double Max) normalRange,
        (double Min, double Max) fullRange,
        double unitCellValue)
    {
        CellType = cellType;
        QualityMeasure = qualityMeasure;
        AcceptableRange = acceptableRange;
        NormalRange = normalRange;
        FullRange = fullRange;
        UnitCellValue = unitCellValue;
    }

    /// <summary>Gets the cell type this info applies to.</summary>
    public CellType CellType { get; }

    /// <summary>Gets the quality measure this info describes.</summary>
    public QualityMeasure QualityMeasure { get; }

    /// <summary>
    /// Gets the acceptable range. Well-behaved cells have values in this range.
    /// </summary>
    public (double Min, double Max) AcceptableRange { get; }

    /// <summary>
    /// Gets the normal range. All cells except those with degeneracies have values
    /// in this range.
    /// </summary>
    public (double Min, double Max) NormalRange { get; }

    /// <summary>
    /// Gets the full range. All cells, including degenerate ones, have values
    /// in this range.
    /// </summary>
    public (double Min, double Max) FullRange { get; }

    /// <summary>
    /// Gets the quality measure value for a reference unit cell (e.g., an equilateral
    /// triangle with edge length one for triangle cells).
    /// </summary>
    public double UnitCellValue { get; }

    /// <inheritdoc />
    public override string ToString() =>
        $"CellQualityInfo(CellType={CellType}, Measure={QualityMeasure}, " +
        $"Acceptable=({AcceptableRange.Min}, {AcceptableRange.Max}), " +
        $"UnitCell={UnitCellValue})";
}

/// <summary>
/// Provides static methods and lookup data for computing and querying cell quality metrics.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.core.utilities.cell_quality</c> module.
/// It contains a pre-populated table of <see cref="CellQualityInfo"/> entries and methods
/// to retrieve quality metadata for supported cell type / measure combinations.
/// </para>
/// </summary>
public static class CellQuality
{
    private static readonly double Sqrt2Over2 = Math.Sqrt(2.0) / 2.0;
    private static readonly double Sqrt3Over3 = Math.Sqrt(3.0) / 3.0;
    private static readonly double TetraAngle = (180.0 / Math.PI) * Math.Acos(1.0 / 3.0);
    private static readonly double Inf = double.PositiveInfinity;
    private static readonly double NInf = double.NegativeInfinity;

    /// <summary>
    /// The complete table of cell quality information entries.
    /// </summary>
    private static readonly CellQualityInfo[] InfoTable = BuildInfoTable();

    /// <summary>
    /// Lookup dictionary keyed by (CellType, QualityMeasure).
    /// </summary>
    private static readonly Dictionary<(CellType, QualityMeasure), CellQualityInfo> Lookup = BuildLookup();

    /// <summary>
    /// Returns information about a cell quality measure for a specific cell type.
    /// <para>
    /// The returned <see cref="CellQualityInfo"/> describes the acceptable, normal, and
    /// full value ranges as well as the reference unit-cell value for the measure.
    /// </para>
    /// </summary>
    /// <param name="cellType">The cell type to query.</param>
    /// <param name="measure">The quality measure to query.</param>
    /// <returns>A <see cref="CellQualityInfo"/> instance with range and unit-cell data.</returns>
    /// <exception cref="ArgumentException">
    /// Thrown when information is not available for the given cell type / measure combination.
    /// </exception>
    /// <example>
    /// <code>
    /// var info = CellQuality.GetCellQualityInfo(CellType.Triangle, QualityMeasure.ScaledJacobian);
    /// Console.WriteLine(info.AcceptableRange); // (0.5, 1.1547...)
    /// Console.WriteLine(info.UnitCellValue);   // 1.0
    /// </code>
    /// </example>
    public static CellQualityInfo GetCellQualityInfo(CellType cellType, QualityMeasure measure)
    {
        if (Lookup.TryGetValue((cellType, measure), out var info))
        {
            return info;
        }

        var validMeasures = GetSupportedMeasures(cellType);
        if (validMeasures.Length == 0)
        {
            throw new ArgumentException(
                $"Cell quality info is not available for cell type '{cellType}'. " +
                $"Valid cell types are: {string.Join(", ", GetSupportedCellTypes())}.");
        }

        throw new ArgumentException(
            $"Cell quality info is not available for '{cellType}' measure '{measure}'. " +
            $"Valid measures are: {string.Join(", ", validMeasures)}.");
    }

    /// <summary>
    /// Returns all <see cref="CellType"/> values that have at least one quality info entry.
    /// </summary>
    /// <returns>An array of supported cell types.</returns>
    public static CellType[] GetSupportedCellTypes()
    {
        var types = new HashSet<CellType>();
        foreach (var entry in InfoTable)
        {
            types.Add(entry.CellType);
        }
        return types.ToArray();
    }

    /// <summary>
    /// Returns all <see cref="QualityMeasure"/> values supported for the specified cell type.
    /// </summary>
    /// <param name="cellType">The cell type to query.</param>
    /// <returns>An array of supported quality measures for the cell type.</returns>
    public static QualityMeasure[] GetSupportedMeasures(CellType cellType)
    {
        var measures = new List<QualityMeasure>();
        foreach (var entry in InfoTable)
        {
            if (entry.CellType == cellType)
            {
                measures.Add(entry.QualityMeasure);
            }
        }
        return measures.ToArray();
    }

    /// <summary>
    /// Returns all <see cref="CellQualityInfo"/> entries for the specified cell type.
    /// </summary>
    /// <param name="cellType">The cell type to query.</param>
    /// <returns>An array of all quality info entries for the cell type.</returns>
    public static CellQualityInfo[] GetAllInfoForCellType(CellType cellType)
    {
        var results = new List<CellQualityInfo>();
        foreach (var entry in InfoTable)
        {
            if (entry.CellType == cellType)
            {
                results.Add(entry);
            }
        }
        return results.ToArray();
    }

    /// <summary>
    /// Computes a simple quality value for a cell defined by its vertices.
    /// <para>
    /// This is a lightweight, VTK-independent implementation that supports
    /// <see cref="QualityMeasure.Area"/> for triangles and
    /// <see cref="QualityMeasure.Volume"/> for tetrahedra. For other measures,
    /// use a full VTK-based pipeline.
    /// </para>
    /// </summary>
    /// <param name="cellType">The type of cell.</param>
    /// <param name="measure">The quality measure to compute.</param>
    /// <param name="vertices">
    /// A flat array of vertex coordinates in (x, y, z) order.
    /// The expected length depends on the cell type (e.g., 9 for triangles, 12 for tetrahedra).
    /// </param>
    /// <returns>The computed quality value.</returns>
    /// <exception cref="ArgumentException">
    /// Thrown when the vertex array length is invalid or the measure is not supported
    /// for the given cell type.
    /// </exception>
    /// <exception cref="NotSupportedException">
    /// Thrown when the cell type / measure combination requires a VTK-based implementation.
    /// </exception>
    public static double ComputeCellQuality(CellType cellType, QualityMeasure measure, double[] vertices)
    {
        ArgumentNullException.ThrowIfNull(vertices);

        return (cellType, measure) switch
        {
            (CellType.Triangle, QualityMeasure.Area) => ComputeTriangleArea(vertices),
            (CellType.Tetra, QualityMeasure.Volume) => ComputeTetraVolume(vertices),
            _ => throw new NotSupportedException(
                $"Computing '{measure}' for '{cellType}' requires a VTK-based pipeline. " +
                "Only Area (Triangle) and Volume (Tetra) are supported without VTK."),
        };
    }

    /// <summary>
    /// Computes the area of a triangle defined by three 3D vertices.
    /// </summary>
    private static double ComputeTriangleArea(double[] v)
    {
        if (v.Length != 9)
        {
            throw new ArgumentException(
                $"Triangle requires 9 coordinate values (3 vertices × 3), got {v.Length}.");
        }

        // Vectors AB and AC
        double abx = v[3] - v[0], aby = v[4] - v[1], abz = v[5] - v[2];
        double acx = v[6] - v[0], acy = v[7] - v[1], acz = v[8] - v[2];

        // Cross product AB × AC
        double cx = aby * acz - abz * acy;
        double cy = abz * acx - abx * acz;
        double cz = abx * acy - aby * acx;

        return 0.5 * Math.Sqrt(cx * cx + cy * cy + cz * cz);
    }

    /// <summary>
    /// Computes the volume of a tetrahedron defined by four 3D vertices.
    /// </summary>
    private static double ComputeTetraVolume(double[] v)
    {
        if (v.Length != 12)
        {
            throw new ArgumentException(
                $"Tetrahedron requires 12 coordinate values (4 vertices × 3), got {v.Length}.");
        }

        // Vectors from vertex 0 to vertices 1, 2, 3
        double ax = v[3] - v[0], ay = v[4] - v[1], az = v[5] - v[2];
        double bx = v[6] - v[0], by = v[7] - v[1], bz = v[8] - v[2];
        double cx = v[9] - v[0], cy = v[10] - v[1], cz = v[11] - v[2];

        // Scalar triple product: a · (b × c)
        double tripleProduct =
            ax * (by * cz - bz * cy) -
            ay * (bx * cz - bz * cx) +
            az * (bx * cy - by * cx);

        return Math.Abs(tripleProduct) / 6.0;
    }

    /// <summary>
    /// Builds the pre-populated table of all cell quality info entries.
    /// </summary>
    private static CellQualityInfo[] BuildInfoTable()
    {
        var i = Inf;
        var ni = NInf;
        var s3o4 = Math.Sqrt(3.0) / 4.0;
        var s6o3 = Math.Sqrt(6.0) / 3.0;
        var s2o12 = Math.Sqrt(2.0) / 12.0;
        var s2o6 = Math.Sqrt(2.0) / 6.0;
        var r33x2 = 2.0 * Sqrt3Over3;

        return new CellQualityInfo[]
        {
            // Triangle measures
            new(CellType.Triangle, QualityMeasure.Area,             (0.0, i),    (0.0, i),    (0.0, i),    s3o4),
            new(CellType.Triangle, QualityMeasure.AspectRatio,      (1.0, 1.3),  (1.0, i),    (1.0, i),    1.0),
            new(CellType.Triangle, QualityMeasure.AspectFrobenius,  (1.0, 1.3),  (1.0, i),    (1.0, i),    1.0),
            new(CellType.Triangle, QualityMeasure.Condition,        (1.0, 1.3),  (1.0, i),    (1.0, i),    1.0),
            new(CellType.Triangle, QualityMeasure.Distortion,       (0.5, 1.0),  (0.0, 1.0),  (ni, i),     1.0),
            new(CellType.Triangle, QualityMeasure.MaxAngle,         (60.0, 90.0),(60.0,180.0),(0.0,180.0),  60.0),
            new(CellType.Triangle, QualityMeasure.MinAngle,         (30.0, 60.0),(0.0, 60.0), (0.0,360.0),  60.0),
            new(CellType.Triangle, QualityMeasure.ScaledJacobian,   (0.5, r33x2),(-r33x2,r33x2),(ni, i),    1.0),
            new(CellType.Triangle, QualityMeasure.RadiusRatio,      (1.0, 3.0),  (1.0, i),    (1.0, i),    1.0),
            new(CellType.Triangle, QualityMeasure.Shape,            (0.25, 1.0), (0.0, 1.0),  (0.0, 1.0),  1.0),
            new(CellType.Triangle, QualityMeasure.ShapeAndSize,     (0.25, 1.0), (0.0, 1.0),  (0.0, 1.0),  1.0),

            // Quad measures
            new(CellType.Quad, QualityMeasure.Area,                 (0.0, i),    (0.0, i),    (ni, i),     1.0),
            new(CellType.Quad, QualityMeasure.AspectRatio,          (1.0, 1.3),  (1.0, i),    (1.0, i),    1.0),
            new(CellType.Quad, QualityMeasure.Condition,            (1.0, 4.0),  (1.0, i),    (1.0, i),    1.0),
            new(CellType.Quad, QualityMeasure.Distortion,           (0.5, 1.0),  (0.0, 1.0),  (ni, i),     1.0),
            new(CellType.Quad, QualityMeasure.Jacobian,             (0.0, i),    (0.0, i),    (ni, i),     1.0),
            new(CellType.Quad, QualityMeasure.MaxAspectFrobenius,   (1.0, 1.3),  (1.0, i),    (1.0, i),    1.0),
            new(CellType.Quad, QualityMeasure.MaxAngle,             (90.0,135.0),(90.0,360.0),(0.0,360.0),  90.0),
            new(CellType.Quad, QualityMeasure.MaxEdgeRatio,         (1.0, 1.3),  (1.0, i),    (1.0, i),    1.0),
            new(CellType.Quad, QualityMeasure.MedAspectFrobenius,   (1.0, 1.3),  (1.0, i),    (1.0, i),    1.0),
            new(CellType.Quad, QualityMeasure.MinAngle,             (45.0, 90.0),(0.0, 90.0), (0.0,360.0),  90.0),
            new(CellType.Quad, QualityMeasure.Oddy,                 (0.0, 0.5),  (0.0, i),    (0.0, i),    0.0),
            new(CellType.Quad, QualityMeasure.RadiusRatio,          (1.0, 1.3),  (1.0, i),    (1.0, i),    1.0),
            new(CellType.Quad, QualityMeasure.RelativeSizeSquared,  (0.3, 1.0),  (0.0, 1.0),  (0.0, 1.0),  1.0),
            new(CellType.Quad, QualityMeasure.ScaledJacobian,       (0.3, 1.0),  (-1.0, 1.0), (-1.0, 1.0), 1.0),
            new(CellType.Quad, QualityMeasure.Shape,                (0.3, 1.0),  (0.0, 1.0),  (0.0, 1.0),  1.0),
            new(CellType.Quad, QualityMeasure.ShapeAndSize,         (0.2, 1.0),  (0.0, 1.0),  (0.0, 1.0),  1.0),
            new(CellType.Quad, QualityMeasure.Shear,                (0.3, 1.0),  (0.0, 1.0),  (0.0, 1.0),  1.0),
            new(CellType.Quad, QualityMeasure.ShearAndSize,         (0.2, 1.0),  (0.0, 1.0),  (0.0, 1.0),  1.0),
            new(CellType.Quad, QualityMeasure.Skew,                 (0.0, 0.5),  (0.0, 1.0),  (0.0, 1.0),  0.0),
            new(CellType.Quad, QualityMeasure.Stretch,              (0.25, 1.0), (0.0, 1.0),  (0.0, i),    1.0),
            new(CellType.Quad, QualityMeasure.Taper,                (0.0, 0.7),  (0.0, i),    (0.0, i),    0.0),
            new(CellType.Quad, QualityMeasure.Warpage,              (0.3, 1.0),  (-1.0, 1.0), (ni, i),     1.0),

            // Tetra measures
            new(CellType.Tetra, QualityMeasure.AspectFrobenius,     (1.0, 1.3),  (1.0, i),    (1.0, i),    1.0),
            new(CellType.Tetra, QualityMeasure.AspectGamma,         (1.0, 3.0),  (1.0, i),    (1.0, i),    1.0),
            new(CellType.Tetra, QualityMeasure.AspectRatio,         (1.0, 3.0),  (1.0, i),    (1.0, i),    1.0),
            new(CellType.Tetra, QualityMeasure.CollapseRatio,       (0.1, i),    (0.0, i),    (0.0, i),    s6o3),
            new(CellType.Tetra, QualityMeasure.Condition,           (1.0, 3.0),  (1.0, i),    (1.0, i),    1.0),
            new(CellType.Tetra, QualityMeasure.Distortion,          (0.5, 1.0),  (0.0, 1.0),  (ni, i),     1.0),
            new(CellType.Tetra, QualityMeasure.Jacobian,            (0.0, i),    (0.0, i),    (ni, i),     Sqrt2Over2),
            new(CellType.Tetra, QualityMeasure.MinAngle,            (40.0, TetraAngle), (0.0, TetraAngle), (0.0, 360.0), TetraAngle),
            new(CellType.Tetra, QualityMeasure.RadiusRatio,         (1.0, 3.0),  (1.0, i),    (1.0, i),    1.0),
            new(CellType.Tetra, QualityMeasure.RelativeSizeSquared, (0.3, 1.0),  (0.0, 1.0),  (0.0, 1.0),  1.0),
            new(CellType.Tetra, QualityMeasure.ScaledJacobian,      (0.5, 1.0),  (-1.0, 1.0), (ni, i),     1.0),
            new(CellType.Tetra, QualityMeasure.Shape,               (0.3, 1.0),  (0.0, 1.0),  (0.0, 1.0),  1.0),
            new(CellType.Tetra, QualityMeasure.ShapeAndSize,        (0.2, 1.0),  (0.0, 1.0),  (0.0, 1.0),  1.0),
            new(CellType.Tetra, QualityMeasure.Volume,              (0.0, i),    (ni, i),     (ni, i),     s2o12),

            // Hexahedron measures
            new(CellType.Hexahedron, QualityMeasure.Diagonal,           (0.65, 1.0),(0.0, 1.0),  (0.0, i),    1.0),
            new(CellType.Hexahedron, QualityMeasure.Dimension,          (0.0, i),   (0.0, i),    (0.0, i),    Sqrt3Over3),
            new(CellType.Hexahedron, QualityMeasure.Distortion,         (0.5, 1.0), (0.0, 1.0),  (ni, i),     1.0),
            new(CellType.Hexahedron, QualityMeasure.Jacobian,           (0.0, i),   (0.0, i),    (ni, i),     1.0),
            new(CellType.Hexahedron, QualityMeasure.MaxEdgeRatio,       (1.0, 1.3), (1.0, i),    (1.0, i),    1.0),
            new(CellType.Hexahedron, QualityMeasure.MaxAspectFrobenius, (1.0, 3.0), (1.0, i),    (1.0, i),    1.0),
            new(CellType.Hexahedron, QualityMeasure.MedAspectFrobenius, (1.0, 3.0), (1.0, i),    (1.0, i),    1.0),
            new(CellType.Hexahedron, QualityMeasure.Oddy,               (0.0, 0.5), (0.0, i),    (0.0, i),    0.0),
            new(CellType.Hexahedron, QualityMeasure.RelativeSizeSquared,(0.5, 1.0), (0.0, 1.0),  (0.0, 1.0),  1.0),
            new(CellType.Hexahedron, QualityMeasure.ScaledJacobian,     (0.5, 1.0), (-1.0, 1.0), (-1.0, i),   1.0),
            new(CellType.Hexahedron, QualityMeasure.Shape,              (0.3, 1.0), (0.0, 1.0),  (0.0, 1.0),  1.0),
            new(CellType.Hexahedron, QualityMeasure.ShapeAndSize,       (0.2, 1.0), (0.0, 1.0),  (0.0, 1.0),  1.0),
            new(CellType.Hexahedron, QualityMeasure.Shear,              (0.3, 1.0), (0.0, 1.0),  (0.0, 1.0),  1.0),
            new(CellType.Hexahedron, QualityMeasure.ShearAndSize,       (0.2, 1.0), (0.0, 1.0),  (0.0, 1.0),  1.0),
            new(CellType.Hexahedron, QualityMeasure.Skew,               (0.0, 0.5), (0.0, 1.0),  (0.0, i),    0.0),
            new(CellType.Hexahedron, QualityMeasure.Stretch,            (0.25, 1.0),(0.0, 1.0),  (0.0, i),    1.0),
            new(CellType.Hexahedron, QualityMeasure.Taper,              (0.0, 0.5), (0.0, i),    (0.0, i),    0.0),
            new(CellType.Hexahedron, QualityMeasure.Volume,             (0.0, i),   (0.0, i),    (ni, i),     1.0),

            // Pyramid measures
            new(CellType.Pyramid, QualityMeasure.Volume, (0.0, i), (ni, i), (ni, i), s2o6),

            // Wedge measures
            new(CellType.Wedge, QualityMeasure.Volume, (0.0, i), (ni, i), (ni, i), s3o4),
        };
    }

    /// <summary>
    /// Builds the lookup dictionary from the info table.
    /// </summary>
    private static Dictionary<(CellType, QualityMeasure), CellQualityInfo> BuildLookup()
    {
        var dict = new Dictionary<(CellType, QualityMeasure), CellQualityInfo>();
        foreach (var entry in InfoTable)
        {
            dict[(entry.CellType, entry.QualityMeasure)] = entry;
        }
        return dict;
    }
}
