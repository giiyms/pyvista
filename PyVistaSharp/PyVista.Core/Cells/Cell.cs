using System.Globalization;
using System.Text;

namespace PyVista.Core.Cells;

/// <summary>
/// Represents a single cell defined by its topology (type and point IDs) and geometry (point
/// coordinates).
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.Cell</c> class. It provides the
/// capability to access a given cell's topology and can be useful when walking through a
/// cell's individual faces or investigating cell properties.
/// </para>
/// <para>
/// The cell object is a standalone copy and is not associated with any parent dataset.
/// Changing its data will not affect the original dataset.
/// </para>
/// </summary>
public sealed class Cell : IEquatable<Cell>
{
    /// <summary>
    /// Float format string used for bounds display, matching <c>pyvista.FLOAT_FORMAT</c>.
    /// </summary>
    private const string FloatFormat = "{0:E3}";

    private readonly CellType _cellType;
    private readonly double[,] _points;
    private readonly int[] _pointIds;

    /// <summary>
    /// Initializes a new instance of the <see cref="Cell"/> class.
    /// </summary>
    /// <param name="cellType">The VTK cell type.</param>
    /// <param name="points">
    /// Point coordinates as an N×3 array where N is the number of points in the cell.
    /// </param>
    /// <param name="pointIds">
    /// The point IDs that define this cell's connectivity in the parent dataset.
    /// </param>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="points"/> or <paramref name="pointIds"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when <paramref name="points"/> does not have exactly 3 columns or when the
    /// number of rows does not match the length of <paramref name="pointIds"/>.
    /// </exception>
    public Cell(CellType cellType, double[,] points, int[] pointIds)
    {
        ArgumentNullException.ThrowIfNull(points);
        ArgumentNullException.ThrowIfNull(pointIds);

        if (points.GetLength(1) != 3)
        {
            throw new ArgumentException("Points array must have exactly 3 columns (x, y, z).", nameof(points));
        }

        if (points.GetLength(0) != pointIds.Length)
        {
            throw new ArgumentException(
                $"Number of points ({points.GetLength(0)}) must match the number of point IDs ({pointIds.Length}).",
                nameof(points));
        }

        _cellType = cellType;
        _points = (double[,])points.Clone();
        _pointIds = (int[])pointIds.Clone();
    }

    /// <summary>
    /// Gets the cell type.
    /// </summary>
    /// <returns>The <see cref="CellType"/> of this cell.</returns>
    public CellType Type => _cellType;

    /// <summary>
    /// Gets a value indicating whether the cell uses linear interpolation.
    /// <para>
    /// Linear cells have type values less than or equal to 20 (the non-linear/quadratic
    /// cell types start at 21).
    /// </para>
    /// </summary>
    public bool IsLinear => (int)_cellType <= 20;

    /// <summary>
    /// Gets the topological dimension of the cell.
    /// <para>
    /// Returns 0 for vertices, 1 for edges/lines, 2 for surfaces, and 3 for volumes.
    /// </para>
    /// </summary>
    public int Dimension => GetCellDimension(_cellType);

    /// <summary>
    /// Gets the number of points composing this cell.
    /// </summary>
    public int NPoints => _pointIds.Length;

    /// <summary>
    /// Gets the number of faces composing this cell.
    /// <para>
    /// Only 3D cells have faces. Returns 0 for 0D, 1D, and 2D cells.
    /// </para>
    /// </summary>
    public int NFaces => GetFaceDefinitions(_cellType)?.Length ?? 0;

    /// <summary>
    /// Gets the number of edges composing this cell.
    /// <para>
    /// Only 1D+ cells have edges. Returns 0 for 0D cells.
    /// </para>
    /// </summary>
    public int NEdges => GetEdgeDefinitions(_cellType)?.Length ?? 0;

    /// <summary>
    /// Gets the point IDs that define this cell's connectivity.
    /// </summary>
    /// <returns>A copy of the point IDs array.</returns>
    public int[] PointIds => (int[])_pointIds.Clone();

    /// <summary>
    /// Gets the point coordinates of the cell as an N×3 array.
    /// </summary>
    /// <returns>A copy of the points array.</returns>
    public double[,] Points => (double[,])_points.Clone();

    /// <summary>
    /// Gets the cell bounds as <c>(XMin, XMax, YMin, YMax, ZMin, ZMax)</c>.
    /// </summary>
    public BoundsTuple Bounds
    {
        get
        {
            int n = _points.GetLength(0);
            if (n == 0)
            {
                return new BoundsTuple(0, 0, 0, 0, 0, 0);
            }

            double xMin = double.MaxValue, xMax = double.MinValue;
            double yMin = double.MaxValue, yMax = double.MinValue;
            double zMin = double.MaxValue, zMax = double.MinValue;

            for (int i = 0; i < n; i++)
            {
                double x = _points[i, 0], y = _points[i, 1], z = _points[i, 2];
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

    /// <summary>
    /// Gets the geometric center of the cell, computed as the average of its point
    /// coordinates.
    /// </summary>
    /// <returns>A tuple of <c>(x, y, z)</c> coordinates.</returns>
    public (double X, double Y, double Z) Center
    {
        get
        {
            int n = _points.GetLength(0);
            if (n == 0)
            {
                return (0.0, 0.0, 0.0);
            }

            double sx = 0, sy = 0, sz = 0;
            for (int i = 0; i < n; i++)
            {
                sx += _points[i, 0];
                sy += _points[i, 1];
                sz += _points[i, 2];
            }

            return (sx / n, sy / n, sz / n);
        }
    }

    /// <summary>
    /// Casts this cell to a <see cref="CellArray"/> and point array suitable for
    /// constructing a PolyData-like surface mesh.
    /// <para>
    /// Only 0D, 1D, and 2D cells can be cast. 3D cells will throw
    /// <see cref="InvalidOperationException"/>.
    /// </para>
    /// </summary>
    /// <returns>
    /// A tuple containing:
    /// <list type="bullet">
    ///   <item><description><c>Points</c> – a copy of the cell point coordinates (N×3).</description></item>
    ///   <item><description><c>Cells</c> – the legacy-format cell connectivity array.</description></item>
    ///   <item><description><c>CellKind</c> – a string indicating the topology kind
    ///     (<c>"verts"</c>, <c>"lines"</c>, <c>"strips"</c>, or <c>"faces"</c>).</description></item>
    /// </list>
    /// </returns>
    /// <exception cref="InvalidOperationException">
    /// Thrown when the cell has dimension 3.
    /// </exception>
    public (double[,] Points, int[] Cells, string CellKind) CastToPolyData()
    {
        int nPts = _pointIds.Length;
        int[] cells = BuildLegacyCells(nPts);
        double[,] pts = (double[,])_points.Clone();

        return Dimension switch
        {
            0 => (pts, cells, "verts"),
            1 => (pts, cells, "lines"),
            2 when _cellType == CellType.TriangleStrip => (pts, cells, "strips"),
            2 => (pts, cells, "faces"),
            _ => throw new InvalidOperationException(
                $"3D cells cannot be cast to PolyData: got cell type {_cellType}."),
        };
    }

    /// <summary>
    /// Casts this cell to arrays suitable for constructing an UnstructuredGrid.
    /// </summary>
    /// <returns>
    /// A tuple containing:
    /// <list type="bullet">
    ///   <item><description><c>CellConnectivity</c> – the legacy-format cell connectivity array.</description></item>
    ///   <item><description><c>CellTypes</c> – an array containing the single cell type.</description></item>
    ///   <item><description><c>Points</c> – a copy of the cell point coordinates (N×3).</description></item>
    /// </list>
    /// </returns>
    public (int[] CellConnectivity, int[] CellTypes, double[,] Points) CastToUnstructuredGrid()
    {
        int nPts = _pointIds.Length;

        int[] cellIds;
        if (_cellType == CellType.Polyhedron)
        {
            // Construct from faces for polyhedron cells.
            var facesDefs = GetFaceDefinitions(_cellType);
            var idList = new List<int> { facesDefs?.Length ?? 0 };
            if (facesDefs != null)
            {
                foreach (int[] faceDef in facesDefs)
                {
                    idList.Add(faceDef.Length);
                    foreach (int localIdx in faceDef)
                    {
                        // Map local face point index back to cell-local index.
                        idList.Add(localIdx);
                    }
                }
            }

            idList.Insert(0, idList.Count);
            cellIds = idList.ToArray();
        }
        else
        {
            cellIds = BuildLegacyCells(nPts);
        }

        return (cellIds, new[] { (int)_cellType }, (double[,])_points.Clone());
    }

    /// <summary>
    /// Gets the edge at the specified index.
    /// </summary>
    /// <param name="index">The zero-based edge index.</param>
    /// <returns>A new <see cref="Cell"/> representing the edge.</returns>
    /// <exception cref="IndexOutOfRangeException">
    /// Thrown when <paramref name="index"/> is out of range.
    /// </exception>
    public Cell GetEdge(int index)
    {
        int[][] edges = GetEdgeDefinitions(_cellType)
            ?? throw new InvalidOperationException(
                $"Cell type {_cellType} has no edges.");

        if (index < 0 || index >= edges.Length)
        {
            throw new IndexOutOfRangeException(
                $"Invalid index {index} for a cell with {edges.Length} edges.");
        }

        return BuildSubCell(edges[index], CellType.Line);
    }

    /// <summary>
    /// Gets the face at the specified index.
    /// </summary>
    /// <param name="index">The zero-based face index.</param>
    /// <returns>A new <see cref="Cell"/> representing the face.</returns>
    /// <exception cref="IndexOutOfRangeException">
    /// Thrown when <paramref name="index"/> is out of range.
    /// </exception>
    public Cell GetFace(int index)
    {
        int[][] faces = GetFaceDefinitions(_cellType)
            ?? throw new InvalidOperationException(
                $"Cell type {_cellType} has no faces.");

        if (index < 0 || index >= faces.Length)
        {
            throw new IndexOutOfRangeException(
                $"Invalid index {index} for a cell with {faces.Length} faces.");
        }

        int[] localIds = faces[index];
        CellType faceType = localIds.Length == 3 ? CellType.Triangle : CellType.Quad;
        return BuildSubCell(localIds, faceType);
    }

    /// <summary>
    /// Gets all edges composing this cell.
    /// </summary>
    /// <returns>A list of <see cref="Cell"/> objects, one per edge.</returns>
    public List<Cell> Edges
    {
        get
        {
            var result = new List<Cell>(NEdges);
            for (int i = 0; i < NEdges; i++)
            {
                result.Add(GetEdge(i));
            }

            return result;
        }
    }

    /// <summary>
    /// Gets all faces composing this cell.
    /// </summary>
    /// <returns>A list of <see cref="Cell"/> objects, one per face.</returns>
    public List<Cell> Faces
    {
        get
        {
            var result = new List<Cell>(NFaces);
            for (int i = 0; i < NFaces; i++)
            {
                result.Add(GetFace(i));
            }

            return result;
        }
    }

    /// <summary>
    /// Returns a copy of this cell.
    /// </summary>
    /// <param name="deep">
    /// When <c>true</c> (the default), a full deep copy is made.
    /// When <c>false</c>, the internal arrays are shared (shallow copy).
    /// </param>
    /// <returns>A new <see cref="Cell"/> instance.</returns>
    public Cell Copy(bool deep = true)
    {
        if (deep)
        {
            return new Cell(_cellType, (double[,])_points.Clone(), (int[])_pointIds.Clone());
        }

        // Shallow copy: share the same backing arrays.
        return new Cell(_cellType, _points, _pointIds, shallow: true);
    }

    /// <summary>
    /// Returns a console-friendly string representation of this cell.
    /// </summary>
    /// <returns>A formatted string describing the cell.</returns>
    public override string ToString()
    {
        var sb = new StringBuilder();
        sb.AppendLine($"Cell (0x{GetHashCode():X})");

        var bds = Bounds;
        var attrs = new (string Name, string Value)[]
        {
            ("Type", _cellType.ToString()),
            ("Linear", IsLinear.ToString()),
            ("Dimension", Dimension.ToString()),
            ("N Points", NPoints.ToString()),
            ("N Faces", NFaces.ToString()),
            ("N Edges", NEdges.ToString()),
            ("X Bounds", string.Format(CultureInfo.InvariantCulture, "{0:E3}, {1:E3}", bds.XMin, bds.XMax)),
            ("Y Bounds", string.Format(CultureInfo.InvariantCulture, "{0:E3}, {1:E3}", bds.YMin, bds.YMax)),
            ("Z Bounds", string.Format(CultureInfo.InvariantCulture, "{0:E3}, {1:E3}", bds.ZMin, bds.ZMax)),
        };

        int maxLen = 0;
        foreach (var attr in attrs)
        {
            if (attr.Name.Length > maxLen)
            {
                maxLen = attr.Name.Length;
            }
        }

        maxLen += 4;

        foreach (var attr in attrs)
        {
            sb.AppendLine($"  {(attr.Name + ":").PadRight(maxLen)}{attr.Value}");
        }

        return sb.ToString().TrimEnd();
    }

    /// <inheritdoc />
    public bool Equals(Cell? other)
    {
        if (other is null) return false;
        if (ReferenceEquals(this, other)) return true;
        if (_cellType != other._cellType) return false;
        if (!_pointIds.AsSpan().SequenceEqual(other._pointIds)) return false;

        int n = _points.GetLength(0);
        if (n != other._points.GetLength(0)) return false;
        for (int i = 0; i < n; i++)
        {
            if (_points[i, 0] != other._points[i, 0] ||
                _points[i, 1] != other._points[i, 1] ||
                _points[i, 2] != other._points[i, 2])
            {
                return false;
            }
        }

        return true;
    }

    /// <inheritdoc />
    public override bool Equals(object? obj) => Equals(obj as Cell);

    /// <inheritdoc />
    public override int GetHashCode() => HashCode.Combine(_cellType, _pointIds.Length);

    /// <summary>
    /// Equality operator.
    /// </summary>
    public static bool operator ==(Cell? left, Cell? right) =>
        left is null ? right is null : left.Equals(right);

    /// <summary>
    /// Inequality operator.
    /// </summary>
    public static bool operator !=(Cell? left, Cell? right) => !(left == right);

    // ------------------------------------------------------------------
    //  Private helpers
    // ------------------------------------------------------------------

    /// <summary>
    /// Private constructor for shallow copies that skips cloning.
    /// </summary>
    private Cell(CellType cellType, double[,] points, int[] pointIds, bool shallow)
    {
        _cellType = cellType;
        _points = points;
        _pointIds = pointIds;
    }

    /// <summary>
    /// Builds a legacy-format cell connectivity array: <c>{ nPts, 0, 1, ..., nPts-1 }</c>.
    /// </summary>
    private static int[] BuildLegacyCells(int nPts)
    {
        int[] cells = new int[nPts + 1];
        cells[0] = nPts;
        for (int i = 0; i < nPts; i++)
        {
            cells[i + 1] = i;
        }

        return cells;
    }

    /// <summary>
    /// Builds a sub-cell (edge or face) from local point indices.
    /// </summary>
    private Cell BuildSubCell(int[] localIds, CellType subCellType)
    {
        var subPoints = new double[localIds.Length, 3];
        var subPointIds = new int[localIds.Length];

        for (int i = 0; i < localIds.Length; i++)
        {
            int li = localIds[i];
            subPoints[i, 0] = _points[li, 0];
            subPoints[i, 1] = _points[li, 1];
            subPoints[i, 2] = _points[li, 2];
            subPointIds[i] = _pointIds[li];
        }

        return new Cell(subCellType, subPoints, subPointIds);
    }

    /// <summary>
    /// Returns the topological dimension for a given cell type.
    /// </summary>
    private static int GetCellDimension(CellType cellType) => cellType switch
    {
        CellType.EmptyCell => 0,
        CellType.Vertex => 0,
        CellType.PolyVertex => 0,

        CellType.Line => 1,
        CellType.PolyLine => 1,
        CellType.QuadraticEdge => 1,
        CellType.CubicLine => 1,
        CellType.LagrangeCurve => 1,
        CellType.BezierCurve => 1,
        CellType.HigherOrderEdge => 1,
        CellType.ParametricCurve => 1,

        CellType.Triangle => 2,
        CellType.TriangleStrip => 2,
        CellType.Polygon => 2,
        CellType.Pixel => 2,
        CellType.Quad => 2,
        CellType.QuadraticTriangle => 2,
        CellType.QuadraticQuad => 2,
        CellType.BiquadraticQuad => 2,
        CellType.QuadraticLinearQuad => 2,
        CellType.BiquadraticTriangle => 2,
        CellType.QuadraticPolygon => 2,
        CellType.LagrangeTriangle => 2,
        CellType.LagrangeQuadrilateral => 2,
        CellType.BezierTriangle => 2,
        CellType.BezierQuadrilateral => 2,
        CellType.HigherOrderTriangle => 2,
        CellType.HigherOrderQuad => 2,
        CellType.HigherOrderPolygon => 2,
        CellType.ParametricSurface => 2,
        CellType.ParametricTriSurface => 2,
        CellType.ParametricQuadSurface => 2,

        _ => 3,
    };

    /// <summary>
    /// Returns the edge definitions for common linear cell types, or <c>null</c> if none.
    /// Each inner array contains the local point indices forming an edge.
    /// </summary>
    private static int[][]? GetEdgeDefinitions(CellType cellType) => cellType switch
    {
        // 0D cells — no edges
        CellType.EmptyCell or CellType.Vertex or CellType.PolyVertex => null,

        // 1D cells — the cell itself is the single edge
        CellType.Line => new[] { new[] { 0, 1 } },

        // 2D cells
        CellType.Triangle => new[] { new[] { 0, 1 }, new[] { 1, 2 }, new[] { 2, 0 } },
        CellType.Quad or CellType.Pixel => new[]
        {
            new[] { 0, 1 }, new[] { 1, 2 }, new[] { 2, 3 }, new[] { 3, 0 },
        },
        CellType.Polygon => null, // variable — not statically known

        // 3D cells
        CellType.Tetra => new[]
        {
            new[] { 0, 1 }, new[] { 1, 2 }, new[] { 2, 0 },
            new[] { 0, 3 }, new[] { 1, 3 }, new[] { 2, 3 },
        },
        CellType.Hexahedron or CellType.Voxel => new[]
        {
            new[] { 0, 1 }, new[] { 1, 2 }, new[] { 2, 3 }, new[] { 3, 0 },
            new[] { 4, 5 }, new[] { 5, 6 }, new[] { 6, 7 }, new[] { 7, 4 },
            new[] { 0, 4 }, new[] { 1, 5 }, new[] { 2, 6 }, new[] { 3, 7 },
        },
        CellType.Wedge => new[]
        {
            new[] { 0, 1 }, new[] { 1, 2 }, new[] { 2, 0 },
            new[] { 3, 4 }, new[] { 4, 5 }, new[] { 5, 3 },
            new[] { 0, 3 }, new[] { 1, 4 }, new[] { 2, 5 },
        },
        CellType.Pyramid => new[]
        {
            new[] { 0, 1 }, new[] { 1, 2 }, new[] { 2, 3 }, new[] { 3, 0 },
            new[] { 0, 4 }, new[] { 1, 4 }, new[] { 2, 4 }, new[] { 3, 4 },
        },
        CellType.PentagonalPrism => new[]
        {
            new[] { 0, 1 }, new[] { 1, 2 }, new[] { 2, 3 }, new[] { 3, 4 }, new[] { 4, 0 },
            new[] { 5, 6 }, new[] { 6, 7 }, new[] { 7, 8 }, new[] { 8, 9 }, new[] { 9, 5 },
            new[] { 0, 5 }, new[] { 1, 6 }, new[] { 2, 7 }, new[] { 3, 8 }, new[] { 4, 9 },
        },
        CellType.HexagonalPrism => new[]
        {
            new[] { 0, 1 }, new[] { 1, 2 }, new[] { 2, 3 }, new[] { 3, 4 }, new[] { 4, 5 }, new[] { 5, 0 },
            new[] { 6, 7 }, new[] { 7, 8 }, new[] { 8, 9 }, new[] { 9, 10 }, new[] { 10, 11 }, new[] { 11, 6 },
            new[] { 0, 6 }, new[] { 1, 7 }, new[] { 2, 8 }, new[] { 3, 9 }, new[] { 4, 10 }, new[] { 5, 11 },
        },
        _ => null,
    };

    /// <summary>
    /// Returns the face definitions for common linear 3D cell types, or <c>null</c> if
    /// the cell has no faces. Each inner array contains the local point indices forming a
    /// face.
    /// </summary>
    private static int[][]? GetFaceDefinitions(CellType cellType) => cellType switch
    {
        CellType.Tetra => new[]
        {
            new[] { 0, 1, 3 },
            new[] { 1, 2, 3 },
            new[] { 2, 0, 3 },
            new[] { 0, 2, 1 },
        },
        CellType.Hexahedron => new[]
        {
            new[] { 0, 4, 7, 3 },
            new[] { 1, 2, 6, 5 },
            new[] { 0, 1, 5, 4 },
            new[] { 3, 7, 6, 2 },
            new[] { 0, 3, 2, 1 },
            new[] { 4, 5, 6, 7 },
        },
        CellType.Voxel => new[]
        {
            new[] { 0, 4, 6, 2 },
            new[] { 1, 3, 7, 5 },
            new[] { 0, 1, 5, 4 },
            new[] { 2, 6, 7, 3 },
            new[] { 0, 2, 3, 1 },
            new[] { 4, 5, 7, 6 },
        },
        CellType.Wedge => new[]
        {
            new[] { 0, 2, 1 },
            new[] { 3, 4, 5 },
            new[] { 0, 1, 4, 3 },
            new[] { 1, 2, 5, 4 },
            new[] { 0, 3, 5, 2 },
        },
        CellType.Pyramid => new[]
        {
            new[] { 0, 3, 2, 1 },
            new[] { 0, 1, 4 },
            new[] { 1, 2, 4 },
            new[] { 2, 3, 4 },
            new[] { 3, 0, 4 },
        },
        CellType.PentagonalPrism => new[]
        {
            new[] { 0, 4, 3, 2, 1 },
            new[] { 5, 6, 7, 8, 9 },
            new[] { 0, 1, 6, 5 },
            new[] { 1, 2, 7, 6 },
            new[] { 2, 3, 8, 7 },
            new[] { 3, 4, 9, 8 },
            new[] { 4, 0, 5, 9 },
        },
        CellType.HexagonalPrism => new[]
        {
            new[] { 0, 5, 4, 3, 2, 1 },
            new[] { 6, 7, 8, 9, 10, 11 },
            new[] { 0, 1, 7, 6 },
            new[] { 1, 2, 8, 7 },
            new[] { 2, 3, 9, 8 },
            new[] { 3, 4, 10, 9 },
            new[] { 4, 5, 11, 10 },
            new[] { 5, 0, 6, 11 },
        },
        _ => null,
    };
}
