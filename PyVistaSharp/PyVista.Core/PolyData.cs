using System.Globalization;

using PyVista.Core.Cells;

namespace PyVista.Core;

/// <summary>
/// Dataset consisting of surface geometry including vertices, lines, polygons, and triangle strips.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.PolyData</c> class.
/// The surface geometry is defined by its <see cref="DataSet.Points"/> and four separate
/// cell connectivity arrays:
/// </para>
/// <list type="bullet">
///   <item><description><see cref="Verts"/> for vertex and poly-vertex cells.</description></item>
///   <item><description><see cref="Lines"/> for line and poly-line cells.</description></item>
///   <item><description><see cref="Faces"/> for triangle, quad, and polygon cells.</description></item>
///   <item><description><see cref="Strips"/> for triangle strip cells.</description></item>
/// </list>
/// </summary>
public class PolyData : PointSet
{
    private CellArray _verts = new();
    private CellArray _lines = new();
    private CellArray _faces = new();
    private CellArray _strips = new();

    /// <summary>
    /// Initializes a new instance of the <see cref="PolyData"/> class (empty).
    /// </summary>
    public PolyData()
    {
    }

    /// <summary>
    /// Initializes a new instance of the <see cref="PolyData"/> class from points.
    /// <para>
    /// When no cell connectivity arrays are provided, each point is automatically
    /// associated with a single vertex cell to create a point cloud.
    /// </para>
    /// </summary>
    /// <param name="points">Flat row-major point coordinates (length must be divisible by 3).</param>
    /// <param name="deep">When <c>true</c>, a copy of the point array is stored.</param>
    public PolyData(double[] points, bool deep = false)
        : base(points, deep)
    {
        _verts = MakeVertexCells(NPoints);
    }

    /// <summary>
    /// Initializes a new instance of the <see cref="PolyData"/> class from points
    /// and optional cell connectivity arrays.
    /// </summary>
    /// <param name="points">Flat row-major point coordinates (length must be divisible by 3).</param>
    /// <param name="faces">
    /// Padded face connectivity array, or <c>null</c>. Format:
    /// <c>[n0, p0_0, …, p0_n, n1, p1_0, …, p1_n, …]</c>.
    /// </param>
    /// <param name="lines">Padded line connectivity array, or <c>null</c>.</param>
    /// <param name="strips">Padded triangle strip connectivity array, or <c>null</c>.</param>
    /// <param name="verts">Padded vertex connectivity array, or <c>null</c>.</param>
    /// <param name="deep">When <c>true</c>, copies the point array.</param>
    public PolyData(
        double[] points,
        int[]? faces = null,
        int[]? lines = null,
        int[]? strips = null,
        int[]? verts = null,
        bool deep = false)
        : base(points, deep)
    {
        bool anyCells = faces is not null || lines is not null || strips is not null || verts is not null;

        if (!anyCells)
        {
            _verts = MakeVertexCells(NPoints);
        }
        else
        {
            if (verts is not null) _verts = new CellArray(verts);
            if (lines is not null) _lines = new CellArray(lines);
            if (faces is not null) _faces = new CellArray(faces);
            if (strips is not null) _strips = new CellArray(strips);
        }
    }

    /// <summary>
    /// Initializes a new instance of the <see cref="PolyData"/> class from a 2-D point array.
    /// </summary>
    /// <param name="points">An (N, 3) array of point coordinates.</param>
    /// <param name="deep">When <c>true</c>, copies the point data.</param>
    public PolyData(double[,] points, bool deep = false)
        : base(points, deep)
    {
        _verts = MakeVertexCells(NPoints);
    }

    // ---------------------------------------------------------------
    //  Verts
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets or sets the vertex padded connectivity array.
    /// <para>
    /// Like all padded VTK connectivity arrays, the format is:
    /// <c>[n0, p0_0, …, p0_n, n1, p1_0, …, p1_n, …]</c>
    /// where <c>n0</c> is the number of points in vertex 0.
    /// </para>
    /// </summary>
    public int[] Verts
    {
        get => _verts.Cells;
        set
        {
            ArgumentNullException.ThrowIfNull(value);
            _verts = new CellArray(value);
        }
    }

    /// <summary>
    /// Gets the number of vertex cells.
    /// </summary>
    public int NVerts => _verts.NCells;

    // ---------------------------------------------------------------
    //  Lines
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets or sets the line padded connectivity array.
    /// <para>
    /// Format: <c>[n0, p0_0, …, p0_n, n1, p1_0, …, p1_n, …]</c>
    /// </para>
    /// </summary>
    public int[] Lines
    {
        get => _lines.Cells;
        set
        {
            ArgumentNullException.ThrowIfNull(value);
            _lines = new CellArray(value);
        }
    }

    /// <summary>
    /// Gets the number of line cells.
    /// </summary>
    public int NLines => _lines.NCells;

    // ---------------------------------------------------------------
    //  Faces
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets or sets the polygonal faces padded connectivity array.
    /// <para>
    /// Format: <c>[n0, p0_0, …, p0_n, n1, p1_0, …, p1_n, …]</c>.
    /// Faces can be triangles, quads, or general polygons.
    /// </para>
    /// </summary>
    public int[] Faces
    {
        get => _faces.Cells;
        set
        {
            ArgumentNullException.ThrowIfNull(value);
            _faces = new CellArray(value);
        }
    }

    /// <summary>
    /// Gets the number of polygonal face cells (triangles, quads, polygons).
    /// </summary>
    public int NFacesStrict => _faces.NCells;

    // ---------------------------------------------------------------
    //  Strips
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets or sets the triangle strip padded connectivity array.
    /// <para>
    /// Format: <c>[n0, p0_0, …, p0_n, n1, p1_0, …, p1_n, …]</c>.
    /// </para>
    /// </summary>
    public int[] Strips
    {
        get => _strips.Cells;
        set
        {
            ArgumentNullException.ThrowIfNull(value);
            _strips = new CellArray(value);
        }
    }

    /// <summary>
    /// Gets the number of triangle strip cells.
    /// </summary>
    public int NStrips => _strips.NCells;

    // ---------------------------------------------------------------
    //  Total cell count
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets the total number of cells (verts + lines + faces + strips).
    /// </summary>
    public int NCellsTotal => NVerts + NLines + NFacesStrict + NStrips;

    // ---------------------------------------------------------------
    //  Regular / Irregular faces
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets or sets the face array when all faces have the same size.
    /// <para>
    /// Returns a 2-D array of shape <c>(NFaces, faceSize)</c>.
    /// </para>
    /// </summary>
    /// <exception cref="InvalidOperationException">
    /// Thrown when faces are not all the same size.
    /// </exception>
    public int[,] RegularFaces
    {
        get => _faces.RegularCells;
        set
        {
            ArgumentNullException.ThrowIfNull(value);
            _faces = CellArray.FromRegularCells(value);
        }
    }

    /// <summary>
    /// Gets or sets the faces as a jagged array of variable-length face arrays.
    /// </summary>
    public int[][] IrregularFaces
    {
        get
        {
            int nFaces = _faces.NCells;
            var offsets = _faces.OffsetArray;
            var connectivity = _faces.ConnectivityArray;
            var result = new int[nFaces][];
            for (int i = 0; i < nFaces; i++)
            {
                int start = offsets[i];
                int end = offsets[i + 1];
                int len = end - start;
                result[i] = new int[len];
                Array.Copy(connectivity, start, result[i], 0, len);
            }

            return result;
        }
        set
        {
            ArgumentNullException.ThrowIfNull(value);
            _faces = CellArray.FromIrregularCells(value);
        }
    }

    /// <summary>
    /// Creates a <see cref="PolyData"/> from a point array and a regular (equal-sized) face array.
    /// </summary>
    /// <param name="points">Flat row-major point coordinates.</param>
    /// <param name="faces">A 2-D array of shape <c>(nFaces, faceSize)</c>.</param>
    /// <param name="deep">When <c>true</c>, deep copies the face connectivity.</param>
    /// <returns>A new <see cref="PolyData"/> instance.</returns>
    public static PolyData FromRegularFaces(double[] points, int[,] faces, bool deep = false)
    {
        var cellArray = CellArray.FromRegularCells(faces, deep);
        var pd = new PolyData();
        pd.Points = deep ? (double[])points.Clone() : points;
        pd._faces = cellArray;
        return pd;
    }

    /// <summary>
    /// Creates a <see cref="PolyData"/> from a point array and a jagged face array.
    /// </summary>
    /// <param name="points">Flat row-major point coordinates.</param>
    /// <param name="faces">A jagged array of face connectivity.</param>
    /// <returns>A new <see cref="PolyData"/> instance.</returns>
    public static PolyData FromIrregularFaces(double[] points, int[][] faces)
    {
        var cellArray = CellArray.FromIrregularCells(faces);
        var pd = new PolyData();
        pd.Points = points;
        pd._faces = cellArray;
        return pd;
    }

    // ---------------------------------------------------------------
    //  IsAllTriangles
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets a value indicating whether all the faces are triangles and the mesh
    /// contains no vertices, lines, or strips.
    /// </summary>
    public bool IsAllTriangles
    {
        get
        {
            if (NFacesStrict == 0 || NLines > 0 || NVerts > 0 || NStrips > 0)
            {
                return false;
            }

            var connectivity = _faces.ConnectivityArray;
            if (connectivity.Length % 3 != 0)
            {
                return false;
            }

            var offsets = _faces.OffsetArray;
            for (int i = 0; i < _faces.NCells; i++)
            {
                if (offsets[i + 1] - offsets[i] != 3)
                {
                    return false;
                }
            }

            return true;
        }
    }

    // ---------------------------------------------------------------
    //  Volume
    // ---------------------------------------------------------------

    /// <summary>
    /// Returns the approximate volume of the dataset.
    /// <para>
    /// Computes the signed volume using the divergence theorem over
    /// triangulated faces. Requires a closed triangular surface for accurate results.
    /// </para>
    /// </summary>
    public new double Volume
    {
        get
        {
            // Triangulate first if needed, then compute signed volume via divergence theorem.
            var faces = GetTriangulatedFaces();
            double vol = 0.0;
            var pts = Points;
            foreach (var (a, b, c) in faces)
            {
                int oa = a * 3, ob = b * 3, oc = c * 3;
                double v321 = pts[oc] * pts[ob + 1] * pts[oa + 2];
                double v231 = pts[ob] * pts[oc + 1] * pts[oa + 2];
                double v312 = pts[oc] * pts[oa + 1] * pts[ob + 2];
                double v132 = pts[oa] * pts[oc + 1] * pts[ob + 2];
                double v213 = pts[ob] * pts[oa + 1] * pts[oc + 2];
                double v123 = pts[oa] * pts[ob + 1] * pts[oc + 2];
                vol += (-v321 + v231 + v312 - v132 - v213 + v123) / 6.0;
            }

            return Math.Abs(vol);
        }
    }

    // ---------------------------------------------------------------
    //  Normals
    // ---------------------------------------------------------------

    /// <summary>
    /// Returns the point normals.
    /// <para>
    /// If active point normals exist they are returned; otherwise normals are
    /// estimated from the face geometry.
    /// </para>
    /// </summary>
    public double[] PointNormals
    {
        get
        {
            var activeNormals = PointData.ActiveNormals;
            if (activeNormals is not null)
            {
                return activeNormals;
            }

            return ComputePointNormals();
        }
    }

    /// <summary>
    /// Returns the cell (face) normals.
    /// <para>
    /// If active cell normals exist they are returned; otherwise normals are
    /// computed from the face geometry.
    /// </para>
    /// </summary>
    public double[] CellNormals
    {
        get
        {
            var activeNormals = CellData.ActiveNormals;
            if (activeNormals is not null)
            {
                return activeNormals;
            }

            return ComputeFaceNormals();
        }
    }

    /// <summary>
    /// Returns the face normals. Alias to <see cref="CellNormals"/>.
    /// </summary>
    public double[] FaceNormals => CellNormals;

    // ---------------------------------------------------------------
    //  IsManifold / NOpenEdges
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets the number of open (boundary) edges on this mesh.
    /// <para>
    /// An edge is "open" if it belongs to only one face. Without VTK filters this
    /// is computed by counting edges that appear exactly once across all faces.
    /// </para>
    /// </summary>
    public int NOpenEdges
    {
        get
        {
            var edgeCounts = new Dictionary<(int, int), int>();
            int nFaces = _faces.NCells;
            var offsets = _faces.OffsetArray;
            var conn = _faces.ConnectivityArray;

            for (int f = 0; f < nFaces; f++)
            {
                int start = offsets[f];
                int end = offsets[f + 1];
                int faceSize = end - start;
                for (int j = 0; j < faceSize; j++)
                {
                    int p0 = conn[start + j];
                    int p1 = conn[start + (j + 1) % faceSize];
                    var edge = p0 < p1 ? (p0, p1) : (p1, p0);
                    edgeCounts.TryGetValue(edge, out int count);
                    edgeCounts[edge] = count + 1;
                }
            }

            int openCount = 0;
            foreach (var count in edgeCounts.Values)
            {
                if (count == 1) openCount++;
            }

            return openCount;
        }
    }

    /// <summary>
    /// Gets a value indicating whether the mesh is manifold (no open edges).
    /// </summary>
    public bool IsManifold => NOpenEdges == 0;

    // ---------------------------------------------------------------
    //  Copy helpers
    // ---------------------------------------------------------------

    /// <inheritdoc />
    public override void DeepCopy(DataObject source)
    {
        base.DeepCopy(source);
        if (source is PolyData pd)
        {
            _verts = new CellArray(pd._verts.Cells);
            _lines = new CellArray(pd._lines.Cells);
            _faces = new CellArray(pd._faces.Cells);
            _strips = new CellArray(pd._strips.Cells);
        }
    }

    /// <inheritdoc />
    public override void ShallowCopy(DataObject source)
    {
        base.ShallowCopy(source);
        if (source is PolyData pd)
        {
            _verts = pd._verts;
            _lines = pd._lines;
            _faces = pd._faces;
            _strips = pd._strips;
        }
    }

    /// <inheritdoc />
    protected override List<(string Name, string Value)> GetAttributes()
    {
        var attrs = base.GetAttributes();
        return attrs;
    }

    // ---------------------------------------------------------------
    //  Private helpers
    // ---------------------------------------------------------------

    /// <summary>
    /// Creates a CellArray of single-point vertex cells.
    /// </summary>
    private static CellArray MakeVertexCells(int nPoints)
    {
        if (nPoints == 0)
        {
            return new CellArray();
        }

        int[] cells = new int[nPoints * 2];
        for (int i = 0; i < nPoints; i++)
        {
            cells[i * 2] = 1;
            cells[i * 2 + 1] = i;
        }

        return new CellArray(cells);
    }

    /// <summary>
    /// Returns triangulated faces as a list of (a, b, c) point-index triples.
    /// Non-triangular faces are fan-triangulated from their first vertex.
    /// </summary>
    private List<(int A, int B, int C)> GetTriangulatedFaces()
    {
        var tris = new List<(int A, int B, int C)>();
        int nFaces = _faces.NCells;
        var offsets = _faces.OffsetArray;
        var conn = _faces.ConnectivityArray;

        for (int f = 0; f < nFaces; f++)
        {
            int start = offsets[f];
            int end = offsets[f + 1];
            int faceSize = end - start;
            if (faceSize < 3) continue;

            int p0 = conn[start];
            for (int j = 1; j < faceSize - 1; j++)
            {
                tris.Add((p0, conn[start + j], conn[start + j + 1]));
            }
        }

        return tris;
    }

    /// <summary>
    /// Computes per-face normals from face geometry.
    /// </summary>
    private double[] ComputeFaceNormals()
    {
        int nFaces = _faces.NCells;
        var normals = new double[nFaces * 3];
        var offsets = _faces.OffsetArray;
        var conn = _faces.ConnectivityArray;
        var pts = Points;

        for (int f = 0; f < nFaces; f++)
        {
            int start = offsets[f];
            int end = offsets[f + 1];
            int faceSize = end - start;

            if (faceSize < 3)
            {
                continue;
            }

            int p0 = conn[start], p1 = conn[start + 1], p2 = conn[start + 2];
            int o0 = p0 * 3, o1 = p1 * 3, o2 = p2 * 3;

            double ux = pts[o1] - pts[o0], uy = pts[o1 + 1] - pts[o0 + 1], uz = pts[o1 + 2] - pts[o0 + 2];
            double vx = pts[o2] - pts[o0], vy = pts[o2 + 1] - pts[o0 + 1], vz = pts[o2 + 2] - pts[o0 + 2];

            double nx = uy * vz - uz * vy;
            double ny = uz * vx - ux * vz;
            double nz = ux * vy - uy * vx;

            double len = Math.Sqrt(nx * nx + ny * ny + nz * nz);
            if (len > double.Epsilon)
            {
                nx /= len;
                ny /= len;
                nz /= len;
            }

            int nOff = f * 3;
            normals[nOff] = nx;
            normals[nOff + 1] = ny;
            normals[nOff + 2] = nz;
        }

        return normals;
    }

    /// <summary>
    /// Computes per-point normals by averaging adjacent face normals.
    /// </summary>
    private double[] ComputePointNormals()
    {
        var faceNorms = ComputeFaceNormals();
        int nPts = NPoints;
        var normals = new double[nPts * 3];
        var counts = new int[nPts];

        int nFaces = _faces.NCells;
        var offsets = _faces.OffsetArray;
        var conn = _faces.ConnectivityArray;

        for (int f = 0; f < nFaces; f++)
        {
            int start = offsets[f];
            int end = offsets[f + 1];
            int fOff = f * 3;

            for (int j = start; j < end; j++)
            {
                int pid = conn[j];
                int pOff = pid * 3;
                normals[pOff] += faceNorms[fOff];
                normals[pOff + 1] += faceNorms[fOff + 1];
                normals[pOff + 2] += faceNorms[fOff + 2];
                counts[pid]++;
            }
        }

        for (int i = 0; i < nPts; i++)
        {
            if (counts[i] > 0)
            {
                int pOff = i * 3;
                double len = Math.Sqrt(
                    normals[pOff] * normals[pOff] +
                    normals[pOff + 1] * normals[pOff + 1] +
                    normals[pOff + 2] * normals[pOff + 2]);
                if (len > double.Epsilon)
                {
                    normals[pOff] /= len;
                    normals[pOff + 1] /= len;
                    normals[pOff + 2] /= len;
                }
            }
        }

        return normals;
    }
}
