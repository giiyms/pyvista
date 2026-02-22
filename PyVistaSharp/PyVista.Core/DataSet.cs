using System.Globalization;
using System.Text;

using PyVista.Core.Cells;

namespace PyVista.Core;

/// <summary>
/// Internal record representing a single cell: its type and the indices into
/// the point array that define its connectivity.
/// </summary>
internal sealed record CellRecord(CellType Type, int[] PointIds);

/// <summary>
/// Abstract base class for spatially-referenced datasets.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.DataSet</c> class. It adds
/// point/cell geometry and topology to the <see cref="DataObject"/> base class.
/// </para>
/// <para>
/// Because VTK C# bindings are not available, points are stored as a flat
/// <see cref="double"/>[] array (N×3 row-major), and cells are stored as a
/// <see cref="List{T}"/> of <see cref="CellRecord"/> instances.
/// </para>
/// </summary>
public abstract class DataSet : DataObject
{
    // ---------------------------------------------------------------
    //  Internal storage
    // ---------------------------------------------------------------

    /// <summary>Flat row-major point coordinates (length = NPoints * 3).</summary>
    private double[] _points = Array.Empty<double>();

    /// <summary>Cell topology.</summary>
    private readonly List<CellRecord> _cells = new();

    /// <summary>Point-associated data arrays.</summary>
    private readonly Dictionary<string, double[]> _pointDataArrays = new();

    /// <summary>Cell-associated data arrays.</summary>
    private readonly Dictionary<string, double[]> _cellDataArrays = new();

    // Active array tracking
    private string? _lastActiveScalarsName;
    private ActiveArrayInfoTuple _activeScalarsInfo = new(FieldAssociation.Point, null);
    private ActiveArrayInfoTuple _activeVectorsInfo = new(FieldAssociation.Point, null);
    private ActiveArrayInfoTuple _activeTensorsInfo = new(FieldAssociation.Point, null);

    /// <summary>
    /// Initializes a new instance of the <see cref="DataSet"/> class.
    /// </summary>
    protected DataSet()
    {
    }

    // ---------------------------------------------------------------
    //  Points
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets or sets the points of this dataset as a flat <see cref="double"/>[] in
    /// row-major order (length = N×3). Use <see cref="GetPoint"/> for convenient
    /// per-point access.
    /// </summary>
    /// <exception cref="ArgumentNullException">
    /// Thrown when the value is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when the array length is not divisible by 3.
    /// </exception>
    public double[] Points
    {
        get => _points;
        set
        {
            ArgumentNullException.ThrowIfNull(value);
            if (value.Length % 3 != 0)
            {
                throw new ArgumentException("Points array length must be divisible by 3.", nameof(value));
            }

            _points = value;
        }
    }

    /// <summary>
    /// Gets or sets the points of this dataset as a 2-D jagged-style
    /// <c>double[,]</c> with shape [N, 3].
    /// </summary>
    public double[,] Points2D
    {
        get
        {
            int n = NPoints;
            var result = new double[n, 3];
            for (int i = 0; i < n; i++)
            {
                int offset = i * 3;
                result[i, 0] = _points[offset];
                result[i, 1] = _points[offset + 1];
                result[i, 2] = _points[offset + 2];
            }

            return result;
        }
        set
        {
            ArgumentNullException.ThrowIfNull(value);
            if (value.GetLength(1) != 3)
            {
                throw new ArgumentException("Points array must have 3 columns.", nameof(value));
            }

            int n = value.GetLength(0);
            _points = new double[n * 3];
            for (int i = 0; i < n; i++)
            {
                int offset = i * 3;
                _points[offset] = value[i, 0];
                _points[offset + 1] = value[i, 1];
                _points[offset + 2] = value[i, 2];
            }
        }
    }

    /// <summary>
    /// Returns the coordinates of the <paramref name="index"/>-th point.
    /// </summary>
    /// <param name="index">Zero-based point index.</param>
    /// <returns>A tuple of (X, Y, Z) coordinates.</returns>
    public (double X, double Y, double Z) GetPoint(int index)
    {
        ValidatePointIndex(index);
        int offset = index * 3;
        return (_points[offset], _points[offset + 1], _points[offset + 2]);
    }

    /// <summary>
    /// Sets the coordinates of the <paramref name="index"/>-th point.
    /// </summary>
    /// <param name="index">Zero-based point index.</param>
    /// <param name="x">X coordinate.</param>
    /// <param name="y">Y coordinate.</param>
    /// <param name="z">Z coordinate.</param>
    public void SetPoint(int index, double x, double y, double z)
    {
        ValidatePointIndex(index);
        int offset = index * 3;
        _points[offset] = x;
        _points[offset + 1] = y;
        _points[offset + 2] = z;
    }

    // ---------------------------------------------------------------
    //  Point data / Cell data
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets the point data associated with this dataset.
    /// </summary>
    public DataSetAttributes PointData => new(_pointDataArrays, FieldAssociation.Point);

    /// <summary>
    /// Gets the cell data associated with this dataset.
    /// </summary>
    public DataSetAttributes CellData => new(_cellDataArrays, FieldAssociation.Cell);

    // ---------------------------------------------------------------
    //  Counts
    // ---------------------------------------------------------------

    /// <summary>Gets the number of points in this dataset.</summary>
    public int NPoints => _points.Length / 3;

    /// <summary>Gets the number of cells in this dataset.</summary>
    public int NCells => _cells.Count;

    /// <inheritdoc />
    public override bool IsEmpty => NPoints == 0;

    // ---------------------------------------------------------------
    //  Bounds / Center / Length
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets the axis-aligned bounding box of this dataset.
    /// </summary>
    public BoundsTuple Bounds
    {
        get
        {
            if (NPoints == 0)
            {
                return new BoundsTuple(0, 0, 0, 0, 0, 0);
            }

            double xMin = double.MaxValue, xMax = double.MinValue;
            double yMin = double.MaxValue, yMax = double.MinValue;
            double zMin = double.MaxValue, zMax = double.MinValue;

            int n = NPoints;
            for (int i = 0; i < n; i++)
            {
                int offset = i * 3;
                double x = _points[offset];
                double y = _points[offset + 1];
                double z = _points[offset + 2];
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
    /// Gets the center of the bounding box.
    /// </summary>
    public (double X, double Y, double Z) Center
    {
        get
        {
            var b = Bounds;
            return ((b.XMin + b.XMax) / 2.0, (b.YMin + b.YMax) / 2.0, (b.ZMin + b.ZMax) / 2.0);
        }
    }

    /// <summary>
    /// Gets the length of the diagonal of the bounding box.
    /// </summary>
    public double Length
    {
        get
        {
            var b = Bounds;
            double dx = b.XMax - b.XMin;
            double dy = b.YMax - b.YMin;
            double dz = b.ZMax - b.ZMin;
            return Math.Sqrt(dx * dx + dy * dy + dz * dz);
        }
    }

    // ---------------------------------------------------------------
    //  Active scalars info
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets the active scalar's association and name.
    /// </summary>
    public ActiveArrayInfoTuple ActiveScalarsInfo
    {
        get
        {
            var (field, name) = (_activeScalarsInfo.Association, _activeScalarsInfo.Name);

            var exclude = new HashSet<string> { "__custom_rgba", "Normals", "vtkOriginalPointIds", "TCoords" };
            if (name is not null && exclude.Contains(name))
            {
                name = _lastActiveScalarsName;
            }

            if (name is not null)
            {
                if (field == FieldAssociation.Cell)
                {
                    if (CellData.ActiveScalarsName != name) name = null;
                }
                else if (field == FieldAssociation.Point)
                {
                    if (PointData.ActiveScalarsName != name) name = null;
                }
            }

            if (name is null)
            {
                _activeScalarsInfo = new ActiveArrayInfoTuple(field, null);
                foreach (var attr in new[] { PointData, CellData })
                {
                    if (attr.ActiveScalarsName is not null)
                    {
                        _activeScalarsInfo = new ActiveArrayInfoTuple(attr.Association, attr.ActiveScalarsName);
                        break;
                    }
                }
            }

            return _activeScalarsInfo;
        }
    }

    /// <summary>
    /// Gets the active vector's association and name.
    /// </summary>
    public ActiveArrayInfoTuple ActiveVectorsInfo
    {
        get
        {
            var (field, name) = (_activeVectorsInfo.Association, _activeVectorsInfo.Name);

            if (name is not null)
            {
                if (field == FieldAssociation.Point && PointData.ActiveVectorsName != name)
                    name = null;
                if (field == FieldAssociation.Cell && CellData.ActiveVectorsName != name)
                    name = null;
            }

            if (name is null)
            {
                _activeVectorsInfo = new ActiveArrayInfoTuple(field, null);
                foreach (var attr in new[] { PointData, CellData })
                {
                    if (attr.ActiveVectorsName is not null)
                    {
                        _activeVectorsInfo = new ActiveArrayInfoTuple(attr.Association, attr.ActiveVectorsName);
                        break;
                    }
                }
            }

            return _activeVectorsInfo;
        }
    }

    /// <summary>
    /// Gets the active tensor's association and name.
    /// </summary>
    public ActiveArrayInfoTuple ActiveTensorsInfo => _activeTensorsInfo;

    // ---------------------------------------------------------------
    //  Active arrays (convenience)
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets the active scalars array, or <c>null</c> if none is active.
    /// </summary>
    public double[]? ActiveScalars
    {
        get
        {
            var info = ActiveScalarsInfo;
            if (info.Name is null) return null;
            try
            {
                return info.Association == FieldAssociation.Point
                    ? PointData[info.Name]
                    : info.Association == FieldAssociation.Cell
                        ? CellData[info.Name]
                        : null;
            }
            catch (KeyNotFoundException)
            {
                return null;
            }
        }
    }

    /// <summary>
    /// Gets the active vectors array, or <c>null</c> if none is active.
    /// </summary>
    public double[]? ActiveVectors
    {
        get
        {
            var info = ActiveVectorsInfo;
            if (info.Name is null) return null;
            try
            {
                return info.Association == FieldAssociation.Point
                    ? PointData[info.Name]
                    : info.Association == FieldAssociation.Cell
                        ? CellData[info.Name]
                        : null;
            }
            catch (KeyNotFoundException)
            {
                return null;
            }
        }
    }

    /// <summary>
    /// Gets the active tensors array, or <c>null</c> if none is active.
    /// </summary>
    public double[]? ActiveTensors
    {
        get
        {
            var info = ActiveTensorsInfo;
            if (info.Name is null) return null;
            try
            {
                return info.Association == FieldAssociation.Point
                    ? PointData[info.Name]
                    : info.Association == FieldAssociation.Cell
                        ? CellData[info.Name]
                        : null;
            }
            catch (KeyNotFoundException)
            {
                return null;
            }
        }
    }

    /// <summary>
    /// Gets the active normals array, or <c>null</c> if none is active.
    /// Returns point normals by default when both are present.
    /// </summary>
    public double[]? ActiveNormals
    {
        get
        {
            if (PointData.ActiveNormals is not null) return PointData.ActiveNormals;
            return CellData.ActiveNormals;
        }
    }

    // ---------------------------------------------------------------
    //  Active scalars / vectors / tensors name properties
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets or sets the name of the active scalars array.
    /// </summary>
    public string? ActiveScalarsName
    {
        get => ActiveScalarsInfo.Name;
        set => SetActiveScalars(value);
    }

    /// <summary>
    /// Gets or sets the name of the active vectors array.
    /// </summary>
    public string? ActiveVectorsName
    {
        get => ActiveVectorsInfo.Name;
        set => SetActiveVectors(value);
    }

    /// <summary>
    /// Gets or sets the name of the active tensors array.
    /// </summary>
    public string? ActiveTensorsName
    {
        get => ActiveTensorsInfo.Name;
        set => SetActiveTensors(value);
    }

    // ---------------------------------------------------------------
    //  Set active scalars / vectors / tensors
    // ---------------------------------------------------------------

    /// <summary>
    /// Finds the scalars by name and sets them as active.
    /// Pass <c>null</c> to deactivate.
    /// </summary>
    /// <param name="name">Name of the scalars array, or <c>null</c> to deactivate.</param>
    /// <param name="preference">Preferred association when both point and cell arrays match.</param>
    /// <returns>
    /// A tuple of the field association and the array (or <c>null</c> when deactivating).
    /// </returns>
    public (FieldAssociation Association, double[]? Array) SetActiveScalars(
        string? name,
        FieldAssociation preference = FieldAssociation.Cell)
    {
        if (name is null)
        {
            PointData.ActiveScalarsName = null;
            CellData.ActiveScalarsName = null;
            _activeScalarsInfo = new ActiveArrayInfoTuple(FieldAssociation.Point, null);
            return (FieldAssociation.None, null);
        }

        var field = FindArrayAssociation(name, preference);
        if (field == FieldAssociation.None)
        {
            if (FieldData.ContainsKey(name))
                throw new ArgumentException($"Data named \"{name}\" is a field array which cannot be active.");
            throw new KeyNotFoundException($"Data named \"{name}\" does not exist in this dataset.");
        }

        _lastActiveScalarsName = ActiveScalarsInfo.Name;

        if (field == FieldAssociation.Point)
            PointData.ActiveScalarsName = name;
        else if (field == FieldAssociation.Cell)
            CellData.ActiveScalarsName = name;

        _activeScalarsInfo = new ActiveArrayInfoTuple(field, name);

        return field == FieldAssociation.Point
            ? (field, PointData.ActiveScalars)
            : (field, CellData.ActiveScalars);
    }

    /// <summary>
    /// Finds the vectors by name and sets them as active.
    /// Pass <c>null</c> to deactivate.
    /// </summary>
    /// <param name="name">Name of the vectors array, or <c>null</c> to deactivate.</param>
    /// <param name="preference">Preferred association when both point and cell arrays match.</param>
    public void SetActiveVectors(string? name, FieldAssociation preference = FieldAssociation.Point)
    {
        if (name is null)
        {
            PointData.ActiveVectorsName = null;
            CellData.ActiveVectorsName = null;
            _activeVectorsInfo = new ActiveArrayInfoTuple(FieldAssociation.Point, null);
            return;
        }

        var field = FindArrayAssociation(name, preference);
        if (field == FieldAssociation.Point)
            PointData.ActiveVectorsName = name;
        else if (field == FieldAssociation.Cell)
            CellData.ActiveVectorsName = name;
        else
            throw new ArgumentException($"Data field ({name}) with type ({field}) not usable.");

        _activeVectorsInfo = new ActiveArrayInfoTuple(field, name);
    }

    /// <summary>
    /// Finds the tensors by name and sets them as active.
    /// Pass <c>null</c> to deactivate.
    /// </summary>
    /// <param name="name">Name of the tensors array, or <c>null</c> to deactivate.</param>
    /// <param name="preference">Preferred association when both point and cell arrays match.</param>
    public void SetActiveTensors(string? name, FieldAssociation preference = FieldAssociation.Point)
    {
        if (name is null)
        {
            _activeTensorsInfo = new ActiveArrayInfoTuple(FieldAssociation.Point, null);
            return;
        }

        var field = FindArrayAssociation(name, preference);
        if (field is not FieldAssociation.Point and not FieldAssociation.Cell)
            throw new ArgumentException($"Data field ({name}) with type ({field}) not usable.");

        _activeTensorsInfo = new ActiveArrayInfoTuple(field, name);
    }

    // ---------------------------------------------------------------
    //  Clear data helpers
    // ---------------------------------------------------------------

    /// <summary>Removes all point data arrays.</summary>
    public void ClearPointData() => PointData.Clear();

    /// <summary>Removes all cell data arrays.</summary>
    public void ClearCellData() => CellData.Clear();

    /// <summary>Removes all point, cell, and field data arrays.</summary>
    public void ClearData()
    {
        ClearPointData();
        ClearCellData();
        ClearFieldData();
    }

    // ---------------------------------------------------------------
    //  Copy from another dataset
    // ---------------------------------------------------------------

    /// <summary>
    /// Overwrites this dataset in-place with the geometry and data from
    /// <paramref name="source"/>.
    /// </summary>
    /// <param name="source">The source dataset.</param>
    /// <param name="deep">When <c>true</c>, performs a deep copy.</param>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="source"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when the source type is not compatible.
    /// </exception>
    public void CopyFrom(DataSet source, bool deep = true)
    {
        ArgumentNullException.ThrowIfNull(source);
        if (!GetType().IsAssignableFrom(source.GetType()))
        {
            throw new ArgumentException(
                $"The input DataSet type {source.GetType().Name} must be compatible with {GetType().Name}.");
        }

        if (deep)
            DeepCopy(source);
        else
            ShallowCopy(source);

        CopyMetaFrom(source, deep);
    }

    /// <inheritdoc />
    public override void DeepCopy(DataObject source)
    {
        base.DeepCopy(source);
        if (source is DataSet ds)
        {
            _points = (double[])ds._points.Clone();
            _cells.Clear();
            foreach (var cell in ds._cells)
                _cells.Add(new CellRecord(cell.Type, (int[])cell.PointIds.Clone()));

            CopyDataDictionary(ds._pointDataArrays, _pointDataArrays, deep: true);
            CopyDataDictionary(ds._cellDataArrays, _cellDataArrays, deep: true);
        }
    }

    /// <inheritdoc />
    public override void ShallowCopy(DataObject source)
    {
        base.ShallowCopy(source);
        if (source is DataSet ds)
        {
            _points = ds._points;
            _cells.Clear();
            _cells.AddRange(ds._cells);

            CopyDataDictionary(ds._pointDataArrays, _pointDataArrays, deep: false);
            CopyDataDictionary(ds._cellDataArrays, _cellDataArrays, deep: false);
        }
    }

    /// <inheritdoc />
    public override void CopyMetaFrom(DataObject source, bool deep = true)
    {
        base.CopyMetaFrom(source, deep);
        if (source is DataSet ds)
        {
            if (deep)
            {
                _activeScalarsInfo = new ActiveArrayInfoTuple(ds._activeScalarsInfo.Association, ds._activeScalarsInfo.Name);
                _activeVectorsInfo = new ActiveArrayInfoTuple(ds._activeVectorsInfo.Association, ds._activeVectorsInfo.Name);
                _activeTensorsInfo = new ActiveArrayInfoTuple(ds._activeTensorsInfo.Association, ds._activeTensorsInfo.Name);
            }
            else
            {
                _activeScalarsInfo = ds._activeScalarsInfo;
                _activeVectorsInfo = ds._activeVectorsInfo;
                _activeTensorsInfo = ds._activeTensorsInfo;
            }
        }
    }

    // ---------------------------------------------------------------
    //  Cast helpers
    // ---------------------------------------------------------------

    /// <summary>
    /// Returns a new <see cref="DataSet"/> instance of this dataset cast as a generic
    /// unstructured grid. Since VTK is not available, this creates a deep copy.
    /// </summary>
    /// <returns>A deep copy of this dataset.</returns>
    public DataSet CastToUnstructuredGrid()
    {
        var copy = (DataSet)Copy(deep: true);
        return copy;
    }

    /// <summary>
    /// Extracts the points of this dataset and returns a new <see cref="DataSet"/>
    /// containing only the points and point data.
    /// </summary>
    /// <param name="passCellData">
    /// When <c>true</c>, cell data is not preserved (no VTK cell_data_to_point_data filter).
    /// </param>
    /// <returns>A new dataset containing only points and point data.</returns>
    public DataSet CastToPointSet(bool passCellData = false)
    {
        var copy = (DataSet)Copy(deep: true);
        copy._cells.Clear();
        copy._cellDataArrays.Clear();
        if (!passCellData)
        {
            // No cell data conversion without VTK; clear cell data.
        }

        return copy;
    }

    // ---------------------------------------------------------------
    //  Cell access
    // ---------------------------------------------------------------

    /// <summary>
    /// Adds a cell to this dataset.
    /// </summary>
    /// <param name="type">The cell type.</param>
    /// <param name="pointIds">The point indices that define the cell connectivity.</param>
    public void AddCell(CellType type, params int[] pointIds)
    {
        ArgumentNullException.ThrowIfNull(pointIds);
        _cells.Add(new CellRecord(type, (int[])pointIds.Clone()));
    }

    /// <summary>
    /// Returns the cell at the specified index.
    /// </summary>
    /// <param name="index">Zero-based cell index.</param>
    /// <returns>A tuple of the cell type and a copy of its point IDs.</returns>
    /// <exception cref="IndexOutOfRangeException">
    /// Thrown when the index is out of range.
    /// </exception>
    public (CellType Type, int[] PointIds) GetCell(int index)
    {
        ValidateCellIndex(index);
        var cell = _cells[index];
        return (cell.Type, (int[])cell.PointIds.Clone());
    }

    /// <summary>
    /// Returns the neighbor cell IDs of the specified cell connected
    /// through shared points.
    /// </summary>
    /// <param name="index">Cell ID.</param>
    /// <param name="connections">
    /// How neighbors are connected: <c>"points"</c> (share any point).
    /// </param>
    /// <returns>List of neighboring cell IDs.</returns>
    public List<int> CellNeighbors(int index, string connections = "points")
    {
        ValidateCellIndex(index);
        if (connections != "points" && connections != "edges" && connections != "faces")
            throw new ArgumentException($"connections must be \"points\", \"edges\", or \"faces\", got \"{connections}\".");

        var cell = _cells[index];
        var cellPointSet = new HashSet<int>(cell.PointIds);

        var neighbors = new HashSet<int>();
        for (int i = 0; i < _cells.Count; i++)
        {
            if (i == index) continue;
            var other = _cells[i];
            int sharedCount = 0;
            foreach (int pid in other.PointIds)
            {
                if (cellPointSet.Contains(pid))
                    sharedCount++;
            }

            int required = connections switch
            {
                "edges" => 2,
                "faces" => 3,
                _ => 1, // "points"
            };

            if (sharedCount >= required)
                neighbors.Add(i);
        }

        return new List<int>(neighbors);
    }

    /// <summary>
    /// Returns the neighbor point IDs of the specified point (i.e. points that
    /// share at least one cell with this point).
    /// </summary>
    /// <param name="index">Point ID.</param>
    /// <returns>List of neighboring point IDs.</returns>
    public List<int> PointNeighbors(int index)
    {
        ValidatePointIndex(index);
        var neighbors = new HashSet<int>();
        for (int ci = 0; ci < _cells.Count; ci++)
        {
            var cell = _cells[ci];
            bool contains = false;
            foreach (int pid in cell.PointIds)
            {
                if (pid == index) { contains = true; break; }
            }

            if (!contains) continue;
            foreach (int pid in cell.PointIds)
            {
                if (pid != index) neighbors.Add(pid);
            }
        }

        return new List<int>(neighbors);
    }

    // ---------------------------------------------------------------
    //  Find closest / containing / along-line / within-bounds
    // ---------------------------------------------------------------

    /// <summary>
    /// Finds the index of the closest point to the given query point.
    /// </summary>
    /// <param name="point">A 3-element array (x, y, z).</param>
    /// <param name="n">Number of closest points to return (default 1).</param>
    /// <returns>
    /// The index of the closest point when <paramref name="n"/> is 1,
    /// or an array of the <paramref name="n"/> closest indices.
    /// </returns>
    public int[] FindClosestPoint(double[] point, int n = 1)
    {
        ArgumentNullException.ThrowIfNull(point);
        if (point.Length != 3)
            throw new ArgumentException("Point must have exactly 3 elements.", nameof(point));
        if (n < 1)
            throw new ArgumentOutOfRangeException(nameof(n), "n must be >= 1.");

        int npts = NPoints;
        if (npts == 0)
            return Array.Empty<int>();

        // Compute distances
        var distances = new (double Dist, int Index)[npts];
        for (int i = 0; i < npts; i++)
        {
            int offset = i * 3;
            double dx = _points[offset] - point[0];
            double dy = _points[offset + 1] - point[1];
            double dz = _points[offset + 2] - point[2];
            distances[i] = (dx * dx + dy * dy + dz * dz, i);
        }

        Array.Sort(distances, (a, b) => a.Dist.CompareTo(b.Dist));
        int count = Math.Min(n, npts);
        var result = new int[count];
        for (int i = 0; i < count; i++)
            result[i] = distances[i].Index;

        return result;
    }

    /// <summary>
    /// Finds the index of the cell whose center is closest to the given point.
    /// </summary>
    /// <param name="point">A 3-element array (x, y, z).</param>
    /// <returns>The index of the closest cell, or -1 if there are no cells.</returns>
    public int FindClosestCell(double[] point)
    {
        ArgumentNullException.ThrowIfNull(point);
        if (point.Length != 3)
            throw new ArgumentException("Point must have exactly 3 elements.", nameof(point));

        if (_cells.Count == 0) return -1;

        int bestIndex = -1;
        double bestDist = double.MaxValue;

        for (int ci = 0; ci < _cells.Count; ci++)
        {
            var (cx, cy, cz) = CellCenter(ci);
            double dx = cx - point[0];
            double dy = cy - point[1];
            double dz = cz - point[2];
            double dist = dx * dx + dy * dy + dz * dz;
            if (dist < bestDist)
            {
                bestDist = dist;
                bestIndex = ci;
            }
        }

        return bestIndex;
    }

    /// <summary>
    /// Finds the index of the cell that contains the given point based on
    /// cell bounding boxes. Returns -1 if no cell contains the point.
    /// </summary>
    /// <param name="point">A 3-element array (x, y, z).</param>
    /// <returns>Cell index, or -1 if not found.</returns>
    public int FindContainingCell(double[] point)
    {
        ArgumentNullException.ThrowIfNull(point);
        if (point.Length != 3)
            throw new ArgumentException("Point must have exactly 3 elements.", nameof(point));

        for (int ci = 0; ci < _cells.Count; ci++)
        {
            var b = CellBounds(ci);
            if (point[0] >= b.XMin && point[0] <= b.XMax &&
                point[1] >= b.YMin && point[1] <= b.YMax &&
                point[2] >= b.ZMin && point[2] <= b.ZMax)
            {
                return ci;
            }
        }

        return -1;
    }

    /// <summary>
    /// Finds cells whose bounding boxes intersect the line from
    /// <paramref name="pointA"/> to <paramref name="pointB"/>.
    /// </summary>
    /// <param name="pointA">Start of the line (length 3).</param>
    /// <param name="pointB">End of the line (length 3).</param>
    /// <param name="tolerance">Tolerance for intersection (default 0).</param>
    /// <returns>Array of cell indices.</returns>
    public int[] FindCellsAlongLine(double[] pointA, double[] pointB, double tolerance = 0.0)
    {
        ArgumentNullException.ThrowIfNull(pointA);
        ArgumentNullException.ThrowIfNull(pointB);
        if (pointA.Length != 3) throw new ArgumentException("pointA must have 3 elements.", nameof(pointA));
        if (pointB.Length != 3) throw new ArgumentException("pointB must have 3 elements.", nameof(pointB));

        var result = new List<int>();
        for (int ci = 0; ci < _cells.Count; ci++)
        {
            var b = CellBounds(ci);
            if (LineIntersectsBox(pointA, pointB, b, tolerance))
                result.Add(ci);
        }

        return result.ToArray();
    }

    /// <summary>
    /// Finds cells whose bounding boxes are within the given bounds.
    /// </summary>
    /// <param name="bounds">
    /// Bounding box as (xMin, xMax, yMin, yMax, zMin, zMax).
    /// </param>
    /// <returns>Array of cell indices within the bounds.</returns>
    public int[] FindCellsWithinBounds(BoundsTuple bounds)
    {
        var result = new List<int>();
        for (int ci = 0; ci < _cells.Count; ci++)
        {
            var cb = CellBounds(ci);
            if (cb.XMin >= bounds.XMin && cb.XMax <= bounds.XMax &&
                cb.YMin >= bounds.YMin && cb.YMax <= bounds.YMax &&
                cb.ZMin >= bounds.ZMin && cb.ZMax <= bounds.ZMax)
            {
                result.Add(ci);
            }
        }

        return result.ToArray();
    }

    // ---------------------------------------------------------------
    //  Per-cell queries
    // ---------------------------------------------------------------

    /// <summary>
    /// Returns the number of points in the specified cell.
    /// </summary>
    /// <param name="index">Zero-based cell index.</param>
    /// <returns>Number of points in the cell.</returns>
    public int CellNPoints(int index)
    {
        ValidateCellIndex(index);
        return _cells[index].PointIds.Length;
    }

    /// <summary>
    /// Returns the point indices of the specified cell.
    /// </summary>
    /// <param name="index">Zero-based cell index.</param>
    /// <returns>A copy of the point ID array for the cell.</returns>
    public int[] CellPoints(int index)
    {
        ValidateCellIndex(index);
        return (int[])_cells[index].PointIds.Clone();
    }

    /// <summary>
    /// Returns the bounding box of the specified cell.
    /// </summary>
    /// <param name="index">Zero-based cell index.</param>
    /// <returns>The bounding box of the cell.</returns>
    public BoundsTuple CellBounds(int index)
    {
        ValidateCellIndex(index);
        var ids = _cells[index].PointIds;
        if (ids.Length == 0)
            return new BoundsTuple(0, 0, 0, 0, 0, 0);

        double xMin = double.MaxValue, xMax = double.MinValue;
        double yMin = double.MaxValue, yMax = double.MinValue;
        double zMin = double.MaxValue, zMax = double.MinValue;

        foreach (int pid in ids)
        {
            int offset = pid * 3;
            double x = _points[offset];
            double y = _points[offset + 1];
            double z = _points[offset + 2];
            if (x < xMin) xMin = x;
            if (x > xMax) xMax = x;
            if (y < yMin) yMin = y;
            if (y > yMax) yMax = y;
            if (z < zMin) zMin = z;
            if (z > zMax) zMax = z;
        }

        return new BoundsTuple(xMin, xMax, yMin, yMax, zMin, zMax);
    }

    /// <summary>
    /// Returns the <see cref="Cells.CellType"/> of the specified cell.
    /// </summary>
    /// <param name="index">Zero-based cell index.</param>
    /// <returns>The cell type.</returns>
    public CellType CellType(int index)
    {
        ValidateCellIndex(index);
        return _cells[index].Type;
    }

    /// <summary>
    /// Returns the set of distinct cell types present in this dataset.
    /// </summary>
    public HashSet<CellType> DistinctCellTypes()
    {
        var types = new HashSet<CellType>();
        foreach (var cell in _cells)
            types.Add(cell.Type);
        return types;
    }

    /// <summary>
    /// Returns <c>true</c> if this dataset contains any non-linear (quadratic or
    /// higher-order) cells.
    /// </summary>
    public bool HasNonlinearCells()
    {
        foreach (var cell in _cells)
        {
            if ((int)cell.Type >= 21)
                return true;
        }

        return false;
    }

    // ---------------------------------------------------------------
    //  GetDataRange
    // ---------------------------------------------------------------

    /// <inheritdoc />
    public override (double Min, double Max) GetDataRange(
        string? name = null,
        FieldAssociation preference = FieldAssociation.Point)
    {
        if (name is null)
        {
            name = ActiveScalarsInfo.Name;
            if (name is null) return (double.NaN, double.NaN);
        }

        double[]? arr = null;
        if (preference == FieldAssociation.Point && PointData.ContainsKey(name))
            arr = PointData[name];
        else if (preference == FieldAssociation.Cell && CellData.ContainsKey(name))
            arr = CellData[name];
        else if (PointData.ContainsKey(name))
            arr = PointData[name];
        else if (CellData.ContainsKey(name))
            arr = CellData[name];
        else if (FieldData.ContainsKey(name))
            arr = FieldData[name];

        if (arr is null || arr.Length == 0) return (double.NaN, double.NaN);

        double min = double.MaxValue;
        double max = double.MinValue;
        foreach (double v in arr)
        {
            if (double.IsNaN(v)) continue;
            if (v < min) min = v;
            if (v > max) max = v;
        }

        return min <= max ? (min, max) : (double.NaN, double.NaN);
    }

    // ---------------------------------------------------------------
    //  Rotation methods
    // ---------------------------------------------------------------

    /// <summary>
    /// Rotates the dataset about the X axis by the given angle (in degrees).
    /// </summary>
    /// <param name="angle">Rotation angle in degrees.</param>
    /// <param name="inplace">When <c>true</c>, modifies this dataset; otherwise returns a copy.</param>
    /// <returns>The rotated dataset.</returns>
    public DataSet RotateX(double angle, bool inplace = false)
    {
        double rad = angle * Math.PI / 180.0;
        double cos = Math.Cos(rad);
        double sin = Math.Sin(rad);
        return ApplyRotation(
            (x, y, z) => (x, y * cos - z * sin, y * sin + z * cos),
            inplace);
    }

    /// <summary>
    /// Rotates the dataset about the Y axis by the given angle (in degrees).
    /// </summary>
    /// <param name="angle">Rotation angle in degrees.</param>
    /// <param name="inplace">When <c>true</c>, modifies this dataset; otherwise returns a copy.</param>
    /// <returns>The rotated dataset.</returns>
    public DataSet RotateY(double angle, bool inplace = false)
    {
        double rad = angle * Math.PI / 180.0;
        double cos = Math.Cos(rad);
        double sin = Math.Sin(rad);
        return ApplyRotation(
            (x, y, z) => (x * cos + z * sin, y, -x * sin + z * cos),
            inplace);
    }

    /// <summary>
    /// Rotates the dataset about the Z axis by the given angle (in degrees).
    /// </summary>
    /// <param name="angle">Rotation angle in degrees.</param>
    /// <param name="inplace">When <c>true</c>, modifies this dataset; otherwise returns a copy.</param>
    /// <returns>The rotated dataset.</returns>
    public DataSet RotateZ(double angle, bool inplace = false)
    {
        double rad = angle * Math.PI / 180.0;
        double cos = Math.Cos(rad);
        double sin = Math.Sin(rad);
        return ApplyRotation(
            (x, y, z) => (x * cos - y * sin, x * sin + y * cos, z),
            inplace);
    }

    // ---------------------------------------------------------------
    //  Flip methods
    // ---------------------------------------------------------------

    /// <summary>
    /// Reflects the dataset about the YZ plane through its center (flips X).
    /// </summary>
    /// <param name="inplace">When <c>true</c>, modifies this dataset; otherwise returns a copy.</param>
    /// <returns>The flipped dataset.</returns>
    public DataSet FlipX(bool inplace = false)
    {
        var c = Center;
        return ApplyRotation(
            (x, y, z) => (2.0 * c.X - x, y, z),
            inplace);
    }

    /// <summary>
    /// Reflects the dataset about the XZ plane through its center (flips Y).
    /// </summary>
    /// <param name="inplace">When <c>true</c>, modifies this dataset; otherwise returns a copy.</param>
    /// <returns>The flipped dataset.</returns>
    public DataSet FlipY(bool inplace = false)
    {
        var c = Center;
        return ApplyRotation(
            (x, y, z) => (x, 2.0 * c.Y - y, z),
            inplace);
    }

    /// <summary>
    /// Reflects the dataset about the XY plane through its center (flips Z).
    /// </summary>
    /// <param name="inplace">When <c>true</c>, modifies this dataset; otherwise returns a copy.</param>
    /// <returns>The flipped dataset.</returns>
    public DataSet FlipZ(bool inplace = false)
    {
        var c = Center;
        return ApplyRotation(
            (x, y, z) => (x, y, 2.0 * c.Z - z),
            inplace);
    }

    // ---------------------------------------------------------------
    //  ToString / Head
    // ---------------------------------------------------------------

    /// <inheritdoc />
    protected override List<(string Name, string Value)> GetAttributes()
    {
        var attrs = new List<(string Name, string Value)>
        {
            ("N Cells", NCells.ToString(CultureInfo.InvariantCulture)),
            ("N Points", NPoints.ToString(CultureInfo.InvariantCulture)),
        };

        var b = Bounds;
        string fmt = "F3";
        attrs.Add(("X Bounds", $"{b.XMin.ToString(fmt, CultureInfo.InvariantCulture)}, {b.XMax.ToString(fmt, CultureInfo.InvariantCulture)}"));
        attrs.Add(("Y Bounds", $"{b.YMin.ToString(fmt, CultureInfo.InvariantCulture)}, {b.YMax.ToString(fmt, CultureInfo.InvariantCulture)}"));
        attrs.Add(("Z Bounds", $"{b.ZMin.ToString(fmt, CultureInfo.InvariantCulture)}, {b.ZMax.ToString(fmt, CultureInfo.InvariantCulture)}"));

        return attrs;
    }

    // ---------------------------------------------------------------
    //  Private helpers
    // ---------------------------------------------------------------

    private void ValidatePointIndex(int index)
    {
        if (index < 0 || index >= NPoints)
            throw new IndexOutOfRangeException(
                $"Invalid index {index} for a dataset with {NPoints} points.");
    }

    private void ValidateCellIndex(int index)
    {
        if (index < 0 || index >= _cells.Count)
            throw new IndexOutOfRangeException(
                $"Invalid index {index} for a dataset with {_cells.Count} cells.");
    }

    private FieldAssociation FindArrayAssociation(string name, FieldAssociation preference)
    {
        bool inPoint = _pointDataArrays.ContainsKey(name);
        bool inCell = _cellDataArrays.ContainsKey(name);

        if (inPoint && inCell)
            return preference == FieldAssociation.Cell ? FieldAssociation.Cell : FieldAssociation.Point;
        if (inPoint) return FieldAssociation.Point;
        if (inCell) return FieldAssociation.Cell;
        return FieldAssociation.None;
    }

    private (double X, double Y, double Z) CellCenter(int index)
    {
        var ids = _cells[index].PointIds;
        if (ids.Length == 0) return (0, 0, 0);

        double sx = 0, sy = 0, sz = 0;
        foreach (int pid in ids)
        {
            int offset = pid * 3;
            sx += _points[offset];
            sy += _points[offset + 1];
            sz += _points[offset + 2];
        }

        double n = ids.Length;
        return (sx / n, sy / n, sz / n);
    }

    private DataSet ApplyRotation(
        Func<double, double, double, (double X, double Y, double Z)> transform,
        bool inplace)
    {
        DataSet target = inplace ? this : (DataSet)Copy(deep: true);
        int n = target.NPoints;
        for (int i = 0; i < n; i++)
        {
            int offset = i * 3;
            var (nx, ny, nz) = transform(
                target._points[offset],
                target._points[offset + 1],
                target._points[offset + 2]);
            target._points[offset] = nx;
            target._points[offset + 1] = ny;
            target._points[offset + 2] = nz;
        }

        return target;
    }

    private static bool LineIntersectsBox(double[] a, double[] b, BoundsTuple box, double tol)
    {
        double tMin = 0.0;
        double tMax = 1.0;

        double[] min = { box.XMin - tol, box.YMin - tol, box.ZMin - tol };
        double[] max = { box.XMax + tol, box.YMax + tol, box.ZMax + tol };

        for (int i = 0; i < 3; i++)
        {
            double d = b[i] - a[i];
            if (Math.Abs(d) < 1e-15)
            {
                if (a[i] < min[i] || a[i] > max[i])
                    return false;
            }
            else
            {
                double t1 = (min[i] - a[i]) / d;
                double t2 = (max[i] - a[i]) / d;
                if (t1 > t2) (t1, t2) = (t2, t1);
                tMin = Math.Max(tMin, t1);
                tMax = Math.Min(tMax, t2);
                if (tMin > tMax) return false;
            }
        }

        return true;
    }

    private static void CopyDataDictionary(
        Dictionary<string, double[]> source,
        Dictionary<string, double[]> dest,
        bool deep)
    {
        dest.Clear();
        foreach (var kvp in source)
        {
            dest[kvp.Key] = deep ? (double[])kvp.Value.Clone() : kvp.Value;
        }
    }
}
