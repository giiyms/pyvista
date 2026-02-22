using System.Collections;
using System.Globalization;
using System.Text;

namespace PyVista.Core;

/// <summary>
/// A composite dataset that stores a collection of <see cref="DataSet"/> partitions.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.PartitionedDataSet</c> class.
/// It wraps the concept of a VTK <c>vtkPartitionedDataSet</c> so that multiple
/// datasets can be stored together and accessed by index.
/// </para>
/// <para>
/// Unlike <see cref="MultiBlock"/>, partitions are not named and the structure
/// is flat (no nesting). Some operations such as deletion and pop are not
/// supported, matching the Python implementation.
/// </para>
/// </summary>
public class PartitionedDataSet : DataObject, IList<DataSet?>
{
    private readonly List<DataSet?> _partitions = new();

    // ---------------------------------------------------------------
    //  Constructors
    // ---------------------------------------------------------------

    /// <summary>
    /// Initializes a new, empty <see cref="PartitionedDataSet"/>.
    /// </summary>
    public PartitionedDataSet()
    {
    }

    /// <summary>
    /// Initializes a new <see cref="PartitionedDataSet"/> from a list of datasets.
    /// </summary>
    /// <param name="datasets">The datasets to include as partitions.</param>
    public PartitionedDataSet(IEnumerable<DataSet?> datasets)
    {
        ArgumentNullException.ThrowIfNull(datasets);
        foreach (var ds in datasets)
        {
            Append(ds);
        }
    }

    /// <summary>
    /// Initializes a new <see cref="PartitionedDataSet"/> as a copy of another instance.
    /// </summary>
    /// <param name="other">The source partitioned dataset to copy.</param>
    /// <param name="deep">When <c>true</c>, performs a deep copy of all partitions.</param>
    /// <exception cref="PartitionedDataSetsNotSupported">
    /// Thrown when <paramref name="deep"/> is <c>false</c> because shallow copy is not supported.
    /// </exception>
    public PartitionedDataSet(PartitionedDataSet other, bool deep = true)
    {
        ArgumentNullException.ThrowIfNull(other);
        if (deep)
        {
            DeepCopy(other);
        }
        else
        {
            throw new PartitionedDataSetsNotSupported(
                "Shallow copy is not supported for PartitionedDataSet.");
        }
    }

    // ---------------------------------------------------------------
    //  Partition count / emptiness
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets or sets the number of partitions.
    /// <para>
    /// Setting a value larger than the current count pads with <c>null</c> partitions.
    /// Setting a smaller value truncates.
    /// </para>
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException">
    /// Thrown when the value is negative.
    /// </exception>
    public int NPartitions
    {
        get => _partitions.Count;
        set
        {
            if (value < 0)
            {
                throw new ArgumentOutOfRangeException(nameof(value), "NPartitions cannot be negative.");
            }

            while (_partitions.Count < value)
            {
                _partitions.Add(null);
            }

            while (_partitions.Count > value)
            {
                _partitions.RemoveAt(_partitions.Count - 1);
            }
        }
    }

    /// <summary>
    /// Gets a value indicating whether there are no partitions.
    /// </summary>
    public override bool IsEmpty => _partitions.Count == 0;

    /// <summary>
    /// Gets the number of partitions. Equivalent to <see cref="NPartitions"/>.
    /// </summary>
    public int Count => _partitions.Count;

    /// <inheritdoc />
    bool ICollection<DataSet?>.IsReadOnly => false;

    // ---------------------------------------------------------------
    //  Bounds / Center / Volume / Length
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets the axis-aligned bounding box that encloses all partitions.
    /// </summary>
    public BoundsTuple Bounds
    {
        get
        {
            double xMin = double.MaxValue, xMax = double.MinValue;
            double yMin = double.MaxValue, yMax = double.MinValue;
            double zMin = double.MaxValue, zMax = double.MinValue;
            bool found = false;

            foreach (var partition in _partitions)
            {
                if (partition is DataSet ds && ds.NPoints > 0)
                {
                    var b = ds.Bounds;
                    if (b.XMin < xMin) xMin = b.XMin;
                    if (b.XMax > xMax) xMax = b.XMax;
                    if (b.YMin < yMin) yMin = b.YMin;
                    if (b.YMax > yMax) yMax = b.YMax;
                    if (b.ZMin < zMin) zMin = b.ZMin;
                    if (b.ZMax > zMax) zMax = b.ZMax;
                    found = true;
                }
            }

            if (!found)
            {
                return new BoundsTuple(0, 0, 0, 0, 0, 0);
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
            return ((b.XMin + b.XMax) / 2.0,
                    (b.YMin + b.YMax) / 2.0,
                    (b.ZMin + b.ZMax) / 2.0);
        }
    }

    /// <summary>
    /// Gets the total volume of all partitions.
    /// </summary>
    public double Volume
    {
        get
        {
            double total = 0;
            foreach (var partition in _partitions)
            {
                if (partition is PointSet ps)
                {
                    total += ps.Volume;
                }
            }

            return total;
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
    //  Indexer
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets or sets the partition at the specified integer index.
    /// <para>
    /// Negative indices are supported (e.g., <c>-1</c> for the last partition).
    /// </para>
    /// </summary>
    /// <param name="index">Zero-based index (negative indices count from the end).</param>
    /// <returns>The partition at the specified index, which may be <c>null</c>.</returns>
    /// <exception cref="IndexOutOfRangeException">Thrown when the index is out of range.</exception>
    public DataSet? this[int index]
    {
        get
        {
            int resolved = ResolveIndex(index);
            return _partitions[resolved];
        }
        set
        {
            int resolved = ResolveIndex(index);
            _partitions[resolved] = value;
        }
    }

    // ---------------------------------------------------------------
    //  Append / Extend / Insert
    // ---------------------------------------------------------------

    /// <summary>
    /// Appends a dataset as a new partition at the end.
    /// </summary>
    /// <param name="dataset">Dataset to append (may be <c>null</c>).</param>
    public void Append(DataSet? dataset)
    {
        _partitions.Add(dataset);
    }

    /// <summary>
    /// Extends this <see cref="PartitionedDataSet"/> with partitions from an iterable.
    /// </summary>
    /// <param name="datasets">The datasets to append.</param>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="datasets"/> is <c>null</c>.
    /// </exception>
    public void Extend(IEnumerable<DataSet?> datasets)
    {
        ArgumentNullException.ThrowIfNull(datasets);
        foreach (var ds in datasets)
        {
            Append(ds);
        }
    }

    /// <summary>
    /// Inserts a dataset before the specified index.
    /// </summary>
    /// <param name="index">Zero-based insertion index.</param>
    /// <param name="dataset">Dataset to insert (may be <c>null</c>).</param>
    public void Insert(int index, DataSet? dataset)
    {
        int resolved = index < 0 ? Math.Max(_partitions.Count + index, 0) : Math.Min(index, _partitions.Count);
        _partitions.Insert(resolved, dataset);
    }

    // ---------------------------------------------------------------
    //  Pop / Delete (not supported)
    // ---------------------------------------------------------------

    /// <summary>
    /// Removing a partition by index is not supported.
    /// </summary>
    /// <param name="index">Unused.</param>
    /// <exception cref="PartitionedDataSetsNotSupported">Always thrown.</exception>
    public DataSet? Pop(int index = -1)
    {
        throw new PartitionedDataSetsNotSupported(
            "Pop is not supported for PartitionedDataSet.");
    }

    /// <summary>
    /// Removing a partition at the specified index is not supported.
    /// </summary>
    /// <param name="index">Unused.</param>
    /// <exception cref="PartitionedDataSetsNotSupported">Always thrown.</exception>
    public void RemoveAt(int index)
    {
        throw new PartitionedDataSetsNotSupported(
            "Deletion is not supported for PartitionedDataSet.");
    }

    /// <summary>
    /// Removing a partition is not supported.
    /// </summary>
    /// <param name="item">Unused.</param>
    /// <returns>Does not return.</returns>
    /// <exception cref="PartitionedDataSetsNotSupported">Always thrown.</exception>
    public bool Remove(DataSet? item)
    {
        throw new PartitionedDataSetsNotSupported(
            "Deletion is not supported for PartitionedDataSet.");
    }

    // ---------------------------------------------------------------
    //  Clean
    // ---------------------------------------------------------------

    /// <summary>
    /// Removes all <c>null</c> partitions in place. Optionally also removes empty datasets.
    /// </summary>
    /// <param name="empty">
    /// When <c>true</c> (default), also removes datasets with zero points.
    /// </param>
    public void Clean(bool empty = true)
    {
        for (int i = _partitions.Count - 1; i >= 0; i--)
        {
            var data = _partitions[i];
            if (data is null)
            {
                _partitions.RemoveAt(i);
            }
            else if (empty && data.NPoints == 0)
            {
                _partitions.RemoveAt(i);
            }
        }
    }

    /// <summary>
    /// Removes all partitions.
    /// </summary>
    public void Clear()
    {
        _partitions.Clear();
    }

    // ---------------------------------------------------------------
    //  Copy
    // ---------------------------------------------------------------

    /// <summary>
    /// Returns a deep copy of this <see cref="PartitionedDataSet"/>.
    /// </summary>
    /// <param name="deep">
    /// When <c>true</c> (default), performs a deep copy of all partitions.
    /// Shallow copy is not supported and throws <see cref="PartitionedDataSetsNotSupported"/>.
    /// </param>
    /// <returns>A new <see cref="PartitionedDataSet"/>.</returns>
    /// <exception cref="PartitionedDataSetsNotSupported">
    /// Thrown when <paramref name="deep"/> is <c>false</c>.
    /// </exception>
    public new PartitionedDataSet Copy(bool deep = true)
    {
        return new PartitionedDataSet(this, deep);
    }

    /// <inheritdoc />
    public override void DeepCopy(DataObject source)
    {
        base.DeepCopy(source);
        if (source is PartitionedDataSet other)
        {
            _partitions.Clear();
            for (int i = 0; i < other._partitions.Count; i++)
            {
                var partition = other._partitions[i];
                _partitions.Add(partition is not null ? (DataSet)partition.Copy(deep: true) : null);
            }
        }
    }

    /// <inheritdoc />
    public override void CopyMetaFrom(DataObject source, bool deep = true)
    {
        // No additional metadata for PartitionedDataSet.
    }

    // ---------------------------------------------------------------
    //  GetDataRange
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets the minimum and maximum values of the named array across all partitions.
    /// </summary>
    /// <param name="name">
    /// Name of the array. When <c>null</c>, an implementation-defined default is used.
    /// </param>
    /// <param name="preference">Preferred field association.</param>
    /// <returns>A (Min, Max) tuple.</returns>
    public override (double Min, double Max) GetDataRange(
        string? name = null,
        FieldAssociation preference = FieldAssociation.Point)
    {
        if (_partitions.Count == 0)
        {
            return (double.NaN, double.NaN);
        }

        double min = double.PositiveInfinity;
        double max = double.NegativeInfinity;

        foreach (var partition in _partitions)
        {
            if (partition is null)
            {
                continue;
            }

            var (partMin, partMax) = partition.GetDataRange(name, preference);
            if (!double.IsNaN(partMin) && partMin < min) min = partMin;
            if (!double.IsNaN(partMax) && partMax > max) max = partMax;
        }

        return (min, max);
    }

    // ---------------------------------------------------------------
    //  Equality
    // ---------------------------------------------------------------

    /// <summary>
    /// Determines whether two <see cref="PartitionedDataSet"/> instances are equal.
    /// Two partitioned datasets are equal if they contain the same number of
    /// partitions and corresponding partitions are equal.
    /// </summary>
    /// <param name="other">The other data object.</param>
    /// <returns><c>true</c> if equal.</returns>
    public new bool Equals(DataObject? other)
    {
        if (other is not PartitionedDataSet otherPds)
        {
            return false;
        }

        if (ReferenceEquals(this, otherPds))
        {
            return true;
        }

        if (_partitions.Count != otherPds._partitions.Count)
        {
            return false;
        }

        for (int i = 0; i < _partitions.Count; i++)
        {
            var a = _partitions[i];
            var b = otherPds._partitions[i];
            if (a is null && b is null) continue;
            if (a is null || b is null) return false;
            if (!a.Equals(b)) return false;
        }

        return true;
    }

    /// <inheritdoc />
    public override bool Equals(object? obj) => Equals(obj as DataObject);

    /// <inheritdoc />
    public override int GetHashCode() => base.GetHashCode();

    // ---------------------------------------------------------------
    //  IList<DataSet?> implementation
    // ---------------------------------------------------------------

    /// <inheritdoc />
    public int IndexOf(DataSet? item) => _partitions.IndexOf(item);

    /// <inheritdoc />
    public bool Contains(DataSet? item) => _partitions.Contains(item);

    /// <inheritdoc />
    void ICollection<DataSet?>.Add(DataSet? item) => Append(item);

    /// <inheritdoc />
    public void CopyTo(DataSet?[] array, int arrayIndex)
    {
        _partitions.CopyTo(array, arrayIndex);
    }

    /// <summary>
    /// Returns an enumerator that iterates through the partitions.
    /// </summary>
    /// <returns>An enumerator of partitions.</returns>
    public IEnumerator<DataSet?> GetEnumerator() => _partitions.GetEnumerator();

    /// <inheritdoc />
    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

    // ---------------------------------------------------------------
    //  String representation
    // ---------------------------------------------------------------

    /// <inheritdoc />
    protected override List<(string Name, string Value)> GetAttributes()
    {
        return new List<(string Name, string Value)>
        {
            ("N Partitions", _partitions.Count.ToString()),
        };
    }

    /// <summary>
    /// Returns a console-friendly representation of this <see cref="PartitionedDataSet"/>.
    /// </summary>
    /// <returns>A formatted string.</returns>
    public override string ToString()
    {
        var sb = new StringBuilder();
        sb.AppendLine($"PartitionedDataSet (0x{GetHashCode():X})");

        var attrs = GetAttributes();
        int maxLen = 0;
        foreach (var a in attrs)
        {
            if (a.Name.Length + 1 > maxLen) maxLen = a.Name.Length + 1;
        }

        maxLen += 3;

        foreach (var a in attrs)
        {
            sb.AppendLine($"  {(a.Name + ":").PadRight(maxLen)}{a.Value}");
        }

        return sb.ToString().TrimEnd();
    }

    // ---------------------------------------------------------------
    //  Private helpers
    // ---------------------------------------------------------------

    /// <summary>
    /// Resolves a possibly negative index to a valid zero-based index.
    /// </summary>
    private int ResolveIndex(int index)
    {
        int count = _partitions.Count;
        if (count == 0)
        {
            throw new IndexOutOfRangeException($"index ({index}) out of range for this dataset.");
        }

        if (index < 0)
        {
            index += count;
        }

        if (index < 0 || index >= count)
        {
            throw new IndexOutOfRangeException($"index ({index}) out of range for this dataset.");
        }

        return index;
    }
}
