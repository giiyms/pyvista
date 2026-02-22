using System.Collections;
using System.Globalization;
using System.Text;

namespace PyVista.Core;

/// <summary>
/// A composite dataset that stores multiple <see cref="DataObject"/> instances as named blocks.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.MultiBlock</c> class.
/// It wraps/extends the concept of a VTK <c>vtkMultiBlockDataSet</c> so that
/// multiple datasets can be stored together and accessed by index or string name.
/// </para>
/// <para>
/// You can think of <see cref="MultiBlock"/> like a list – it can be iterated by index –
/// but it also has dictionary-like features since blocks can be accessed by their
/// string name.
/// </para>
/// </summary>
public class MultiBlock : DataObject, IList<DataObject?>
{
    private readonly List<DataObject?> _blocks = new();
    private readonly List<string> _names = new();

    // ---------------------------------------------------------------
    //  Constructors
    // ---------------------------------------------------------------

    /// <summary>
    /// Initializes a new, empty <see cref="MultiBlock"/>.
    /// </summary>
    public MultiBlock()
    {
    }

    /// <summary>
    /// Initializes a new <see cref="MultiBlock"/> from a list of datasets.
    /// </summary>
    /// <param name="datasets">The datasets to include as blocks.</param>
    public MultiBlock(IEnumerable<DataObject?> datasets)
    {
        ArgumentNullException.ThrowIfNull(datasets);
        foreach (var ds in datasets)
        {
            Append(ds);
        }
    }

    /// <summary>
    /// Initializes a new <see cref="MultiBlock"/> from a dictionary of named datasets.
    /// </summary>
    /// <param name="namedDatasets">A mapping of block names to datasets.</param>
    public MultiBlock(IDictionary<string, DataObject?> namedDatasets)
    {
        ArgumentNullException.ThrowIfNull(namedDatasets);
        foreach (var kvp in namedDatasets)
        {
            Append(kvp.Value, kvp.Key);
        }
    }

    /// <summary>
    /// Initializes a new <see cref="MultiBlock"/> as a copy of another <see cref="MultiBlock"/>.
    /// </summary>
    /// <param name="other">The source multiblock to copy.</param>
    /// <param name="deep">When <c>true</c>, performs a deep copy of all blocks.</param>
    public MultiBlock(MultiBlock other, bool deep = true)
    {
        ArgumentNullException.ThrowIfNull(other);
        if (deep)
        {
            DeepCopy(other);
        }
        else
        {
            ShallowCopy(other);
        }
    }

    // ---------------------------------------------------------------
    //  Block count / emptiness
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets or sets the number of blocks.
    /// <para>
    /// Setting a value larger than the current count pads with <c>null</c> blocks.
    /// Setting a smaller value truncates.
    /// </para>
    /// </summary>
    public int NBlocks
    {
        get => _blocks.Count;
        set
        {
            if (value < 0)
            {
                throw new ArgumentOutOfRangeException(nameof(value), "NBlocks cannot be negative.");
            }

            while (_blocks.Count < value)
            {
                _blocks.Add(null);
                _names.Add(DefaultBlockName(_blocks.Count - 1));
            }

            while (_blocks.Count > value)
            {
                _blocks.RemoveAt(_blocks.Count - 1);
                _names.RemoveAt(_names.Count - 1);
            }
        }
    }

    /// <inheritdoc />
    public override bool IsEmpty => _blocks.Count == 0;

    /// <summary>
    /// Gets the number of blocks. Equivalent to <see cref="NBlocks"/>.
    /// </summary>
    public int Count => _blocks.Count;

    /// <inheritdoc />
    bool ICollection<DataObject?>.IsReadOnly => false;

    // ---------------------------------------------------------------
    //  Nesting queries
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets a value indicating whether any of the blocks are themselves a <see cref="MultiBlock"/>.
    /// </summary>
    public bool IsNested
    {
        get
        {
            foreach (var block in _blocks)
            {
                if (block is MultiBlock)
                {
                    return true;
                }
            }

            return false;
        }
    }

    // ---------------------------------------------------------------
    //  Bounds / Center / Volume / Length
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets the axis-aligned bounding box that encloses all blocks.
    /// </summary>
    public BoundsTuple Bounds
    {
        get
        {
            double xMin = double.MaxValue, xMax = double.MinValue;
            double yMin = double.MaxValue, yMax = double.MinValue;
            double zMin = double.MaxValue, zMax = double.MinValue;
            bool found = false;

            foreach (var block in _blocks)
            {
                if (block is DataSet ds && ds.NPoints > 0)
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
                else if (block is MultiBlock mb && !mb.IsEmpty)
                {
                    var b = mb.Bounds;
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
    /// Gets the total volume of all mesh blocks.
    /// </summary>
    public double Volume
    {
        get
        {
            double total = 0;
            foreach (var block in _blocks)
            {
                if (block is MultiBlock mb)
                {
                    total += mb.Volume;
                }
                else if (block is PointSet ps)
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
    //  Indexers
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets or sets the block at the specified integer index.
    /// <para>
    /// Negative indices are supported (e.g., <c>-1</c> for the last block).
    /// When setting, the block name defaults to <c>"Block-{i:00}"</c> unless
    /// previously set.
    /// </para>
    /// </summary>
    /// <param name="index">Zero-based index (negative indices count from the end).</param>
    /// <returns>The block at the specified index, which may be <c>null</c>.</returns>
    /// <exception cref="IndexOutOfRangeException">Thrown when the index is out of range.</exception>
    public DataObject? this[int index]
    {
        get
        {
            int resolved = ResolveIndex(index);
            return _blocks[resolved];
        }
        set
        {
            int resolved = ResolveIndex(index);
            _blocks[resolved] = value;
        }
    }

    /// <summary>
    /// Gets or sets the block with the specified name.
    /// <para>
    /// If the name is not found during a set, a new block is appended.
    /// </para>
    /// </summary>
    /// <param name="name">Name of the block.</param>
    /// <returns>The block with the specified name, which may be <c>null</c>.</returns>
    /// <exception cref="KeyNotFoundException">Thrown when the name is not found during a get.</exception>
    public DataObject? this[string name]
    {
        get
        {
            int idx = GetIndexByName(name);
            return _blocks[idx];
        }
        set
        {
            int idx = IndexOfName(name);
            if (idx < 0)
            {
                Append(value, name);
            }
            else
            {
                _blocks[idx] = value;
            }
        }
    }

    // ---------------------------------------------------------------
    //  Name management
    // ---------------------------------------------------------------

    /// <summary>
    /// Sets the name of the block at the specified index.
    /// </summary>
    /// <param name="index">Zero-based block index.</param>
    /// <param name="name">
    /// Name to assign. If <c>null</c>, the name is not changed.
    /// </param>
    public void SetBlockName(int index, string? name)
    {
        if (name is null)
        {
            return;
        }

        int resolved = ResolveIndex(index);
        _names[resolved] = name;
    }

    /// <summary>
    /// Gets the name of the block at the specified index.
    /// </summary>
    /// <param name="index">Zero-based block index (negative indices supported).</param>
    /// <returns>The name of the block.</returns>
    public string GetBlockName(int index)
    {
        int resolved = ResolveIndex(index);
        return _names[resolved];
    }

    /// <summary>
    /// Gets all block names.
    /// </summary>
    /// <returns>A list of block names.</returns>
    public List<string> Keys()
    {
        return new List<string>(_names);
    }

    /// <summary>
    /// Finds the index of the first block with the specified name.
    /// </summary>
    /// <param name="name">The name to search for.</param>
    /// <returns>The zero-based index of the block.</returns>
    /// <exception cref="KeyNotFoundException">Thrown when no block has the specified name.</exception>
    public int GetIndexByName(string name)
    {
        int idx = IndexOfName(name);
        if (idx < 0)
        {
            throw new KeyNotFoundException($"Block name ({name}) not found");
        }

        return idx;
    }

    // ---------------------------------------------------------------
    //  Append / Extend / Insert / Pop / Remove
    // ---------------------------------------------------------------

    /// <summary>
    /// Appends a dataset as a new block.
    /// </summary>
    /// <param name="dataset">Dataset to append (may be <c>null</c>).</param>
    /// <param name="name">
    /// Optional name. When <c>null</c>, a default name <c>"Block-{i:00}"</c> is assigned.
    /// </param>
    /// <exception cref="InvalidOperationException">
    /// Thrown when attempting to nest a <see cref="MultiBlock"/> inside itself.
    /// </exception>
    public void Append(DataObject? dataset, string? name = null)
    {
        if (ReferenceEquals(dataset, this))
        {
            throw new InvalidOperationException("Cannot nest a composite dataset in itself.");
        }

        int index = _blocks.Count;
        _blocks.Add(dataset);
        _names.Add(name ?? DefaultBlockName(index));
    }

    /// <summary>
    /// Extends this <see cref="MultiBlock"/> with the blocks from an iterable.
    /// <para>
    /// If the iterable is itself a <see cref="MultiBlock"/>, the block names are preserved.
    /// </para>
    /// </summary>
    /// <param name="datasets">The datasets to append.</param>
    public void Extend(IEnumerable<DataObject?> datasets)
    {
        ArgumentNullException.ThrowIfNull(datasets);
        if (datasets is MultiBlock mb)
        {
            for (int i = 0; i < mb.NBlocks; i++)
            {
                Append(mb[i], mb.GetBlockName(i));
            }
        }
        else
        {
            foreach (var ds in datasets)
            {
                Append(ds);
            }
        }
    }

    /// <summary>
    /// Inserts a dataset at the given index.
    /// </summary>
    /// <param name="index">Zero-based insertion index.</param>
    /// <param name="dataset">Dataset to insert (may be <c>null</c>).</param>
    /// <param name="name">
    /// Optional name. When <c>null</c>, a default name is assigned.
    /// </param>
    public void Insert(int index, DataObject? dataset, string? name = null)
    {
        int resolved = index < 0 ? Math.Max(_blocks.Count + index, 0) : Math.Min(index, _blocks.Count);
        _blocks.Insert(resolved, dataset);
        _names.Insert(resolved, name ?? DefaultBlockName(resolved));
    }

    /// <inheritdoc />
    void IList<DataObject?>.Insert(int index, DataObject? item) => Insert(index, item);

    /// <summary>
    /// Removes and returns the block at the specified index or name.
    /// </summary>
    /// <param name="index">Zero-based index (negative indices supported). Defaults to <c>-1</c> (last).</param>
    /// <returns>The removed block.</returns>
    public DataObject? Pop(int index = -1)
    {
        int resolved = ResolveIndex(index);
        var data = _blocks[resolved];
        _blocks.RemoveAt(resolved);
        _names.RemoveAt(resolved);
        return data;
    }

    /// <summary>
    /// Removes and returns the block with the specified name.
    /// </summary>
    /// <param name="name">Name of the block to remove.</param>
    /// <returns>The removed block.</returns>
    /// <exception cref="KeyNotFoundException">Thrown when no block has the specified name.</exception>
    public DataObject? Pop(string name)
    {
        int idx = GetIndexByName(name);
        return Pop(idx);
    }

    /// <summary>
    /// Removes the block at the specified index.
    /// </summary>
    /// <param name="index">Zero-based index (negative indices supported).</param>
    public void RemoveAt(int index)
    {
        int resolved = ResolveIndex(index);
        _blocks.RemoveAt(resolved);
        _names.RemoveAt(resolved);
    }

    /// <summary>
    /// Removes the first occurrence of the specified block.
    /// </summary>
    /// <param name="item">The block to remove.</param>
    /// <returns><c>true</c> if the block was found and removed.</returns>
    public bool Remove(DataObject? item)
    {
        int idx = _blocks.IndexOf(item);
        if (idx < 0)
        {
            return false;
        }

        _blocks.RemoveAt(idx);
        _names.RemoveAt(idx);
        return true;
    }

    /// <summary>
    /// Removes the block with the specified name.
    /// </summary>
    /// <param name="name">Name of the block to remove.</param>
    /// <exception cref="KeyNotFoundException">Thrown when no block has the specified name.</exception>
    public void RemoveByName(string name)
    {
        int idx = GetIndexByName(name);
        _blocks.RemoveAt(idx);
        _names.RemoveAt(idx);
    }

    /// <summary>
    /// Replaces the block at the specified index while preserving the block name.
    /// </summary>
    /// <param name="index">Zero-based index of the block to replace.</param>
    /// <param name="dataset">The replacement dataset.</param>
    public void Replace(int index, DataObject? dataset)
    {
        int resolved = ResolveIndex(index);
        _blocks[resolved] = dataset;
    }

    /// <summary>
    /// Replaces the block with the specified name while preserving the name.
    /// </summary>
    /// <param name="name">Name of the block to replace.</param>
    /// <param name="dataset">The replacement dataset.</param>
    public void Replace(string name, DataObject? dataset)
    {
        int idx = GetIndexByName(name);
        _blocks[idx] = dataset;
    }

    // ---------------------------------------------------------------
    //  Clean / Reverse
    // ---------------------------------------------------------------

    /// <summary>
    /// Removes all <c>null</c> blocks in place. Optionally also removes empty meshes.
    /// </summary>
    /// <param name="empty">
    /// When <c>true</c> (default), also removes datasets with zero points.
    /// </param>
    public void Clean(bool empty = true)
    {
        for (int i = _blocks.Count - 1; i >= 0; i--)
        {
            var data = _blocks[i];
            if (data is MultiBlock nested)
            {
                nested.Clean(empty);
                if (nested.NBlocks == 0)
                {
                    _blocks.RemoveAt(i);
                    _names.RemoveAt(i);
                }
            }
            else if (data is null)
            {
                _blocks.RemoveAt(i);
                _names.RemoveAt(i);
            }
            else if (empty && data is DataSet ds && ds.NPoints == 0)
            {
                _blocks.RemoveAt(i);
                _names.RemoveAt(i);
            }
        }
    }

    /// <summary>
    /// Reverses the order of blocks in place.
    /// </summary>
    public void Reverse()
    {
        _blocks.Reverse();
        _names.Reverse();
    }

    /// <summary>
    /// Removes all blocks.
    /// </summary>
    public void Clear()
    {
        _blocks.Clear();
        _names.Clear();
    }

    // ---------------------------------------------------------------
    //  Type queries
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets a value indicating whether all recursively nested leaf blocks are <see cref="PolyData"/>.
    /// </summary>
    public bool IsAllPolyData
    {
        get
        {
            foreach (var block in _blocks)
            {
                if (block is MultiBlock mb)
                {
                    if (!mb.IsAllPolyData)
                    {
                        return false;
                    }
                }
                else if (block is not PolyData)
                {
                    return false;
                }
            }

            return _blocks.Count > 0;
        }
    }

    /// <summary>
    /// Gets a value indicating whether all recursively nested leaf blocks have the same runtime type.
    /// </summary>
    public bool IsHomogeneous
    {
        get
        {
            var types = GetNestedBlockTypes();
            return types.Count == 1;
        }
    }

    /// <summary>
    /// Gets a value indicating whether any two recursively nested leaf blocks have different runtime types.
    /// </summary>
    public bool IsHeterogeneous
    {
        get
        {
            var types = GetNestedBlockTypes();
            return types.Count > 1;
        }
    }

    /// <summary>
    /// Returns the set of all runtime types across all recursively nested leaf blocks.
    /// </summary>
    /// <returns>A set of <see cref="Type"/> instances.</returns>
    public HashSet<Type> GetNestedBlockTypes()
    {
        var types = new HashSet<Type>();
        CollectNestedTypes(types);
        return types;
    }

    // ---------------------------------------------------------------
    //  Conversion helpers
    // ---------------------------------------------------------------

    /// <summary>
    /// Returns a new <see cref="MultiBlock"/> where all leaf blocks are <see cref="PolyData"/>.
    /// <para>
    /// Blocks that are already <see cref="PolyData"/> are reused as-is (unless <paramref name="copy"/>
    /// is <c>true</c>). <c>null</c> blocks are converted to empty <see cref="PolyData"/> instances.
    /// Other <see cref="DataSet"/> blocks are not converted (stored as-is) because VTK surface
    /// extraction is not available in pure managed code.
    /// </para>
    /// </summary>
    /// <param name="copy">When <c>true</c>, copies blocks that are already <see cref="PolyData"/>.</param>
    /// <returns>A new <see cref="MultiBlock"/> containing <see cref="PolyData"/> blocks where possible.</returns>
    public MultiBlock AsPolyDataBlocks(bool copy = false)
    {
        var result = new MultiBlock();
        for (int i = 0; i < _blocks.Count; i++)
        {
            var block = _blocks[i];
            DataObject? converted;
            if (block is null)
            {
                converted = new PolyData();
            }
            else if (block is MultiBlock mb)
            {
                converted = mb.AsPolyDataBlocks(copy);
            }
            else if (block is PolyData pd)
            {
                converted = copy ? (PolyData)pd.Copy(deep: true) : pd;
            }
            else
            {
                // Without VTK, we cannot extract surface; store the block as-is.
                converted = block;
            }

            result.Append(converted, _names[i]);
        }

        return result;
    }

    // ---------------------------------------------------------------
    //  Copy
    // ---------------------------------------------------------------

    /// <summary>
    /// Returns a copy of this <see cref="MultiBlock"/>.
    /// </summary>
    /// <param name="deep">When <c>true</c> (default), performs a deep copy of all blocks.</param>
    /// <returns>A new <see cref="MultiBlock"/>.</returns>
    public new MultiBlock Copy(bool deep = true)
    {
        return new MultiBlock(this, deep);
    }

    /// <inheritdoc />
    public override void ShallowCopy(DataObject source)
    {
        base.ShallowCopy(source);
        if (source is MultiBlock other)
        {
            _blocks.Clear();
            _names.Clear();
            _blocks.AddRange(other._blocks);
            _names.AddRange(other._names);
        }
    }

    /// <inheritdoc />
    public override void DeepCopy(DataObject source)
    {
        base.DeepCopy(source);
        if (source is MultiBlock other)
        {
            _blocks.Clear();
            _names.Clear();
            for (int i = 0; i < other._blocks.Count; i++)
            {
                var block = other._blocks[i];
                _blocks.Add(block?.Copy(deep: true));
                _names.Add(other._names[i]);
            }
        }
    }

    // ---------------------------------------------------------------
    //  GetDataRange
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets the minimum and maximum values of the named array across all blocks.
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
        if (_blocks.Count == 0)
        {
            return (double.NaN, double.NaN);
        }

        double min = double.PositiveInfinity;
        double max = double.NegativeInfinity;

        foreach (var block in _blocks)
        {
            if (block is null)
            {
                continue;
            }

            var (blockMin, blockMax) = block.GetDataRange(name, preference);
            if (!double.IsNaN(blockMin) && blockMin < min) min = blockMin;
            if (!double.IsNaN(blockMax) && blockMax > max) max = blockMax;
        }

        return (min, max);
    }

    // ---------------------------------------------------------------
    //  Equality
    // ---------------------------------------------------------------

    /// <summary>
    /// Determines whether two <see cref="MultiBlock"/> instances are equal.
    /// Two multiblocks are equal if they contain the same number of blocks,
    /// have the same names, and corresponding blocks are equal.
    /// </summary>
    /// <param name="other">The other data object.</param>
    /// <returns><c>true</c> if equal.</returns>
    public new bool Equals(DataObject? other)
    {
        if (other is not MultiBlock otherMb)
        {
            return false;
        }

        if (ReferenceEquals(this, otherMb))
        {
            return true;
        }

        if (_blocks.Count != otherMb._blocks.Count)
        {
            return false;
        }

        for (int i = 0; i < _blocks.Count; i++)
        {
            if (_names[i] != otherMb._names[i])
            {
                return false;
            }

            var a = _blocks[i];
            var b = otherMb._blocks[i];
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
    //  IList<DataObject?> implementation
    // ---------------------------------------------------------------

    /// <inheritdoc />
    public int IndexOf(DataObject? item) => _blocks.IndexOf(item);

    /// <inheritdoc />
    public bool Contains(DataObject? item) => _blocks.Contains(item);

    /// <inheritdoc />
    void ICollection<DataObject?>.Add(DataObject? item) => Append(item);

    /// <inheritdoc />
    public void CopyTo(DataObject?[] array, int arrayIndex)
    {
        _blocks.CopyTo(array, arrayIndex);
    }

    /// <summary>
    /// Returns an enumerator that iterates through the blocks.
    /// </summary>
    /// <returns>An enumerator of blocks.</returns>
    public IEnumerator<DataObject?> GetEnumerator() => _blocks.GetEnumerator();

    /// <inheritdoc />
    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

    // ---------------------------------------------------------------
    //  String representation
    // ---------------------------------------------------------------

    /// <inheritdoc />
    protected override List<(string Name, string Value)> GetAttributes()
    {
        var attrs = new List<(string Name, string Value)>
        {
            ("N Blocks", _blocks.Count.ToString()),
        };

        if (_blocks.Count > 0)
        {
            var b = Bounds;
            attrs.Add(("X Bounds", string.Format(CultureInfo.InvariantCulture,
                "{0:E3}, {1:E3}", b.XMin, b.XMax)));
            attrs.Add(("Y Bounds", string.Format(CultureInfo.InvariantCulture,
                "{0:E3}, {1:E3}", b.YMin, b.YMax)));
            attrs.Add(("Z Bounds", string.Format(CultureInfo.InvariantCulture,
                "{0:E3}, {1:E3}", b.ZMin, b.ZMax)));
        }

        return attrs;
    }

    /// <summary>
    /// Returns a console-friendly representation of this <see cref="MultiBlock"/>.
    /// </summary>
    /// <returns>A formatted string.</returns>
    public override string ToString()
    {
        var sb = new StringBuilder();
        sb.AppendLine($"MultiBlock (0x{GetHashCode():X})");

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
        int count = _blocks.Count;
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

    /// <summary>
    /// Returns the index of the first block with the given name, or <c>-1</c> if not found.
    /// </summary>
    private int IndexOfName(string name)
    {
        for (int i = 0; i < _names.Count; i++)
        {
            if (_names[i] == name)
            {
                return i;
            }
        }

        return -1;
    }

    /// <summary>
    /// Generates a default block name for the given index.
    /// </summary>
    private static string DefaultBlockName(int index) => $"Block-{index:D2}";

    /// <summary>
    /// Recursively collects leaf block types into the provided set.
    /// </summary>
    private void CollectNestedTypes(HashSet<Type> types)
    {
        foreach (var block in _blocks)
        {
            if (block is MultiBlock mb)
            {
                mb.CollectNestedTypes(types);
            }
            else
            {
                types.Add(block?.GetType() ?? typeof(object));
            }
        }
    }
}
