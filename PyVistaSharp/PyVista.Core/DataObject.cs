using System.Text;

namespace PyVista.Core;

/// <summary>
/// Abstract base class for all PyVista data objects.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.DataObject</c> class.
/// It provides common functionality shared by all wrapped data objects, including
/// field data management, copy operations, and string representations.
/// </para>
/// <para>
/// Since VTK C# bindings are not available, data is stored internally using
/// <see cref="Dictionary{TKey, TValue}"/> collections.
/// </para>
/// </summary>
public abstract class DataObject : IDisposable, IEquatable<DataObject>
{
    /// <summary>Key used to store the default vector array name in field data.</summary>
    public const string DefaultVectorKey = "_vectors";

    /// <summary>Key used to store the user dictionary as a serialized JSON string.</summary>
    public const string UserDictKey = "_PYVISTA_USER_DICT";

    private readonly Dictionary<string, double[]> _fieldDataArrays = new();
    private bool _disposed;

    /// <summary>
    /// Initializes a new instance of the <see cref="DataObject"/> class.
    /// </summary>
    protected DataObject()
    {
    }

    /// <summary>
    /// Gets the field data associated with this data object.
    /// <para>
    /// Use field data when the size of the data you wish to associate with the
    /// dataset does not match the number of points or cells.
    /// </para>
    /// </summary>
    /// <returns>A <see cref="DataSetAttributes"/> view over the field data.</returns>
    public DataSetAttributes FieldData => new(_fieldDataArrays, FieldAssociation.None);

    /// <summary>
    /// Gets a value indicating whether this data object contains no data.
    /// </summary>
    public abstract bool IsEmpty { get; }

    /// <summary>
    /// Gets the memory address of this managed object, formatted as a hex string.
    /// <para>
    /// Because this is a managed (.NET) object rather than a native VTK object, the
    /// address is derived from <see cref="object.GetHashCode"/>.
    /// </para>
    /// </summary>
    public string MemoryAddress => $"Addr=0x{GetHashCode():X}";

    /// <summary>
    /// Gets the approximate memory size of this object in kibibytes (1024 bytes).
    /// <para>
    /// The estimate accounts for all field data arrays currently stored.
    /// </para>
    /// </summary>
    public long ActualMemorySize
    {
        get
        {
            long bytes = 0;
            foreach (var array in _fieldDataArrays.Values)
            {
                bytes += array.Length * sizeof(double);
            }
            // Return size in kibibytes, with a minimum of 1 when data exists.
            return bytes > 0 ? (bytes / 1024 is 0 ? 1 : bytes / 1024) : 0;
        }
    }

    /// <summary>
    /// Performs a shallow copy from the specified source data object.
    /// <para>
    /// After this call the current instance shares the same array references
    /// as <paramref name="source"/>.
    /// </para>
    /// </summary>
    /// <param name="source">The data object to copy from.</param>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="source"/> is <c>null</c>.
    /// </exception>
    public virtual void ShallowCopy(DataObject source)
    {
        ArgumentNullException.ThrowIfNull(source);
        _fieldDataArrays.Clear();
        foreach (var kvp in source._fieldDataArrays)
        {
            _fieldDataArrays[kvp.Key] = kvp.Value;
        }
    }

    /// <summary>
    /// Performs a deep copy from the specified source data object.
    /// <para>
    /// All arrays are cloned so the two objects are fully independent.
    /// </para>
    /// </summary>
    /// <param name="source">The data object to copy from.</param>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="source"/> is <c>null</c>.
    /// </exception>
    public virtual void DeepCopy(DataObject source)
    {
        ArgumentNullException.ThrowIfNull(source);
        _fieldDataArrays.Clear();
        foreach (var kvp in source._fieldDataArrays)
        {
            _fieldDataArrays[kvp.Key] = (double[])kvp.Value.Clone();
        }
    }

    /// <summary>
    /// Returns a copy of the object.
    /// </summary>
    /// <param name="deep">
    /// When <c>true</c> (default), performs a full deep copy.
    /// When <c>false</c>, performs a shallow copy where data arrays are shared references.
    /// </param>
    /// <returns>A new <see cref="DataObject"/> of the same runtime type.</returns>
    /// <exception cref="InvalidOperationException">
    /// Thrown when the derived type cannot be instantiated.
    /// </exception>
    public DataObject Copy(bool deep = true)
    {
        var clone = (DataObject)(Activator.CreateInstance(GetType())
            ?? throw new InvalidOperationException(
                $"Unable to create an instance of {GetType().Name}. Ensure the type has a public parameterless constructor."));

        if (deep)
        {
            clone.DeepCopy(this);
        }
        else
        {
            clone.ShallowCopy(this);
        }

        clone.CopyMetaFrom(this, deep);
        return clone;
    }

    /// <summary>
    /// Copies PyVista metadata from another data object onto this instance.
    /// <para>
    /// Intended to be overridden by subclasses that carry additional metadata.
    /// </para>
    /// </summary>
    /// <param name="source">The data object to copy metadata from.</param>
    /// <param name="deep">Whether to perform a deep copy of metadata.</param>
    public virtual void CopyMetaFrom(DataObject source, bool deep = true)
    {
        // Base implementation is intentionally empty.
        // Subclasses override to copy type-specific metadata.
    }

    /// <summary>
    /// Adds a field data array.
    /// <para>
    /// Use field data when the size of the data does not match the number of
    /// points or cells of the dataset.
    /// </para>
    /// </summary>
    /// <param name="array">The data to add.</param>
    /// <param name="name">The name to assign to the field array.</param>
    /// <param name="deep">
    /// When <c>true</c> (default), a deep copy of the data is stored.
    /// </param>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="array"/> or <paramref name="name"/> is <c>null</c>.
    /// </exception>
    public void AddFieldData(double[] array, string name, bool deep = true)
    {
        var data = deep ? (double[])array.Clone() : array;
        FieldData.SetArray(data, name);
    }

    /// <summary>
    /// Removes all field data arrays.
    /// </summary>
    public void ClearFieldData()
    {
        FieldData.Clear();
    }

    /// <summary>
    /// Copies the structure (geometry and topology) from the specified dataset.
    /// <para>
    /// The base implementation copies field data arrays. Subclasses should override
    /// to copy type-specific structural information (e.g., points, cells).
    /// </para>
    /// </summary>
    /// <param name="dataset">The data object to copy structure from.</param>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="dataset"/> is <c>null</c>.
    /// </exception>
    public virtual void CopyStructure(DataObject dataset)
    {
        ArgumentNullException.ThrowIfNull(dataset);
        if (ReferenceEquals(this, dataset))
        {
            return;
        }

        _fieldDataArrays.Clear();
        foreach (var kvp in dataset._fieldDataArrays)
        {
            _fieldDataArrays[kvp.Key] = (double[])kvp.Value.Clone();
        }
    }

    /// <summary>
    /// Copies all data attributes from the specified dataset.
    /// <para>
    /// The base implementation copies field data arrays. Subclasses should override
    /// to copy additional attribute types (e.g., point data, cell data).
    /// </para>
    /// </summary>
    /// <param name="dataset">The data object to copy attributes from.</param>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="dataset"/> is <c>null</c>.
    /// </exception>
    public virtual void CopyAttributes(DataObject dataset)
    {
        ArgumentNullException.ThrowIfNull(dataset);
        _fieldDataArrays.Clear();
        foreach (var kvp in dataset._fieldDataArrays)
        {
            _fieldDataArrays[kvp.Key] = (double[])kvp.Value.Clone();
        }
    }

    /// <summary>
    /// Converts this data object to a <see cref="MultiBlock"/> containing a single block.
    /// </summary>
    /// <returns>A <see cref="MultiBlock"/> wrapping this data object.</returns>
    public MultiBlock CastToMultiBlock()
    {
        var mb = new MultiBlock();
        mb.Append(this);
        return mb;
    }

    /// <summary>
    /// Gets the non-NaN minimum and maximum values of the named array.
    /// </summary>
    /// <param name="name">The name of the array. When <c>null</c>, an implementation-defined default is used.</param>
    /// <param name="preference">The preferred field association to search.</param>
    /// <returns>A tuple of (min, max) values.</returns>
    public abstract (double Min, double Max) GetDataRange(string? name = null, FieldAssociation preference = FieldAssociation.Point);

    /// <summary>
    /// Returns the header statistics of this dataset as a console-friendly string.
    /// </summary>
    /// <returns>A formatted string describing this data object.</returns>
    public string Head()
    {
        var sb = new StringBuilder();
        sb.AppendLine($"{GetType().Name} ({MemoryAddress})");

        var attrs = GetAttributes();
        if (attrs.Count == 0)
        {
            return sb.ToString().TrimEnd();
        }

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

        sb.AppendLine($"  {"N Arrays:".PadRight(maxLen)}{_fieldDataArrays.Count}");

        return sb.ToString().TrimEnd();
    }

    /// <summary>
    /// Returns a console-friendly representation of this data object,
    /// equivalent to <see cref="Head"/>.
    /// </summary>
    /// <returns>A formatted string describing this data object.</returns>
    public override string ToString() => Head();

    /// <inheritdoc />
    public bool Equals(DataObject? other)
    {
        if (other is null) return false;
        if (ReferenceEquals(this, other)) return true;
        if (GetType() != other.GetType()) return false;

        if (_fieldDataArrays.Count != other._fieldDataArrays.Count) return false;

        foreach (var kvp in _fieldDataArrays)
        {
            if (!other._fieldDataArrays.TryGetValue(kvp.Key, out var otherArray))
            {
                return false;
            }
            if (!kvp.Value.AsSpan().SequenceEqual(otherArray))
            {
                return false;
            }
        }

        return true;
    }

    /// <inheritdoc />
    public override bool Equals(object? obj) => Equals(obj as DataObject);

    /// <inheritdoc />
    public override int GetHashCode()
    {
        // Mutable object – return a constant to satisfy the contract without
        // violating the rule that equal objects must have equal hash codes.
        return GetType().GetHashCode();
    }

    /// <summary>
    /// Returns descriptive attribute tuples for the <see cref="Head"/> representation.
    /// <para>
    /// Subclasses should override to supply type-specific attributes
    /// (e.g., number of points, number of cells, bounds).
    /// </para>
    /// </summary>
    /// <returns>A list of name/value pairs describing this object.</returns>
    protected virtual List<(string Name, string Value)> GetAttributes() => new();

    /// <inheritdoc />
    public void Dispose()
    {
        Dispose(disposing: true);
        GC.SuppressFinalize(this);
    }

    /// <summary>
    /// Releases the resources used by this <see cref="DataObject"/>.
    /// </summary>
    /// <param name="disposing">
    /// <c>true</c> to release both managed and unmanaged resources;
    /// <c>false</c> to release only unmanaged resources.
    /// </param>
    protected virtual void Dispose(bool disposing)
    {
        if (!_disposed)
        {
            if (disposing)
            {
                _fieldDataArrays.Clear();
            }
            _disposed = true;
        }
    }
}

/// <summary>
/// A composite dataset that stores multiple <see cref="DataObject"/> instances as blocks.
/// <para>
/// This is a minimal analogue of the Python <c>pyvista.MultiBlock</c> class, provided so
/// that <see cref="DataObject.CastToMultiBlock"/> can return a concrete type.
/// </para>
/// </summary>
public sealed class MultiBlock : DataObject
{
    private readonly List<DataObject> _blocks = new();

    /// <summary>
    /// Initializes a new instance of the <see cref="MultiBlock"/> class.
    /// </summary>
    public MultiBlock()
    {
    }

    /// <inheritdoc />
    public override bool IsEmpty => _blocks.Count == 0;

    /// <summary>Gets the number of blocks.</summary>
    public int Count => _blocks.Count;

    /// <summary>
    /// Gets or sets the block at the specified <paramref name="index"/>.
    /// </summary>
    /// <param name="index">The zero-based index of the block.</param>
    /// <returns>The data object at the specified index.</returns>
    public DataObject this[int index]
    {
        get => _blocks[index];
        set => _blocks[index] = value ?? throw new ArgumentNullException(nameof(value));
    }

    /// <summary>
    /// Appends a data object as a new block.
    /// </summary>
    /// <param name="block">The data object to append.</param>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="block"/> is <c>null</c>.
    /// </exception>
    public void Append(DataObject block)
    {
        ArgumentNullException.ThrowIfNull(block);
        _blocks.Add(block);
    }

    /// <inheritdoc />
    public override (double Min, double Max) GetDataRange(string? name = null, FieldAssociation preference = FieldAssociation.Point)
    {
        if (_blocks.Count == 0)
        {
            return (double.NaN, double.NaN);
        }

        double min = double.MaxValue;
        double max = double.MinValue;

        foreach (var block in _blocks)
        {
            var (blockMin, blockMax) = block.GetDataRange(name, preference);
            if (blockMin < min) min = blockMin;
            if (blockMax > max) max = blockMax;
        }

        return (min, max);
    }

    /// <inheritdoc />
    protected override List<(string Name, string Value)> GetAttributes()
    {
        return new List<(string Name, string Value)>
        {
            ("N Blocks", _blocks.Count.ToString()),
        };
    }
}
