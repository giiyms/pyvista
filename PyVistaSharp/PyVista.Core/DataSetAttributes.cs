using System.Collections;
using System.Diagnostics.CodeAnalysis;

namespace PyVista.Core;

/// <summary>
/// A dictionary-like container for arrays associated with a dataset.
/// Provides the ability to pick one of the present arrays as the currently active
/// array for each attribute type (scalars, vectors, normals, texture coordinates).
/// </summary>
/// <remarks>
/// This is a C# port of the Python <c>DataSetAttributes</c> class from PyVista.
/// Since there is no VTK dependency, a <see cref="Dictionary{TKey,TValue}"/> is used
/// as internal storage with insertion-order preserved via a separate key list.
/// </remarks>
public class DataSetAttributes : IDictionary<string, double[]>, IEnumerable
{
    private readonly Dictionary<string, double[]> _arrays = new();
    private readonly List<string> _orderedKeys = new();

    private string? _activeScalarsName;
    private string? _activeVectorsName;
    private string? _activeNormalsName;
    private string? _activeTextureCoordinatesName;

    /// <summary>
    /// Initializes a new instance of the <see cref="DataSetAttributes"/> class.
    /// </summary>
    /// <param name="association">The field association type for this attribute set.</param>
    public DataSetAttributes(FieldAssociation association)
    {
        Association = association;
    }

    /// <summary>
    /// Initializes a new instance of the <see cref="DataSetAttributes"/> class
    /// that wraps the supplied backing dictionary.
    /// </summary>
    /// <param name="arrays">The backing dictionary of named arrays.</param>
    /// <param name="association">The field association type for this attribute set.</param>
    internal DataSetAttributes(Dictionary<string, double[]> arrays, FieldAssociation association)
    {
        ArgumentNullException.ThrowIfNull(arrays);
        _arrays = arrays;
        _orderedKeys = new List<string>(arrays.Keys);
        Association = association;
    }

    /// <summary>
    /// Gets the field association for this attribute set.
    /// </summary>
    public FieldAssociation Association { get; }

    /// <summary>
    /// Gets the number of arrays stored in this object.
    /// </summary>
    public int Count => _arrays.Count;

    /// <inheritdoc />
    bool ICollection<KeyValuePair<string, double[]>>.IsReadOnly => false;

    /// <summary>
    /// Gets or sets the active scalars array, or <c>null</c> if no active scalars are set.
    /// </summary>
    /// <exception cref="InvalidOperationException">
    /// Thrown when accessing on <see cref="FieldAssociation.None"/> association.
    /// </exception>
    public double[]? ActiveScalars
    {
        get
        {
            RaiseFieldDataNoScalarsVectorsNormals();
            return _activeScalarsName is not null && _arrays.TryGetValue(_activeScalarsName, out var arr)
                ? arr
                : null;
        }
    }

    /// <summary>
    /// Gets the active vectors array, or <c>null</c> if no active vectors are set.
    /// </summary>
    /// <exception cref="InvalidOperationException">
    /// Thrown when accessing on <see cref="FieldAssociation.None"/> association.
    /// </exception>
    public double[]? ActiveVectors
    {
        get
        {
            RaiseFieldDataNoScalarsVectorsNormals();
            return _activeVectorsName is not null && _arrays.TryGetValue(_activeVectorsName, out var arr)
                ? arr
                : null;
        }
    }

    /// <summary>
    /// Gets the active normals array, or <c>null</c> if no active normals are set.
    /// </summary>
    /// <exception cref="InvalidOperationException">
    /// Thrown when accessing on <see cref="FieldAssociation.None"/> association.
    /// </exception>
    public double[]? ActiveNormals
    {
        get
        {
            RaiseFieldDataNoScalarsVectorsNormals();
            return _activeNormalsName is not null && _arrays.TryGetValue(_activeNormalsName, out var arr)
                ? arr
                : null;
        }
    }

    /// <summary>
    /// Gets the active texture coordinates array, or <c>null</c> if none are set.
    /// </summary>
    /// <exception cref="InvalidOperationException">
    /// Thrown when accessing on <see cref="FieldAssociation.None"/> association.
    /// </exception>
    public double[]? ActiveTextureCoordinates
    {
        get
        {
            RaiseNoTextureCoordinates();
            return _activeTextureCoordinatesName is not null
                && _arrays.TryGetValue(_activeTextureCoordinatesName, out var arr)
                ? arr
                : null;
        }
    }

    /// <summary>
    /// Gets or sets the name of the active scalars array.
    /// </summary>
    /// <exception cref="InvalidOperationException">
    /// Thrown when setting on <see cref="FieldAssociation.None"/> association.
    /// </exception>
    /// <exception cref="KeyNotFoundException">
    /// Thrown when the specified name does not exist in this attribute set.
    /// </exception>
    public string? ActiveScalarsName
    {
        get => _activeScalarsName is not null && _arrays.ContainsKey(_activeScalarsName)
            ? _activeScalarsName
            : null;
        set
        {
            if (value is null)
            {
                _activeScalarsName = null;
                return;
            }

            RaiseFieldDataNoScalarsVectorsNormals();
            if (!_arrays.ContainsKey(value))
            {
                throw new KeyNotFoundException($"DataSetAttributes does not contain \"{value}\"");
            }

            _activeScalarsName = value;
        }
    }

    /// <summary>
    /// Gets or sets the name of the active vectors array.
    /// </summary>
    /// <exception cref="InvalidOperationException">
    /// Thrown when setting on <see cref="FieldAssociation.None"/> association.
    /// </exception>
    /// <exception cref="KeyNotFoundException">
    /// Thrown when the specified name does not exist in this attribute set.
    /// </exception>
    public string? ActiveVectorsName
    {
        get => _activeVectorsName is not null && _arrays.ContainsKey(_activeVectorsName)
            ? _activeVectorsName
            : null;
        set
        {
            if (value is null)
            {
                _activeVectorsName = null;
                return;
            }

            RaiseFieldDataNoScalarsVectorsNormals();
            if (!_arrays.ContainsKey(value))
            {
                throw new KeyNotFoundException($"DataSetAttributes does not contain \"{value}\"");
            }

            _activeVectorsName = value;
        }
    }

    /// <summary>
    /// Gets or sets the name of the active normals array.
    /// </summary>
    /// <exception cref="InvalidOperationException">
    /// Thrown when setting on <see cref="FieldAssociation.None"/> association.
    /// </exception>
    /// <exception cref="KeyNotFoundException">
    /// Thrown when the specified name does not exist in this attribute set.
    /// </exception>
    public string? ActiveNormalsName
    {
        get => _activeNormalsName is not null && _arrays.ContainsKey(_activeNormalsName)
            ? _activeNormalsName
            : null;
        set
        {
            if (value is null)
            {
                _activeNormalsName = null;
                return;
            }

            RaiseFieldDataNoScalarsVectorsNormals();
            if (!_arrays.ContainsKey(value))
            {
                throw new KeyNotFoundException($"DataSetAttributes does not contain \"{value}\"");
            }

            _activeNormalsName = value;
        }
    }

    /// <summary>
    /// Gets or sets the name of the active texture coordinates array.
    /// </summary>
    /// <exception cref="InvalidOperationException">
    /// Thrown when setting on <see cref="FieldAssociation.None"/> association.
    /// </exception>
    /// <exception cref="KeyNotFoundException">
    /// Thrown when the specified name does not exist in this attribute set.
    /// </exception>
    public string? ActiveTextureCoordinatesName
    {
        get => _activeTextureCoordinatesName is not null
            && _arrays.ContainsKey(_activeTextureCoordinatesName)
            ? _activeTextureCoordinatesName
            : null;
        set
        {
            if (value is null)
            {
                _activeTextureCoordinatesName = null;
                return;
            }

            RaiseNoTextureCoordinates();
            if (!_arrays.ContainsKey(value))
            {
                throw new KeyNotFoundException($"DataSetAttributes does not contain \"{value}\"");
            }

            _activeTextureCoordinatesName = value;
        }
    }

    /// <summary>
    /// Gets a collection of all array names in insertion order.
    /// </summary>
    public ICollection<string> Keys => _orderedKeys.AsReadOnly();

    /// <summary>
    /// Gets a collection of all arrays in insertion order.
    /// </summary>
    public ICollection<double[]> Values
    {
        get
        {
            var result = new List<double[]>(_orderedKeys.Count);
            foreach (var key in _orderedKeys)
            {
                result.Add(_arrays[key]);
            }

            return result.AsReadOnly();
        }
    }

    /// <summary>
    /// Gets or sets the array associated with the specified name.
    /// </summary>
    /// <param name="key">The name of the array.</param>
    /// <returns>The array associated with the specified name.</returns>
    /// <exception cref="ArgumentNullException">Thrown when <paramref name="key"/> is <c>null</c>.</exception>
    /// <exception cref="KeyNotFoundException">Thrown when the key does not exist (getter only).</exception>
    /// <remarks>
    /// When setting, if the key is new and there is no active scalar name,
    /// the new array is automatically made the active scalars for
    /// <see cref="FieldAssociation.Point"/> and <see cref="FieldAssociation.Cell"/> associations.
    /// </remarks>
    public double[] this[string key]
    {
        get
        {
            ArgumentNullException.ThrowIfNull(key);
            return GetArray(key);
        }
        set
        {
            ArgumentNullException.ThrowIfNull(key);
            ArgumentNullException.ThrowIfNull(value);

            bool existed = _arrays.ContainsKey(key);
            SetArray(value, key);

            if (existed)
            {
                return;
            }

            if (Association is FieldAssociation.Point or FieldAssociation.Cell
                && ActiveScalarsName is null)
            {
                _activeScalarsName = key;
            }
        }
    }

    /// <summary>
    /// Returns the value associated with the specified key, or a default value
    /// if the key does not exist.
    /// </summary>
    /// <param name="key">The name of the array to retrieve.</param>
    /// <param name="defaultValue">
    /// The value to return if <paramref name="key"/> is not found. Defaults to <c>null</c>.
    /// </param>
    /// <returns>The array if found; otherwise <paramref name="defaultValue"/>.</returns>
    public double[]? Get(string key, double[]? defaultValue = null)
    {
        ArgumentNullException.ThrowIfNull(key);
        return _arrays.TryGetValue(key, out var arr) ? arr : defaultValue;
    }

    /// <summary>
    /// Gets an array by name.
    /// </summary>
    /// <param name="key">The name of the array to retrieve.</param>
    /// <returns>The array associated with the specified name.</returns>
    /// <exception cref="KeyNotFoundException">Thrown when the key does not exist.</exception>
    public double[] GetArray(string key)
    {
        ArgumentNullException.ThrowIfNull(key);
        if (!_arrays.TryGetValue(key, out var arr))
        {
            throw new KeyNotFoundException(key);
        }

        return arr;
    }

    /// <summary>
    /// Gets an array by its insertion-order index.
    /// </summary>
    /// <param name="index">The zero-based index of the array.</param>
    /// <returns>The array at the specified index.</returns>
    /// <exception cref="ArgumentOutOfRangeException">
    /// Thrown when <paramref name="index"/> is out of range.
    /// </exception>
    public double[] GetArray(int index)
    {
        if (index < 0 || index >= _orderedKeys.Count)
        {
            throw new ArgumentOutOfRangeException(
                nameof(index),
                index,
                $"Array index ({index}) out of range [0, {_orderedKeys.Count - 1}]");
        }

        return _arrays[_orderedKeys[index]];
    }

    /// <summary>
    /// Adds an array to this object. If an array with the same name already exists,
    /// it is replaced.
    /// </summary>
    /// <param name="data">The array data.</param>
    /// <param name="name">The name to assign to the array.</param>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="data"/> or <paramref name="name"/> is <c>null</c>.
    /// </exception>
    public void SetArray(double[] data, string name)
    {
        ArgumentNullException.ThrowIfNull(data);
        ArgumentNullException.ThrowIfNull(name);

        if (_arrays.ContainsKey(name))
        {
            _arrays[name] = data;
        }
        else
        {
            _arrays[name] = data;
            _orderedKeys.Add(name);
        }
    }

    /// <summary>
    /// Sets the active scalars of the dataset with an array.
    /// </summary>
    /// <param name="scalars">The scalar data array.</param>
    /// <param name="name">The name to assign to the scalars. Defaults to <c>"scalars"</c>.</param>
    /// <exception cref="InvalidOperationException">
    /// Thrown when association is <see cref="FieldAssociation.None"/>.
    /// </exception>
    public void SetScalars(double[] scalars, string name = "scalars")
    {
        RaiseFieldDataNoScalarsVectorsNormals();
        SetArray(scalars, name);
        _activeScalarsName = name;
    }

    /// <summary>
    /// Sets the active vectors of this data attribute.
    /// Vectors are a quantity that has magnitude and direction.
    /// </summary>
    /// <param name="vectors">The vector data array (should have 3 components per element).</param>
    /// <param name="name">The name to assign to the vectors.</param>
    /// <exception cref="InvalidOperationException">
    /// Thrown when association is <see cref="FieldAssociation.None"/>.
    /// </exception>
    public void SetVectors(double[] vectors, string name)
    {
        RaiseFieldDataNoScalarsVectorsNormals();
        SetArray(vectors, name);
        _activeVectorsName = name;
    }

    /// <summary>
    /// Sets the active normals of this data attribute.
    /// </summary>
    /// <param name="normals">The normals data array (should have 3 components per element).</param>
    /// <param name="name">The name to assign to the normals. Defaults to <c>"Normals"</c>.</param>
    /// <exception cref="InvalidOperationException">
    /// Thrown when association is <see cref="FieldAssociation.None"/>.
    /// </exception>
    public void SetNormals(double[] normals, string name = "Normals")
    {
        RaiseFieldDataNoScalarsVectorsNormals();
        SetArray(normals, name);
        _activeNormalsName = name;
    }

    /// <summary>
    /// Removes the array with the specified name.
    /// </summary>
    /// <param name="key">The name of the array to remove.</param>
    /// <exception cref="KeyNotFoundException">Thrown when the key does not exist.</exception>
    public void Remove(string key)
    {
        ArgumentNullException.ThrowIfNull(key);
        if (!_arrays.ContainsKey(key))
        {
            throw new KeyNotFoundException($"{key} not present.");
        }

        _arrays.Remove(key);
        _orderedKeys.Remove(key);

        // Clear active references if they pointed to the removed array.
        if (_activeScalarsName == key) _activeScalarsName = null;
        if (_activeVectorsName == key) _activeVectorsName = null;
        if (_activeNormalsName == key) _activeNormalsName = null;
        if (_activeTextureCoordinatesName == key) _activeTextureCoordinatesName = null;
    }

    /// <summary>
    /// Removes the array with the specified name and returns it.
    /// </summary>
    /// <param name="key">The name of the array to remove and return.</param>
    /// <returns>The removed array.</returns>
    /// <exception cref="KeyNotFoundException">Thrown when the key does not exist.</exception>
    public double[] Pop(string key)
    {
        var arr = GetArray(key);
        Remove(key);
        return arr;
    }

    /// <summary>
    /// Removes the array with the specified name and returns it,
    /// or returns <paramref name="defaultValue"/> if the key is not found.
    /// </summary>
    /// <param name="key">The name of the array to remove and return.</param>
    /// <param name="defaultValue">The value to return if the key is not found.</param>
    /// <returns>The removed array, or <paramref name="defaultValue"/>.</returns>
    public double[]? Pop(string key, double[]? defaultValue)
    {
        ArgumentNullException.ThrowIfNull(key);
        if (!_arrays.ContainsKey(key))
        {
            return defaultValue;
        }

        return Pop(key);
    }

    /// <summary>
    /// Returns a list of (array name, array value) tuples in insertion order.
    /// </summary>
    /// <returns>A list of key-value pairs.</returns>
    public List<KeyValuePair<string, double[]>> Items()
    {
        var result = new List<KeyValuePair<string, double[]>>(_orderedKeys.Count);
        foreach (var key in _orderedKeys)
        {
            result.Add(new KeyValuePair<string, double[]>(key, _arrays[key]));
        }

        return result;
    }

    /// <summary>
    /// Removes all arrays from this object and clears all active array references.
    /// </summary>
    public void Clear()
    {
        _arrays.Clear();
        _orderedKeys.Clear();
        _activeScalarsName = null;
        _activeVectorsName = null;
        _activeNormalsName = null;
        _activeTextureCoordinatesName = null;
    }

    /// <summary>
    /// Updates this object with arrays from the given dictionary.
    /// Existing arrays with matching names are replaced.
    /// </summary>
    /// <param name="arrayDict">A dictionary of array name to data mappings.</param>
    public void Update(IDictionary<string, double[]> arrayDict)
    {
        ArgumentNullException.ThrowIfNull(arrayDict);
        foreach (var (name, array) in arrayDict)
        {
            this[name] = (double[])array.Clone();
        }
    }

    /// <summary>
    /// Determines whether an array with the specified name exists.
    /// </summary>
    /// <param name="key">The array name to locate.</param>
    /// <returns><c>true</c> if the array exists; otherwise, <c>false</c>.</returns>
    public bool ContainsKey(string key) => _arrays.ContainsKey(key);

    /// <summary>
    /// Adds an array. Throws if the name already exists.
    /// </summary>
    /// <param name="key">The name of the array.</param>
    /// <param name="value">The array data.</param>
    /// <exception cref="ArgumentException">
    /// Thrown when an array with the same name already exists.
    /// </exception>
    public void Add(string key, double[] value)
    {
        ArgumentNullException.ThrowIfNull(key);
        ArgumentNullException.ThrowIfNull(value);
        if (_arrays.ContainsKey(key))
        {
            throw new ArgumentException($"An array with name '{key}' already exists.", nameof(key));
        }

        this[key] = value;
    }

    /// <summary>
    /// Attempts to get the array associated with the specified name.
    /// </summary>
    /// <param name="key">The array name.</param>
    /// <param name="value">
    /// When this method returns, contains the array if found; otherwise, <c>null</c>.
    /// </param>
    /// <returns><c>true</c> if the array was found; otherwise, <c>false</c>.</returns>
    public bool TryGetValue(string key, [MaybeNullWhen(false)] out double[] value) =>
        _arrays.TryGetValue(key, out value);

    /// <inheritdoc />
    bool IDictionary<string, double[]>.Remove(string key)
    {
        if (!_arrays.ContainsKey(key))
        {
            return false;
        }

        Remove(key);
        return true;
    }

    /// <inheritdoc />
    void ICollection<KeyValuePair<string, double[]>>.Add(KeyValuePair<string, double[]> item) =>
        Add(item.Key, item.Value);

    /// <inheritdoc />
    bool ICollection<KeyValuePair<string, double[]>>.Contains(KeyValuePair<string, double[]> item) =>
        _arrays.TryGetValue(item.Key, out var val) && ReferenceEquals(val, item.Value);

    /// <inheritdoc />
    void ICollection<KeyValuePair<string, double[]>>.CopyTo(KeyValuePair<string, double[]>[] array, int arrayIndex)
    {
        ArgumentNullException.ThrowIfNull(array);
        if (arrayIndex < 0 || arrayIndex + _orderedKeys.Count > array.Length)
        {
            throw new ArgumentOutOfRangeException(nameof(arrayIndex));
        }

        foreach (var key in _orderedKeys)
        {
            array[arrayIndex++] = new KeyValuePair<string, double[]>(key, _arrays[key]);
        }
    }

    /// <inheritdoc />
    bool ICollection<KeyValuePair<string, double[]>>.Remove(KeyValuePair<string, double[]> item)
    {
        if (_arrays.TryGetValue(item.Key, out var val) && ReferenceEquals(val, item.Value))
        {
            Remove(item.Key);
            return true;
        }

        return false;
    }

    /// <summary>
    /// Returns an enumerator that iterates through the array names in insertion order.
    /// </summary>
    /// <returns>An enumerator of key-value pairs.</returns>
    public IEnumerator<KeyValuePair<string, double[]>> GetEnumerator()
    {
        foreach (var key in _orderedKeys)
        {
            yield return new KeyValuePair<string, double[]>(key, _arrays[key]);
        }
    }

    /// <inheritdoc />
    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

    /// <summary>
    /// Returns a string representation of this <see cref="DataSetAttributes"/> instance,
    /// including the association, active arrays, and contained array names.
    /// </summary>
    /// <returns>A human-readable summary.</returns>
    public override string ToString()
    {
        var lines = new List<string> { "PyVista DataSetAttributes" };
        lines.Add($"Association     : {Association}");

        if (Association is FieldAssociation.Point or FieldAssociation.Cell)
        {
            lines.Add($"Active Scalars  : {ActiveScalarsName ?? "None"}");
            lines.Add($"Active Vectors  : {ActiveVectorsName ?? "None"}");
            lines.Add($"Active Texture  : {ActiveTextureCoordinatesName ?? "None"}");
            lines.Add($"Active Normals  : {ActiveNormalsName ?? "None"}");
        }

        if (_orderedKeys.Count == 0)
        {
            lines.Add("Contains arrays : None");
        }
        else
        {
            lines.Add("Contains arrays :");
            foreach (var key in _orderedKeys)
            {
                var arr = _arrays[key];
                var displayName = key.Length > 23 ? string.Concat(key.AsSpan(0, 20), "...") : key;
                lines.Add($"    {displayName,-24}float64    ({arr.Length})");
            }
        }

        return string.Join(Environment.NewLine, lines);
    }

    private void RaiseFieldDataNoScalarsVectorsNormals()
    {
        if (Association == FieldAssociation.None)
        {
            throw new InvalidOperationException(
                "FieldData does not have active scalars or vectors or normals.");
        }
    }

    private void RaiseNoTextureCoordinates()
    {
        if (Association == FieldAssociation.None)
        {
            throw new InvalidOperationException(
                "FieldData does not have active texture coordinates.");
        }
    }
}
