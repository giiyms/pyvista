using System.Collections;
using System.Globalization;
using System.Text;

namespace PyVista.Core;

/// <summary>
/// A tabular data container for named column arrays.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.Table</c> class
/// (a wrapper for <c>vtkTable</c>). It stores row-oriented data as named
/// <see cref="double"/>[] columns and provides dictionary-like access.
/// </para>
/// </summary>
public class Table : DataObject, IEnumerable<double[]>
{
    private readonly DataSetAttributes _rowArrays = new(FieldAssociation.Row);

    // ---------------------------------------------------------------
    //  Constructors
    // ---------------------------------------------------------------

    /// <summary>
    /// Initializes a new, empty <see cref="Table"/>.
    /// </summary>
    public Table()
    {
    }

    /// <summary>
    /// Initializes a new <see cref="Table"/> from a 1-D or 2-D array.
    /// <para>
    /// A 1-D array is treated as a single column.
    /// A 2-D array with shape <c>(nRows, nColumns)</c> produces one column per second-dimension index.
    /// </para>
    /// </summary>
    /// <param name="arrays">The array data.</param>
    /// <param name="deep">When <c>true</c> (default), a deep copy of the data is stored.</param>
    public Table(double[] arrays, bool deep = true)
    {
        ArgumentNullException.ThrowIfNull(arrays);
        var data = deep ? (double[])arrays.Clone() : arrays;
        _rowArrays.SetArray(data, "Array 0");
    }

    /// <summary>
    /// Initializes a new <see cref="Table"/> from a 2-D array with shape <c>(nRows, nColumns)</c>.
    /// Each column of the input becomes a named array (<c>"Array 0"</c>, <c>"Array 1"</c>, …).
    /// </summary>
    /// <param name="arrays">The 2-D array data.</param>
    /// <param name="deep">When <c>true</c> (default), deep copies are stored.</param>
    public Table(double[,] arrays, bool deep = true)
    {
        ArgumentNullException.ThrowIfNull(arrays);
        int nRows = arrays.GetLength(0);
        int nCols = arrays.GetLength(1);
        for (int c = 0; c < nCols; c++)
        {
            var col = new double[nRows];
            for (int r = 0; r < nRows; r++)
            {
                col[r] = arrays[r, c];
            }

            _rowArrays.SetArray(col, $"Array {c}");
        }
    }

    /// <summary>
    /// Initializes a new <see cref="Table"/> from a dictionary of named column arrays.
    /// </summary>
    /// <param name="arrayDict">A mapping of column names to data arrays.</param>
    /// <param name="deep">When <c>true</c> (default), deep copies are stored.</param>
    public Table(IDictionary<string, double[]> arrayDict, bool deep = true)
    {
        ArgumentNullException.ThrowIfNull(arrayDict);
        foreach (var kvp in arrayDict)
        {
            var data = deep ? (double[])kvp.Value.Clone() : kvp.Value;
            _rowArrays.SetArray(data, kvp.Key);
        }
    }

    /// <summary>
    /// Initializes a new <see cref="Table"/> as a copy of another <see cref="Table"/>.
    /// </summary>
    /// <param name="other">The source table.</param>
    /// <param name="deep">When <c>true</c> (default), performs a deep copy.</param>
    public Table(Table other, bool deep = true)
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
    //  Row / Column counts
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets or sets the number of rows.
    /// <para>
    /// When setting, all existing column arrays are resized (truncated or zero-padded).
    /// </para>
    /// </summary>
    public int NRows
    {
        get
        {
            if (_rowArrays.Count == 0)
            {
                return 0;
            }

            // All columns should have the same length; return the first.
            foreach (var arr in _rowArrays.Values)
            {
                return arr.Length;
            }

            return 0;
        }
        set
        {
            if (value < 0)
            {
                throw new ArgumentOutOfRangeException(nameof(value), "NRows cannot be negative.");
            }

            var keys = new List<string>(_rowArrays.Keys);
            foreach (var key in keys)
            {
                var old = _rowArrays[key];
                if (old.Length == value)
                {
                    continue;
                }

                var resized = new double[value];
                Array.Copy(old, resized, Math.Min(old.Length, value));
                _rowArrays.SetArray(resized, key);
            }
        }
    }

    /// <summary>
    /// Gets the number of columns (arrays).
    /// </summary>
    public int NColumns => _rowArrays.Count;

    /// <summary>
    /// Gets the number of arrays. Alias for <see cref="NColumns"/>.
    /// </summary>
    public int NArrays => NColumns;

    /// <inheritdoc />
    public override bool IsEmpty => NRows == 0 && NColumns == 0;

    // ---------------------------------------------------------------
    //  Row arrays (DataSetAttributes view)
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets the row-associated data arrays.
    /// </summary>
    public DataSetAttributes RowArrays => _rowArrays;

    // ---------------------------------------------------------------
    //  Dictionary-like access
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets all column names.
    /// </summary>
    /// <returns>A collection of column names.</returns>
    public ICollection<string> Keys() => _rowArrays.Keys;

    /// <summary>
    /// Gets all column arrays as a list of (name, array) pairs.
    /// </summary>
    /// <returns>A list of key-value pairs.</returns>
    public List<KeyValuePair<string, double[]>> Items() => _rowArrays.Items();

    /// <summary>
    /// Gets all column arrays.
    /// </summary>
    /// <returns>A collection of arrays.</returns>
    public ICollection<double[]> Values() => _rowArrays.Values;

    /// <summary>
    /// Gets or sets the column array with the specified name.
    /// </summary>
    /// <param name="name">Column name.</param>
    /// <returns>The column data array.</returns>
    /// <exception cref="KeyNotFoundException">
    /// Thrown when the column name does not exist (getter only).
    /// </exception>
    public double[] this[string name]
    {
        get => _rowArrays.GetArray(name);
        set
        {
            ArgumentNullException.ThrowIfNull(name);
            _rowArrays[name] = value;
        }
    }

    /// <summary>
    /// Returns the column array with the specified name.
    /// </summary>
    /// <param name="name">Column name.</param>
    /// <returns>The column data array.</returns>
    public double[] Get(string name) => _rowArrays.GetArray(name);

    /// <summary>
    /// Updates this table with the entries from the given dictionary.
    /// Existing columns with matching names are replaced.
    /// </summary>
    /// <param name="data">A dictionary of column names to arrays.</param>
    public void Update(IDictionary<string, double[]> data)
    {
        _rowArrays.Update(data);
    }

    /// <summary>
    /// Updates this table from a 2-D array. Each column is named
    /// <c>"Array 0"</c>, <c>"Array 1"</c>, etc.
    /// </summary>
    /// <param name="arrays">A 2-D array of shape <c>(nRows, nColumns)</c>.</param>
    public void Update(double[,] arrays)
    {
        ArgumentNullException.ThrowIfNull(arrays);
        int nRows = arrays.GetLength(0);
        int nCols = arrays.GetLength(1);
        var dict = new Dictionary<string, double[]>(nCols);
        for (int c = 0; c < nCols; c++)
        {
            var col = new double[nRows];
            for (int r = 0; r < nRows; r++)
            {
                col[r] = arrays[r, c];
            }

            dict[$"Array {c}"] = col;
        }

        _rowArrays.Update(dict);
    }

    /// <summary>
    /// Removes the column with the specified name and returns its data.
    /// </summary>
    /// <param name="name">Column name.</param>
    /// <returns>The removed column data.</returns>
    /// <exception cref="KeyNotFoundException">
    /// Thrown when the column name does not exist.
    /// </exception>
    public double[] Pop(string name) => _rowArrays.Pop(name);

    /// <summary>
    /// Removes the column with the specified name.
    /// </summary>
    /// <param name="name">Column name.</param>
    /// <exception cref="KeyNotFoundException">
    /// Thrown when the column name does not exist.
    /// </exception>
    public void RemoveColumn(string name) => _rowArrays.Remove(name);

    // ---------------------------------------------------------------
    //  GetDataRange
    // ---------------------------------------------------------------

    /// <summary>
    /// Gets the minimum and maximum values of the named column array.
    /// </summary>
    /// <param name="name">
    /// Column name. When <c>null</c>, the first column is used.
    /// </param>
    /// <param name="preference">Unused; present for API compatibility.</param>
    /// <returns>A (Min, Max) tuple. Returns <c>(NaN, NaN)</c> if the array is empty.</returns>
    public override (double Min, double Max) GetDataRange(
        string? name = null,
        FieldAssociation preference = FieldAssociation.Row)
    {
        if (_rowArrays.Count == 0)
        {
            return (double.NaN, double.NaN);
        }

        double[] arr;
        if (name is null)
        {
            // Use the first column.
            var enumerator = _rowArrays.Values.GetEnumerator();
            if (!enumerator.MoveNext())
            {
                return (double.NaN, double.NaN);
            }

            arr = enumerator.Current;
        }
        else
        {
            arr = _rowArrays.GetArray(name);
        }

        if (arr.Length == 0)
        {
            return (double.NaN, double.NaN);
        }

        double min = double.PositiveInfinity;
        double max = double.NegativeInfinity;
        foreach (double v in arr)
        {
            if (double.IsNaN(v))
            {
                continue;
            }

            if (v < min) min = v;
            if (v > max) max = v;
        }

        if (double.IsPositiveInfinity(min))
        {
            return (double.NaN, double.NaN);
        }

        return (min, max);
    }

    // ---------------------------------------------------------------
    //  Copy
    // ---------------------------------------------------------------

    /// <summary>
    /// Returns a copy of this <see cref="Table"/>.
    /// </summary>
    /// <param name="deep">When <c>true</c> (default), performs a deep copy.</param>
    /// <returns>A new <see cref="Table"/>.</returns>
    public new Table Copy(bool deep = true)
    {
        return new Table(this, deep);
    }

    /// <inheritdoc />
    public override void ShallowCopy(DataObject source)
    {
        base.ShallowCopy(source);
        if (source is Table other)
        {
            _rowArrays.Clear();
            foreach (var kvp in other._rowArrays)
            {
                _rowArrays.SetArray(kvp.Value, kvp.Key);
            }
        }
    }

    /// <inheritdoc />
    public override void DeepCopy(DataObject source)
    {
        base.DeepCopy(source);
        if (source is Table other)
        {
            _rowArrays.Clear();
            foreach (var kvp in other._rowArrays)
            {
                _rowArrays.SetArray((double[])kvp.Value.Clone(), kvp.Key);
            }
        }
    }

    // ---------------------------------------------------------------
    //  IEnumerable<double[]>
    // ---------------------------------------------------------------

    /// <summary>
    /// Returns an enumerator that iterates through all column arrays.
    /// </summary>
    /// <returns>An enumerator of column arrays.</returns>
    public IEnumerator<double[]> GetEnumerator()
    {
        foreach (var kvp in _rowArrays)
        {
            yield return kvp.Value;
        }
    }

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
            ("N Rows", NRows.ToString()),
            ("N Columns", NColumns.ToString()),
        };
    }

    /// <summary>
    /// Returns a console-friendly representation of this <see cref="Table"/>.
    /// </summary>
    /// <returns>A formatted string.</returns>
    public override string ToString()
    {
        var sb = new StringBuilder();
        sb.AppendLine($"Table (0x{GetHashCode():X})");

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

        if (_rowArrays.Count > 0)
        {
            sb.AppendLine($"  {"N Arrays:".PadRight(maxLen)}{_rowArrays.Count}");
        }

        return sb.ToString().TrimEnd();
    }
}
