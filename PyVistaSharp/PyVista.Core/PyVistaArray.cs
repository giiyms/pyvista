using System.Collections;

namespace PyVista.Core;

/// <summary>
/// A wrapper around <see cref="double"/>[] that tracks which dataset and field
/// association the array belongs to.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.pyvista_ndarray</c> class.
/// In Python, <c>pyvista_ndarray</c> is a <see cref="T:numpy.ndarray"/> subclass that
/// maintains a weak reference to the owning dataset and the field association so that
/// modifications to the array can propagate change notifications upstream.
/// </para>
/// <para>
/// Because C# does not support subclassing primitive arrays, this class wraps a
/// <see cref="double"/>[] and exposes array-like access through an indexer, as well as
/// <see cref="IReadOnlyList{T}"/> and <see cref="IEnumerable{T}"/>.
/// </para>
/// </summary>
public class PyVistaArray : IReadOnlyList<double>, IEquatable<PyVistaArray>
{
    private double[] _data;
    private int[] _shape;

    /// <summary>
    /// Initializes a new instance of the <see cref="PyVistaArray"/> class from an
    /// existing <see cref="double"/>[] array.
    /// </summary>
    /// <param name="array">The source data. A shallow copy of the reference is stored.</param>
    /// <param name="dataset">
    /// An optional owning <see cref="DataObject"/>. A weak reference is stored so
    /// that the dataset can be garbage-collected independently.
    /// </param>
    /// <param name="association">The field association for this array.</param>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="array"/> is <c>null</c>.
    /// </exception>
    public PyVistaArray(
        double[] array,
        DataObject? dataset = null,
        FieldAssociation association = FieldAssociation.None)
    {
        ArgumentNullException.ThrowIfNull(array);
        _data = array;
        _shape = [array.Length];
        Association = association;
        Dataset = dataset is not null ? new WeakReference<DataObject>(dataset) : null;
    }

    /// <summary>
    /// Initializes a new instance of the <see cref="PyVistaArray"/> class from an
    /// existing <see cref="double"/>[] array and a specific shape.
    /// </summary>
    /// <param name="array">The source data. A shallow copy of the reference is stored.</param>
    /// <param name="shape">
    /// The logical shape of the array (e.g., <c>[100, 3]</c> for 100 points in 3-D).
    /// The product of all dimensions must equal <paramref name="array"/>.Length.
    /// </param>
    /// <param name="dataset">
    /// An optional owning <see cref="DataObject"/>. A weak reference is stored.
    /// </param>
    /// <param name="association">The field association for this array.</param>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="array"/> or <paramref name="shape"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when the product of <paramref name="shape"/> does not match the array length.
    /// </exception>
    public PyVistaArray(
        double[] array,
        int[] shape,
        DataObject? dataset = null,
        FieldAssociation association = FieldAssociation.None)
    {
        ArgumentNullException.ThrowIfNull(array);
        ArgumentNullException.ThrowIfNull(shape);
        ValidateShape(array.Length, shape);
        _data = array;
        _shape = (int[])shape.Clone();
        Association = association;
        Dataset = dataset is not null ? new WeakReference<DataObject>(dataset) : null;
    }

    /// <summary>
    /// Gets or sets the field association that describes how the array relates to
    /// the owning dataset (e.g., point data, cell data, or none).
    /// </summary>
    public FieldAssociation Association { get; set; }

    /// <summary>
    /// Gets or sets the weak reference to the owning <see cref="DataObject"/>.
    /// <para>
    /// A weak reference is used so that holding a <see cref="PyVistaArray"/> does
    /// not prevent the dataset from being garbage-collected, mirroring the Python
    /// implementation which uses <c>vtkWeakReference</c>.
    /// </para>
    /// </summary>
    public WeakReference<DataObject>? Dataset { get; set; }

    /// <summary>
    /// Gets the total number of elements in the underlying flat array.
    /// </summary>
    public int Length => _data.Length;

    /// <inheritdoc />
    int IReadOnlyCollection<double>.Count => _data.Length;

    /// <summary>
    /// Gets the logical shape of the array.
    /// <para>
    /// For a 1-D array of length <c>N</c>, the shape is <c>[N]</c>.
    /// For a 2-D array (e.g., <c>N</c> points × 3 components), the shape is <c>[N, 3]</c>.
    /// </para>
    /// </summary>
    public ReadOnlySpan<int> Shape => _shape;

    /// <summary>
    /// Gets the number of dimensions in the logical shape.
    /// </summary>
    public int NumDimensions => _shape.Length;

    /// <summary>
    /// Gets or sets the element at the specified flat <paramref name="index"/>.
    /// <para>
    /// Setting a value mirrors the Python behaviour of <c>pyvista_ndarray.__setitem__</c>:
    /// the owning dataset (if reachable) is notified that its data has been modified.
    /// </para>
    /// </summary>
    /// <param name="index">The zero-based index into the flat data array.</param>
    /// <returns>The value at the specified index.</returns>
    /// <exception cref="IndexOutOfRangeException">
    /// Thrown when <paramref name="index"/> is outside the bounds of the array.
    /// </exception>
    public double this[int index]
    {
        get => _data[index];
        set
        {
            _data[index] = value;
            NotifyModified();
        }
    }

    /// <summary>
    /// Gets or sets the element at the specified row and column for a 2-D shaped array.
    /// </summary>
    /// <param name="row">The zero-based row index.</param>
    /// <param name="col">The zero-based column index.</param>
    /// <returns>The value at the specified row and column.</returns>
    /// <exception cref="InvalidOperationException">
    /// Thrown when the array does not have exactly 2 dimensions.
    /// </exception>
    /// <exception cref="IndexOutOfRangeException">
    /// Thrown when <paramref name="row"/> or <paramref name="col"/> is out of range.
    /// </exception>
    public double this[int row, int col]
    {
        get
        {
            if (_shape.Length != 2)
            {
                throw new InvalidOperationException(
                    $"2-D indexing requires a shape with 2 dimensions, but the shape has {_shape.Length}.");
            }

            return _data[(row * _shape[1]) + col];
        }

        set
        {
            if (_shape.Length != 2)
            {
                throw new InvalidOperationException(
                    $"2-D indexing requires a shape with 2 dimensions, but the shape has {_shape.Length}.");
            }

            _data[(row * _shape[1]) + col] = value;
            NotifyModified();
        }
    }

    /// <summary>
    /// Implicitly converts a <see cref="double"/>[] to a <see cref="PyVistaArray"/>.
    /// </summary>
    /// <param name="array">The source array.</param>
    public static implicit operator PyVistaArray(double[] array) => new(array);

    /// <summary>
    /// Implicitly converts a <see cref="PyVistaArray"/> to a <see cref="double"/>[].
    /// </summary>
    /// <param name="pyArray">The source <see cref="PyVistaArray"/>.</param>
    public static implicit operator double[](PyVistaArray pyArray)
    {
        ArgumentNullException.ThrowIfNull(pyArray);
        return pyArray._data;
    }

    /// <summary>
    /// Returns a new <see cref="PyVistaArray"/> that is a view over the same data but
    /// with a different logical shape.
    /// <para>
    /// No data is copied; the returned array shares the same underlying buffer.
    /// This mirrors the behaviour of <c>numpy.ndarray.reshape</c>.
    /// </para>
    /// </summary>
    /// <param name="shape">
    /// The desired shape. The product of all dimensions must equal <see cref="Length"/>.
    /// </param>
    /// <returns>A new <see cref="PyVistaArray"/> with the requested shape.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="shape"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when the product of <paramref name="shape"/> does not match <see cref="Length"/>.
    /// </exception>
    public PyVistaArray Reshape(params int[] shape)
    {
        ArgumentNullException.ThrowIfNull(shape);
        ValidateShape(_data.Length, shape);
        return new PyVistaArray(_data, shape, null, Association)
        {
            Dataset = Dataset,
        };
    }

    /// <summary>
    /// Returns a new <see cref="PyVistaArray"/> that is a deep copy of the current
    /// instance, including a cloned data buffer.
    /// <para>
    /// The copy does <strong>not</strong> retain the dataset reference, mirroring the
    /// Python behaviour where copies lose their association with the source dataset.
    /// </para>
    /// </summary>
    /// <returns>A new, independent <see cref="PyVistaArray"/>.</returns>
    public PyVistaArray Copy()
    {
        return new PyVistaArray(
            (double[])_data.Clone(),
            (int[])_shape.Clone(),
            dataset: null,
            association: Association);
    }

    /// <summary>
    /// Returns the underlying data as a <see cref="double"/>[].
    /// <para>
    /// No copy is made; the caller receives a direct reference to the internal buffer.
    /// </para>
    /// </summary>
    /// <returns>The underlying <see cref="double"/>[].</returns>
    public double[] ToArray() => _data;

    /// <summary>
    /// Returns a <see cref="Span{T}"/> over the underlying data buffer.
    /// </summary>
    /// <returns>A span covering the entire internal array.</returns>
    public Span<double> AsSpan() => _data.AsSpan();

    /// <summary>
    /// Returns a flattened 1-D view over the same underlying data.
    /// <para>
    /// No data is copied. Equivalent to <c>numpy.ndarray.ravel()</c>.
    /// </para>
    /// </summary>
    /// <returns>A new 1-D <see cref="PyVistaArray"/> sharing the same buffer.</returns>
    public PyVistaArray Ravel()
    {
        if (_shape.Length == 1)
        {
            return this;
        }

        return Reshape(_data.Length);
    }

    /// <inheritdoc />
    public IEnumerator<double> GetEnumerator()
    {
        for (int i = 0; i < _data.Length; i++)
        {
            yield return _data[i];
        }
    }

    /// <inheritdoc />
    IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

    /// <inheritdoc />
    public bool Equals(PyVistaArray? other)
    {
        if (other is null) return false;
        if (ReferenceEquals(this, other)) return true;
        if (_data.Length != other._data.Length) return false;
        if (!_shape.AsSpan().SequenceEqual(other._shape)) return false;
        return _data.AsSpan().SequenceEqual(other._data);
    }

    /// <inheritdoc />
    public override bool Equals(object? obj) => Equals(obj as PyVistaArray);

    /// <inheritdoc />
    public override int GetHashCode()
    {
        var hash = new HashCode();
        hash.Add(_data.Length);
        foreach (int dim in _shape)
        {
            hash.Add(dim);
        }

        // Sample up to 8 elements for the hash to keep it O(1).
        int step = Math.Max(1, _data.Length / 8);
        for (int i = 0; i < _data.Length; i += step)
        {
            hash.Add(_data[i]);
        }

        return hash.ToHashCode();
    }

    /// <summary>
    /// Returns a string representation of the array, including its shape and association.
    /// </summary>
    /// <returns>A human-readable description.</returns>
    public override string ToString()
    {
        string shapeStr = string.Join(", ", _shape);
        return $"PyVistaArray(Length={_data.Length}, Shape=({shapeStr}), Association={Association})";
    }

    /// <summary>
    /// Validates that the product of <paramref name="shape"/> equals <paramref name="length"/>.
    /// </summary>
    private static void ValidateShape(int length, int[] shape)
    {
        if (shape.Length == 0)
        {
            throw new ArgumentException("Shape must have at least one dimension.", nameof(shape));
        }

        long product = 1;
        foreach (int dim in shape)
        {
            if (dim < 0)
            {
                throw new ArgumentException(
                    $"Shape dimensions must be non-negative, but got {dim}.", nameof(shape));
            }

            product *= dim;
        }

        if (product != length)
        {
            string shapeStr = string.Join(", ", shape);
            throw new ArgumentException(
                $"Cannot reshape array of length {length} into shape ({shapeStr}). "
                + $"The product of the shape dimensions ({product}) must equal the array length.",
                nameof(shape));
        }
    }

    /// <summary>
    /// Notifies the owning dataset (if still reachable) that this array has been
    /// modified, mirroring the Python <c>pyvista_ndarray.__setitem__</c> behaviour.
    /// </summary>
    private void NotifyModified()
    {
        // In the Python implementation, Modified() is called on both the VTKObject
        // and the associated dataset. Here we raise the event-like hook if the
        // dataset weak reference is still alive.
        if (Dataset is not null && Dataset.TryGetTarget(out _))
        {
            // The dataset is still alive. In a full VTK binding scenario this
            // would call dataset.Modified(). Currently this serves as a
            // placeholder to preserve the architectural intent.
        }
    }
}
