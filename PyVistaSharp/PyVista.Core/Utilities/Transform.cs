using System;
using System.Text;

namespace PyVista.Core.Utilities;

/// <summary>
/// Specifies the multiplication mode used when composing transformations.
/// </summary>
public enum MultiplyMode
{
    /// <summary>
    /// Additional transformations occur <em>before</em> the current matrix.
    /// </summary>
    Pre,

    /// <summary>
    /// Additional transformations occur <em>after</em> the current matrix (default).
    /// </summary>
    Post,
}

/// <summary>
/// Describes linear (affine) transformations via a 4×4 homogeneous matrix.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.Transform</c> class.
/// It supports translation, scaling, rotation, reflection, concatenation,
/// and inversion.  All operations use a right-handed coordinate system
/// with right-handed rotations.
/// </para>
/// <para>
/// Transformations can be chained fluently:
/// <code>
/// var t = new Transform()
///     .Translate(1, 2, 3)
///     .RotateZ(45)
///     .Scale(2);
/// </code>
/// </para>
/// </summary>
public class Transform
{
    private double[,] _matrix;
    private MultiplyMode _multiplyMode;

    /// <summary>
    /// Initializes a new instance of the <see cref="Transform"/> class
    /// set to the 4×4 identity matrix.
    /// </summary>
    public Transform()
    {
        _matrix = Transformations.Identity4x4();
        _multiplyMode = MultiplyMode.Post;
    }

    /// <summary>
    /// Initializes a new instance of the <see cref="Transform"/> class
    /// from an existing 4×4 matrix.
    /// </summary>
    /// <param name="matrix">A 4×4 transformation matrix to copy.</param>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="matrix"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when <paramref name="matrix"/> is not 4×4.
    /// </exception>
    public Transform(double[,] matrix)
    {
        ArgumentNullException.ThrowIfNull(matrix);
        ValidateMatrix(matrix);
        _matrix = (double[,])matrix.Clone();
        _multiplyMode = MultiplyMode.Post;
    }

    /// <summary>
    /// Gets or sets the 4×4 homogeneous transformation matrix.
    /// </summary>
    /// <exception cref="ArgumentException">
    /// Thrown when the assigned matrix is not 4×4.
    /// </exception>
    public double[,] Matrix
    {
        get => (double[,])_matrix.Clone();
        set
        {
            ArgumentNullException.ThrowIfNull(value);
            ValidateMatrix(value);
            _matrix = (double[,])value.Clone();
        }
    }

    /// <summary>
    /// Gets the inverse of the current transformation matrix.
    /// </summary>
    /// <returns>A 4×4 matrix that is the inverse of <see cref="Matrix"/>.</returns>
    /// <exception cref="InvalidOperationException">
    /// Thrown when the matrix is singular and cannot be inverted.
    /// </exception>
    public double[,] InverseMatrix => Transformations.InvertMatrix(Matrix);

    /// <summary>
    /// Gets or sets the multiplication mode used when composing transformations.
    /// </summary>
    public MultiplyMode MultiplyMode
    {
        get => _multiplyMode;
        set => _multiplyMode = value;
    }

    /// <summary>
    /// Sets the multiplication mode to <see cref="Utilities.MultiplyMode.Pre"/>.
    /// </summary>
    /// <returns>This <see cref="Transform"/> for fluent chaining.</returns>
    public Transform PreMultiply()
    {
        _multiplyMode = MultiplyMode.Pre;
        return this;
    }

    /// <summary>
    /// Sets the multiplication mode to <see cref="Utilities.MultiplyMode.Post"/>.
    /// </summary>
    /// <returns>This <see cref="Transform"/> for fluent chaining.</returns>
    public Transform PostMultiply()
    {
        _multiplyMode = MultiplyMode.Post;
        return this;
    }

    /// <summary>
    /// Composes a translation into the current transformation.
    /// </summary>
    /// <param name="x">Translation along the X axis.</param>
    /// <param name="y">Translation along the Y axis.</param>
    /// <param name="z">Translation along the Z axis.</param>
    /// <returns>This <see cref="Transform"/> for fluent chaining.</returns>
    public Transform Translate(double x, double y, double z)
    {
        var t = Transformations.Identity4x4();
        t[0, 3] = x;
        t[1, 3] = y;
        t[2, 3] = z;
        Concatenate(t);
        return this;
    }

    /// <summary>
    /// Composes a translation into the current transformation.
    /// </summary>
    /// <param name="vector">A 3-element translation vector (x, y, z).</param>
    /// <returns>This <see cref="Transform"/> for fluent chaining.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="vector"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when <paramref name="vector"/> does not have exactly 3 elements.
    /// </exception>
    public Transform Translate(double[] vector)
    {
        ArgumentNullException.ThrowIfNull(vector);
        if (vector.Length != 3)
        {
            throw new ArgumentException("Translation vector must have exactly 3 elements.", nameof(vector));
        }

        return Translate(vector[0], vector[1], vector[2]);
    }

    /// <summary>
    /// Composes a uniform or non-uniform scale into the current transformation.
    /// </summary>
    /// <param name="sx">Scale factor along the X axis.</param>
    /// <param name="sy">Scale factor along the Y axis.</param>
    /// <param name="sz">Scale factor along the Z axis.</param>
    /// <returns>This <see cref="Transform"/> for fluent chaining.</returns>
    public Transform Scale(double sx, double sy, double sz)
    {
        var s = Transformations.Identity4x4();
        s[0, 0] = sx;
        s[1, 1] = sy;
        s[2, 2] = sz;
        Concatenate(s);
        return this;
    }

    /// <summary>
    /// Composes a uniform scale into the current transformation.
    /// </summary>
    /// <param name="factor">Uniform scale factor applied to all axes.</param>
    /// <returns>This <see cref="Transform"/> for fluent chaining.</returns>
    public Transform Scale(double factor) => Scale(factor, factor, factor);

    /// <summary>
    /// Composes a rotation about the X axis.
    /// </summary>
    /// <param name="angleDegrees">Rotation angle in degrees (counterclockwise when
    /// looking from +X toward the origin).</param>
    /// <returns>This <see cref="Transform"/> for fluent chaining.</returns>
    public Transform RotateX(double angleDegrees)
    {
        var r = Transformations.AxisAngleRotation(new[] { 1.0, 0.0, 0.0 }, angleDegrees);
        Concatenate(r);
        return this;
    }

    /// <summary>
    /// Composes a rotation about the Y axis.
    /// </summary>
    /// <param name="angleDegrees">Rotation angle in degrees (counterclockwise when
    /// looking from +Y toward the origin).</param>
    /// <returns>This <see cref="Transform"/> for fluent chaining.</returns>
    public Transform RotateY(double angleDegrees)
    {
        var r = Transformations.AxisAngleRotation(new[] { 0.0, 1.0, 0.0 }, angleDegrees);
        Concatenate(r);
        return this;
    }

    /// <summary>
    /// Composes a rotation about the Z axis.
    /// </summary>
    /// <param name="angleDegrees">Rotation angle in degrees (counterclockwise when
    /// looking from +Z toward the origin).</param>
    /// <returns>This <see cref="Transform"/> for fluent chaining.</returns>
    public Transform RotateZ(double angleDegrees)
    {
        var r = Transformations.AxisAngleRotation(new[] { 0.0, 0.0, 1.0 }, angleDegrees);
        Concatenate(r);
        return this;
    }

    /// <summary>
    /// Composes a rotation about an arbitrary axis passing through the origin.
    /// </summary>
    /// <param name="axis">The direction vector of the rotation axis.</param>
    /// <param name="angleDegrees">Rotation angle in degrees.</param>
    /// <returns>This <see cref="Transform"/> for fluent chaining.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="axis"/> is <c>null</c>.
    /// </exception>
    public Transform Rotate(double[] axis, double angleDegrees)
    {
        ArgumentNullException.ThrowIfNull(axis);
        var r = Transformations.AxisAngleRotation(axis, angleDegrees);
        Concatenate(r);
        return this;
    }

    /// <summary>
    /// Composes a rotation about an arbitrary axis passing through the given point.
    /// </summary>
    /// <param name="axis">The direction vector of the rotation axis.</param>
    /// <param name="angleDegrees">Rotation angle in degrees.</param>
    /// <param name="point">A point on the rotation axis.</param>
    /// <returns>This <see cref="Transform"/> for fluent chaining.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="axis"/> or <paramref name="point"/> is <c>null</c>.
    /// </exception>
    public Transform Rotate(double[] axis, double angleDegrees, double[] point)
    {
        ArgumentNullException.ThrowIfNull(axis);
        ArgumentNullException.ThrowIfNull(point);
        var r = Transformations.AxisAngleRotation(axis, angleDegrees, point);
        Concatenate(r);
        return this;
    }

    /// <summary>
    /// Composes a reflection across the plane defined by the given normal through the origin.
    /// </summary>
    /// <param name="normal">The normal vector of the reflection plane.</param>
    /// <returns>This <see cref="Transform"/> for fluent chaining.</returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="normal"/> is <c>null</c>.
    /// </exception>
    public Transform Reflect(double[] normal)
    {
        ArgumentNullException.ThrowIfNull(normal);
        var r = Transformations.ReflectionMatrix(normal);
        Concatenate(r);
        return this;
    }

    /// <summary>
    /// Composes a reflection across the plane defined by the given normal through a point.
    /// </summary>
    /// <param name="normal">The normal vector of the reflection plane.</param>
    /// <param name="point">A point on the reflection plane.</param>
    /// <returns>This <see cref="Transform"/> for fluent chaining.</returns>
    public Transform Reflect(double[] normal, double[] point)
    {
        ArgumentNullException.ThrowIfNull(normal);
        ArgumentNullException.ThrowIfNull(point);
        var r = Transformations.ReflectionMatrix(normal, point);
        Concatenate(r);
        return this;
    }

    /// <summary>
    /// Concatenates the specified 4×4 matrix into the current transformation,
    /// respecting the current <see cref="MultiplyMode"/>.
    /// </summary>
    /// <param name="other">A 4×4 transformation matrix to compose.</param>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="other"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when <paramref name="other"/> is not 4×4.
    /// </exception>
    public void Concatenate(double[,] other)
    {
        ArgumentNullException.ThrowIfNull(other);
        ValidateMatrix(other);

        _matrix = _multiplyMode == MultiplyMode.Post
            ? Transformations.MultiplyMatrices(other, _matrix)
            : Transformations.MultiplyMatrices(_matrix, other);
    }

    /// <summary>
    /// Concatenates another <see cref="Transform"/> into this one.
    /// </summary>
    /// <param name="other">The transform to compose.</param>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="other"/> is <c>null</c>.
    /// </exception>
    public void Concatenate(Transform other)
    {
        ArgumentNullException.ThrowIfNull(other);
        Concatenate(other._matrix);
    }

    /// <summary>
    /// Replaces the current transformation with its inverse.
    /// </summary>
    /// <returns>This <see cref="Transform"/> for fluent chaining.</returns>
    /// <exception cref="InvalidOperationException">
    /// Thrown when the matrix is singular and cannot be inverted.
    /// </exception>
    public Transform Invert()
    {
        _matrix = Transformations.InvertMatrix(_matrix);
        return this;
    }

    /// <summary>
    /// Resets this transform to the identity matrix.
    /// </summary>
    /// <returns>This <see cref="Transform"/> for fluent chaining.</returns>
    public Transform Identity()
    {
        _matrix = Transformations.Identity4x4();
        return this;
    }

    /// <summary>
    /// Returns a deep copy of this transform.
    /// </summary>
    /// <returns>A new <see cref="Transform"/> with an independent copy of the matrix.</returns>
    public Transform Copy()
    {
        return new Transform((double[,])_matrix.Clone())
        {
            _multiplyMode = _multiplyMode,
        };
    }

    /// <summary>
    /// Applies this transformation to a set of 3D points.
    /// </summary>
    /// <param name="points">
    /// An array of shape [N, 3] stored as a flat array of length <c>N * 3</c>.
    /// </param>
    /// <returns>
    /// A new flat array of length <c>N * 3</c> containing the transformed points.
    /// </returns>
    /// <exception cref="ArgumentNullException">
    /// Thrown when <paramref name="points"/> is <c>null</c>.
    /// </exception>
    /// <exception cref="ArgumentException">
    /// Thrown when the length of <paramref name="points"/> is not a multiple of 3.
    /// </exception>
    public double[] ApplyToPoints(double[] points)
    {
        ArgumentNullException.ThrowIfNull(points);
        return Transformations.ApplyTransformationToPoints(_matrix, points);
    }

    /// <inheritdoc />
    public override string ToString()
    {
        var sb = new StringBuilder();
        sb.AppendLine($"Transform (0x{GetHashCode():X})");
        sb.AppendLine($"  MultiplyMode: {_multiplyMode}");
        sb.Append("  Matrix: ");
        for (int i = 0; i < 4; i++)
        {
            if (i > 0) sb.Append("          ");
            sb.Append('[');
            for (int j = 0; j < 4; j++)
            {
                if (j > 0) sb.Append(", ");
                sb.Append(_matrix[i, j].ToString("F4"));
            }
            sb.AppendLine("]");
        }

        return sb.ToString().TrimEnd();
    }

    private static void ValidateMatrix(double[,] matrix)
    {
        if (matrix.GetLength(0) != 4 || matrix.GetLength(1) != 4)
        {
            throw new ArgumentException(
                $"Matrix must be 4x4, but was {matrix.GetLength(0)}x{matrix.GetLength(1)}.",
                nameof(matrix));
        }
    }
}
