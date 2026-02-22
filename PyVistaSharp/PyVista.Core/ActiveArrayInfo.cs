namespace PyVista.Core;

/// <summary>
/// Immutable record containing information about an active array,
/// including its field association and name.
/// </summary>
/// <param name="Association">The field association of the array.</param>
/// <param name="Name">The name of the array, or <c>null</c> if not set.</param>
public sealed record ActiveArrayInfoTuple(FieldAssociation Association, string? Name);

/// <summary>
/// Active array info class with support for serialization.
/// </summary>
public sealed class ActiveArrayInfo : IEquatable<ActiveArrayInfo>
{
    /// <summary>
    /// Initializes a new instance of the <see cref="ActiveArrayInfo"/> class.
    /// </summary>
    /// <param name="association">The field association of the array.</param>
    /// <param name="name">The name of the array, or <c>null</c> if not set.</param>
    public ActiveArrayInfo(FieldAssociation association, string? name)
    {
        Association = association;
        Name = name;
    }

    /// <summary>Gets the field association of the array.</summary>
    public FieldAssociation Association { get; }

    /// <summary>Gets the name of the array.</summary>
    public string? Name { get; }

    /// <summary>
    /// Returns a copy of this instance.
    /// </summary>
    /// <returns>A new <see cref="ActiveArrayInfo"/> with the same values.</returns>
    public ActiveArrayInfo Copy() => new(Association, Name);

    /// <summary>
    /// Converts this instance to an <see cref="ActiveArrayInfoTuple"/>.
    /// </summary>
    /// <returns>An equivalent <see cref="ActiveArrayInfoTuple"/>.</returns>
    public ActiveArrayInfoTuple ToTuple() => new(Association, Name);

    /// <inheritdoc />
    public bool Equals(ActiveArrayInfo? other)
    {
        if (other is null) return false;
        if (ReferenceEquals(this, other)) return true;
        return Association == other.Association && Name == other.Name;
    }

    /// <inheritdoc />
    public override bool Equals(object? obj) => Equals(obj as ActiveArrayInfo);

    /// <inheritdoc />
    public override int GetHashCode() => HashCode.Combine(Association, Name);

    /// <inheritdoc />
    public override string ToString() => $"ActiveArrayInfo(Association={Association}, Name={Name})";
}
