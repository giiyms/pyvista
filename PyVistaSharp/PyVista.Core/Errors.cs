// PyVista specific errors and warnings.

using System;

namespace PyVista.Core
{
    /// <summary>
    /// Exception when a mesh does not contain all triangles.
    /// </summary>
    public class NotAllTrianglesError : ArgumentException
    {
        private const string DefaultMessage = "Mesh must consist of only triangles";

        /// <summary>
        /// Initializes a new instance of the <see cref="NotAllTrianglesError"/> class.
        /// </summary>
        /// <param name="message">Error message.</param>
        public NotAllTrianglesError(string message = DefaultMessage)
            : base(message) { }

        /// <summary>
        /// Initializes a new instance of the <see cref="NotAllTrianglesError"/> class.
        /// </summary>
        /// <param name="message">Error message.</param>
        /// <param name="innerException">The inner exception.</param>
        public NotAllTrianglesError(string message, Exception innerException)
            : base(message, innerException) { }
    }

    /// <summary>
    /// Used for deprecated methods and functions.
    /// </summary>
    public class DeprecationError : InvalidOperationException
    {
        private const string DefaultMessage = "This feature has been deprecated";

        /// <summary>
        /// Initializes a new instance of the <see cref="DeprecationError"/> class.
        /// </summary>
        /// <param name="message">Error message.</param>
        public DeprecationError(string message = DefaultMessage)
            : base(message) { }

        /// <summary>
        /// Initializes a new instance of the <see cref="DeprecationError"/> class.
        /// </summary>
        /// <param name="message">Error message.</param>
        /// <param name="innerException">The inner exception.</param>
        public DeprecationError(string message, Exception innerException)
            : base(message, innerException) { }
    }

    /// <summary>
    /// Requested feature is not supported by the installed VTK version.
    /// </summary>
    public class VTKVersionError : InvalidOperationException
    {
        private const string DefaultMessage = "The requested feature is not supported by the installed VTK version.";

        /// <summary>
        /// Initializes a new instance of the <see cref="VTKVersionError"/> class.
        /// </summary>
        /// <param name="message">Error message.</param>
        public VTKVersionError(string message = DefaultMessage)
            : base(message) { }

        /// <summary>
        /// Initializes a new instance of the <see cref="VTKVersionError"/> class.
        /// </summary>
        /// <param name="message">Error message.</param>
        /// <param name="innerException">The inner exception.</param>
        public VTKVersionError(string message, Exception innerException)
            : base(message, innerException) { }
    }

    /// <summary>
    /// Requested filter or property is not supported by the PointSet class.
    /// </summary>
    public class PointSetNotSupported : InvalidCastException
    {
        private const string DefaultMessage = "The requested operation is not supported for PointSets.";

        /// <summary>
        /// Initializes a new instance of the <see cref="PointSetNotSupported"/> class.
        /// </summary>
        /// <param name="message">Error message.</param>
        public PointSetNotSupported(string message = DefaultMessage)
            : base(message) { }

        /// <summary>
        /// Initializes a new instance of the <see cref="PointSetNotSupported"/> class.
        /// </summary>
        /// <param name="message">Error message.</param>
        /// <param name="innerException">The inner exception.</param>
        public PointSetNotSupported(string message, Exception innerException)
            : base(message, innerException) { }
    }

    /// <summary>
    /// Cell operations are not supported on PointSets because they contain no cells.
    /// </summary>
    public class PointSetCellOperationError : PointSetNotSupported
    {
        private const string DefaultMessage = "Cell operations are not supported. PointSets contain no cells.";

        /// <summary>
        /// Initializes a new instance of the <see cref="PointSetCellOperationError"/> class.
        /// </summary>
        /// <param name="message">Error message.</param>
        public PointSetCellOperationError(string message = DefaultMessage)
            : base(message) { }

        /// <summary>
        /// Initializes a new instance of the <see cref="PointSetCellOperationError"/> class.
        /// </summary>
        /// <param name="message">Error message.</param>
        /// <param name="innerException">The inner exception.</param>
        public PointSetCellOperationError(string message, Exception innerException)
            : base(message, innerException) { }
    }

    /// <summary>
    /// Slice and other dimension reducing filters are not supported on PointSets.
    /// </summary>
    public class PointSetDimensionReductionError : PointSetNotSupported
    {
        private const string DefaultMessage = "Slice and other dimension reducing filters are not supported on PointSets.";

        /// <summary>
        /// Initializes a new instance of the <see cref="PointSetDimensionReductionError"/> class.
        /// </summary>
        /// <param name="message">Error message.</param>
        public PointSetDimensionReductionError(string message = DefaultMessage)
            : base(message) { }

        /// <summary>
        /// Initializes a new instance of the <see cref="PointSetDimensionReductionError"/> class.
        /// </summary>
        /// <param name="message">Error message.</param>
        /// <param name="innerException">The inner exception.</param>
        public PointSetDimensionReductionError(string message, Exception innerException)
            : base(message, innerException) { }
    }

    /// <summary>
    /// Requested filter or property is not supported by the PartitionedDataSets class.
    /// </summary>
    public class PartitionedDataSetsNotSupported : InvalidCastException
    {
        private const string DefaultMessage = "The requested operation is not supported for PartitionedDataSets.";

        /// <summary>
        /// Initializes a new instance of the <see cref="PartitionedDataSetsNotSupported"/> class.
        /// </summary>
        /// <param name="message">Error message.</param>
        public PartitionedDataSetsNotSupported(string message = DefaultMessage)
            : base(message) { }

        /// <summary>
        /// Initializes a new instance of the <see cref="PartitionedDataSetsNotSupported"/> class.
        /// </summary>
        /// <param name="message">Error message.</param>
        /// <param name="innerException">The inner exception.</param>
        public PartitionedDataSetsNotSupported(string message, Exception innerException)
            : base(message, innerException) { }
    }

    /// <summary>
    /// Exception when data is missing, e.g. no active scalars can be set.
    /// </summary>
    public class MissingDataError : ArgumentException
    {
        private const string DefaultMessage = "No data available.";

        /// <summary>
        /// Initializes a new instance of the <see cref="MissingDataError"/> class.
        /// </summary>
        /// <param name="message">Error message.</param>
        public MissingDataError(string message = DefaultMessage)
            : base(message) { }

        /// <summary>
        /// Initializes a new instance of the <see cref="MissingDataError"/> class.
        /// </summary>
        /// <param name="message">Error message.</param>
        /// <param name="innerException">The inner exception.</param>
        public MissingDataError(string message, Exception innerException)
            : base(message, innerException) { }
    }

    /// <summary>
    /// Exception when data is ambiguous, e.g. multiple active scalars can be set.
    /// </summary>
    public class AmbiguousDataError : ArgumentException
    {
        private const string DefaultMessage = "Multiple data available.";

        /// <summary>
        /// Initializes a new instance of the <see cref="AmbiguousDataError"/> class.
        /// </summary>
        /// <param name="message">Error message.</param>
        public AmbiguousDataError(string message = DefaultMessage)
            : base(message) { }

        /// <summary>
        /// Initializes a new instance of the <see cref="AmbiguousDataError"/> class.
        /// </summary>
        /// <param name="message">Error message.</param>
        /// <param name="innerException">The inner exception.</param>
        public AmbiguousDataError(string message, Exception innerException)
            : base(message, innerException) { }
    }

    /// <summary>
    /// Exception when a cell array size is invalid.
    /// </summary>
    public class CellSizeError : ArgumentException
    {
        private const string DefaultMessage = "Cell array size is invalid.";

        /// <summary>
        /// Initializes a new instance of the <see cref="CellSizeError"/> class.
        /// </summary>
        /// <param name="message">Error message.</param>
        public CellSizeError(string message = DefaultMessage)
            : base(message) { }

        /// <summary>
        /// Initializes a new instance of the <see cref="CellSizeError"/> class.
        /// </summary>
        /// <param name="message">Error message.</param>
        /// <param name="innerException">The inner exception.</param>
        public CellSizeError(string message, Exception innerException)
            : base(message, innerException) { }
    }

    /// <summary>
    /// Exception when a VTK pipeline runs into an issue.
    /// </summary>
    public class PyVistaPipelineError : InvalidOperationException
    {
        private const string DefaultMessage = "VTK pipeline issue detected by PyVista.";

        /// <summary>
        /// Initializes a new instance of the <see cref="PyVistaPipelineError"/> class.
        /// </summary>
        /// <param name="message">Error message.</param>
        public PyVistaPipelineError(string message = DefaultMessage)
            : base(message) { }

        /// <summary>
        /// Initializes a new instance of the <see cref="PyVistaPipelineError"/> class.
        /// </summary>
        /// <param name="message">Error message.</param>
        /// <param name="innerException">The inner exception.</param>
        public PyVistaPipelineError(string message, Exception innerException)
            : base(message, innerException) { }
    }

    /// <summary>
    /// Exception when accessing an attribute that is not part of the PyVista API.
    /// </summary>
    public class PyVistaAttributeError : MemberAccessException
    {
        private const string DefaultMessage = "The attribute is not part of the PyVista API.";

        /// <summary>
        /// Initializes a new instance of the <see cref="PyVistaAttributeError"/> class.
        /// </summary>
        /// <param name="message">Error message.</param>
        public PyVistaAttributeError(string message = DefaultMessage)
            : base(message) { }

        /// <summary>
        /// Initializes a new instance of the <see cref="PyVistaAttributeError"/> class.
        /// </summary>
        /// <param name="message">Error message.</param>
        /// <param name="innerException">The inner exception.</param>
        public PyVistaAttributeError(string message, Exception innerException)
            : base(message, innerException) { }
    }

    /// <summary>
    /// Error for invalid mesh properties.
    /// </summary>
    public class InvalidMeshError : ArgumentException
    {
        private const string DefaultMessage = "Invalid mesh.";

        /// <summary>
        /// Initializes a new instance of the <see cref="InvalidMeshError"/> class.
        /// </summary>
        /// <param name="message">Error message.</param>
        public InvalidMeshError(string message = DefaultMessage)
            : base(message) { }

        /// <summary>
        /// Initializes a new instance of the <see cref="InvalidMeshError"/> class.
        /// </summary>
        /// <param name="message">Error message.</param>
        /// <param name="innerException">The inner exception.</param>
        public InvalidMeshError(string message, Exception innerException)
            : base(message, innerException) { }
    }

    /// <summary>
    /// Exception when a VTK output message is detected.
    /// </summary>
    public class VTKExecutionError : InvalidOperationException
    {
        private const string DefaultMessage = "VTK output message was detected by PyVista.";

        /// <summary>
        /// Initializes a new instance of the <see cref="VTKExecutionError"/> class.
        /// </summary>
        /// <param name="message">Error message.</param>
        public VTKExecutionError(string message = DefaultMessage)
            : base(message) { }

        /// <summary>
        /// Initializes a new instance of the <see cref="VTKExecutionError"/> class.
        /// </summary>
        /// <param name="message">Error message.</param>
        /// <param name="innerException">The inner exception.</param>
        public VTKExecutionError(string message, Exception innerException)
            : base(message, innerException) { }
    }

    // --- Warning classes ---

    /// <summary>
    /// Base class for PyVista warnings.
    /// </summary>
    public class PyVistaWarning : Exception
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="PyVistaWarning"/> class.
        /// </summary>
        /// <param name="message">Warning message.</param>
        public PyVistaWarning(string? message = null)
            : base(message) { }

        /// <summary>
        /// Initializes a new instance of the <see cref="PyVistaWarning"/> class.
        /// </summary>
        /// <param name="message">Warning message.</param>
        /// <param name="innerException">The inner exception.</param>
        public PyVistaWarning(string message, Exception innerException)
            : base(message, innerException) { }
    }

    /// <summary>
    /// Non-suppressed deprecation warning.
    /// </summary>
    public class PyVistaDeprecationWarning : PyVistaWarning
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="PyVistaDeprecationWarning"/> class.
        /// </summary>
        /// <param name="message">Warning message.</param>
        public PyVistaDeprecationWarning(string? message = null)
            : base(message) { }

        /// <summary>
        /// Initializes a new instance of the <see cref="PyVistaDeprecationWarning"/> class.
        /// </summary>
        /// <param name="message">Warning message.</param>
        /// <param name="innerException">The inner exception.</param>
        public PyVistaDeprecationWarning(string message, Exception innerException)
            : base(message, innerException) { }
    }

    /// <summary>
    /// Non-suppressed future warning.
    /// </summary>
    public class PyVistaFutureWarning : PyVistaWarning
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="PyVistaFutureWarning"/> class.
        /// </summary>
        /// <param name="message">Warning message.</param>
        public PyVistaFutureWarning(string? message = null)
            : base(message) { }

        /// <summary>
        /// Initializes a new instance of the <see cref="PyVistaFutureWarning"/> class.
        /// </summary>
        /// <param name="message">Warning message.</param>
        /// <param name="innerException">The inner exception.</param>
        public PyVistaFutureWarning(string message, Exception innerException)
            : base(message, innerException) { }
    }

    /// <summary>
    /// Efficiency warning.
    /// </summary>
    public class PyVistaEfficiencyWarning : PyVistaWarning
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="PyVistaEfficiencyWarning"/> class.
        /// </summary>
        /// <param name="message">Warning message.</param>
        public PyVistaEfficiencyWarning(string? message = null)
            : base(message) { }

        /// <summary>
        /// Initializes a new instance of the <see cref="PyVistaEfficiencyWarning"/> class.
        /// </summary>
        /// <param name="message">Warning message.</param>
        /// <param name="innerException">The inner exception.</param>
        public PyVistaEfficiencyWarning(string message, Exception innerException)
            : base(message, innerException) { }
    }

    /// <summary>
    /// Warning when a VTK output message is detected.
    /// </summary>
    public class VTKExecutionWarning : PyVistaWarning
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="VTKExecutionWarning"/> class.
        /// </summary>
        /// <param name="message">Warning message.</param>
        public VTKExecutionWarning(string? message = null)
            : base(message) { }

        /// <summary>
        /// Initializes a new instance of the <see cref="VTKExecutionWarning"/> class.
        /// </summary>
        /// <param name="message">Warning message.</param>
        /// <param name="innerException">The inner exception.</param>
        public VTKExecutionWarning(string message, Exception innerException)
            : base(message, innerException) { }
    }

    /// <summary>
    /// Warning for invalid mesh properties.
    /// </summary>
    public class InvalidMeshWarning : PyVistaWarning
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="InvalidMeshWarning"/> class.
        /// </summary>
        /// <param name="message">Warning message.</param>
        public InvalidMeshWarning(string? message = null)
            : base(message) { }

        /// <summary>
        /// Initializes a new instance of the <see cref="InvalidMeshWarning"/> class.
        /// </summary>
        /// <param name="message">Warning message.</param>
        /// <param name="innerException">The inner exception.</param>
        public InvalidMeshWarning(string message, Exception innerException)
            : base(message, innerException) { }
    }
}
