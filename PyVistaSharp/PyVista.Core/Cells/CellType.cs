// Define types of cells.
// Ported from pyvista/core/celltype.py and VTK's vtkCellType.h.

namespace PyVista.Core.Cells;

/// <summary>
/// Metadata for a cell type, including documentation and override information.
/// </summary>
/// <param name="Doc">Short description of this cell type.</param>
/// <param name="PointsOverride">
/// Override the displayed number of points. Use <c>"variable"</c> for composite cells
/// whose point count is not fixed (e.g. <see cref="CellType.PolyLine"/>), or <c>"n/a"</c>
/// when the concept does not apply.
/// </param>
/// <param name="EdgesOverride">
/// Override the displayed number of edges. Use <c>"variable"</c> for composite cells
/// whose edge count is not fixed (e.g. <see cref="CellType.TriangleStrip"/>), or
/// <c>"n/a"</c> when the concept does not apply.
/// </param>
/// <param name="FacesOverride">
/// Override the displayed number of faces. Use <c>"variable"</c> for cells whose face
/// count is not fixed (e.g. <see cref="CellType.Polyhedron"/>), or <c>"n/a"</c> when
/// the concept does not apply.
/// </param>
public readonly record struct CellTypeInfo(
    string? Doc = null,
    string? PointsOverride = null,
    string? EdgesOverride = null,
    string? FacesOverride = null);

/// <summary>
/// Defines types of cells.
/// <para>
/// Cells are defined by specifying a type in combination with an ordered list of points.
/// The ordered list, often referred to as the connectivity list, combined with the type
/// specification, implicitly defines the topology of the cell. The x-y-z point coordinates
/// define the cell geometry.
/// </para>
/// <para>
/// Although point coordinates are defined in three dimensions, the cell topology can be
/// 0, 1, 2, or 3-dimensional.
/// </para>
/// <para>
/// Cells can be primary (e.g. triangle) or composite (e.g. triangle strip). Composite cells
/// consist of one or more primary cells, while primary cells cannot be decomposed.
/// </para>
/// <para>
/// Cells can also be characterized as linear or non-linear. Linear cells use linear or constant
/// interpolation. Non-linear cells may use quadratic, cubic, or some other interpolation.
/// </para>
/// <para>
/// This enumeration defines all cell types used in VTK and supported by PyVista. The values
/// correspond to the VTK cell type constants defined in <c>vtkCellType.h</c>.
/// </para>
/// </summary>
/// <seealso href="https://book.vtk.org/en/latest/VTKBook/05Chapter5.html#cell-types">
/// VTK Book: Cell Types
/// </seealso>
/// <seealso href="https://vtk.org/doc/nightly/html/vtkCellType_8h_source.html">
/// vtkCellType.h — List of all cell types defined in VTK.
/// </seealso>
public enum CellType : int
{
    // ---------------------------------------------------------------------------------
    // Linear cells
    // ---------------------------------------------------------------------------------

    /// <summary>
    /// Used as a place-holder during processing.
    /// </summary>
    EmptyCell = 0,

    /// <summary>
    /// Represents a point in 3D space. The vertex is a primary zero-dimensional cell
    /// defined by a single point.
    /// </summary>
    Vertex = 1,

    /// <summary>
    /// Represents a set of points in 3D space. The polyvertex is a composite
    /// zero-dimensional cell defined by an arbitrarily ordered list of points.
    /// </summary>
    PolyVertex = 2,

    /// <summary>
    /// Represents a 1D line. The line is a primary one-dimensional cell defined by two
    /// points. The direction along the line is from the first point to the second point.
    /// </summary>
    Line = 3,

    /// <summary>
    /// Represents a set of 1D lines. The polyline is a composite one-dimensional cell
    /// consisting of one or more connected lines.
    /// </summary>
    PolyLine = 4,

    /// <summary>
    /// Represents a 2D triangle. The triangle is a primary two-dimensional cell defined
    /// by a counter-clockwise ordered list of three points.
    /// </summary>
    Triangle = 5,

    /// <summary>
    /// Represents a 2D triangle strip. The triangle strip is a composite two-dimensional
    /// cell consisting of one or more triangles. It is a compact representation of
    /// triangles connected edge-to-edge.
    /// </summary>
    TriangleStrip = 6,

    /// <summary>
    /// Represents a 2D n-sided polygon. The polygon is a primary two-dimensional cell
    /// defined by an ordered list of three or more points lying in a plane.
    /// </summary>
    Polygon = 7,

    /// <summary>
    /// Represents a 2D orthogonal quadrilateral. The pixel is a primary two-dimensional
    /// cell defined by an ordered list of four points.
    /// </summary>
    Pixel = 8,

    /// <summary>
    /// Represents a 2D quadrilateral. The quadrilateral is a primary two-dimensional cell
    /// defined by an ordered list of four points lying in a plane.
    /// </summary>
    Quad = 9,

    /// <summary>
    /// Represents a 3D tetrahedron. The tetrahedron is a primary three-dimensional cell
    /// defined by a list of four non-planar points. It has six edges and four triangular
    /// faces.
    /// </summary>
    Tetra = 10,

    /// <summary>
    /// Represents a 3D orthogonal parallelepiped. The voxel is a primary three-dimensional
    /// cell defined by an ordered list of eight points.
    /// </summary>
    Voxel = 11,

    /// <summary>
    /// Represents a 3D rectangular hexahedron. The hexahedron is a primary
    /// three-dimensional cell consisting of six quadrilateral faces, twelve edges,
    /// and eight vertices.
    /// </summary>
    Hexahedron = 12,

    /// <summary>
    /// Represents a linear 3D wedge. The wedge is a primary three-dimensional cell
    /// consisting of two triangular and three quadrilateral faces.
    /// </summary>
    Wedge = 13,

    /// <summary>
    /// Represents a 3D pyramid. The pyramid is a primary three-dimensional cell consisting
    /// of a rectangular base with four triangular faces, defined by an ordered list of
    /// five points.
    /// </summary>
    Pyramid = 14,

    /// <summary>
    /// Represents a convex 3D prism with a pentagonal base and five quadrilateral faces.
    /// The pentagonal prism is a primary three-dimensional cell defined by an ordered list
    /// of ten points.
    /// </summary>
    PentagonalPrism = 15,

    /// <summary>
    /// Represents a 3D prism with hexagonal base and six quadrilateral faces. The hexagonal
    /// prism is a primary three-dimensional cell defined by an ordered list of twelve points.
    /// </summary>
    HexagonalPrism = 16,

    // ---------------------------------------------------------------------------------
    // Quadratic, isoparametric cells
    // ---------------------------------------------------------------------------------

    /// <summary>
    /// Represents a 1D, 3-node, iso-parametric parabolic line. The cell includes a
    /// mid-edge node.
    /// </summary>
    QuadraticEdge = 21,

    /// <summary>
    /// Represents a 2D, 6-node, iso-parametric parabolic triangle. The cell includes a
    /// mid-edge node for each of the three edges of the cell.
    /// </summary>
    QuadraticTriangle = 22,

    /// <summary>
    /// Represents a 2D, 8-node iso-parametric parabolic quadrilateral element. The cell
    /// includes a mid-edge node for each of the four edges of the cell.
    /// </summary>
    QuadraticQuad = 23,

    /// <summary>
    /// Represents a 3D, 10-node, iso-parametric parabolic tetrahedron. The cell includes
    /// a mid-edge node on each of the side edges of the tetrahedron.
    /// </summary>
    QuadraticTetra = 24,

    /// <summary>
    /// Represents a 3D, 20-node iso-parametric parabolic hexahedron. The cell includes a
    /// mid-edge node.
    /// </summary>
    QuadraticHexahedron = 25,

    /// <summary>
    /// Represents a 3D, 15-node iso-parametric parabolic wedge. The cell includes a
    /// mid-edge node.
    /// </summary>
    QuadraticWedge = 26,

    /// <summary>
    /// Represents a 3D, 13-node iso-parametric parabolic pyramid. The cell includes a
    /// mid-edge node.
    /// </summary>
    QuadraticPyramid = 27,

    /// <summary>
    /// Represents a 2D, 9-node iso-parametric parabolic quadrilateral element with a
    /// center-point. The cell includes a mid-edge node for each of the four edges and
    /// a center node at the surface.
    /// </summary>
    BiquadraticQuad = 28,

    /// <summary>
    /// Represents a 3D, 27-node iso-parametric triquadratic hexahedron. The cell includes
    /// 8 edge nodes, 12 mid-edge nodes, 6 mid-face nodes and one mid-volume node.
    /// </summary>
    TriquadraticHexahedron = 29,

    /// <summary>
    /// Represents a 2D, 6-node iso-parametric quadratic-linear quadrilateral element.
    /// The cell includes a mid-edge node for two of the four edges.
    /// </summary>
    QuadraticLinearQuad = 30,

    /// <summary>
    /// Represents a 3D, 12-node iso-parametric linear quadratic wedge. The cell includes
    /// mid-edge nodes in the triangle edges.
    /// </summary>
    QuadraticLinearWedge = 31,

    /// <summary>
    /// Represents a 3D, 18-node iso-parametric bi-quadratic wedge. The cell includes a
    /// mid-edge node.
    /// </summary>
    BiquadraticQuadraticWedge = 32,

    /// <summary>
    /// Represents a 3D, 24-node iso-parametric biquadratic hexahedron. The cell includes
    /// mid-edge and center-face nodes.
    /// </summary>
    BiquadraticQuadraticHexahedron = 33,

    /// <summary>
    /// Represents a 2D, 7-node, iso-parametric parabolic triangle. The cell includes three
    /// mid-edge nodes besides the three triangle vertices and a center node.
    /// </summary>
    BiquadraticTriangle = 34,

    // ---------------------------------------------------------------------------------
    // Cubic, isoparametric cell
    // ---------------------------------------------------------------------------------

    /// <summary>
    /// Represents a 1D iso-parametric cubic line. The cell includes two mid-edge nodes.
    /// </summary>
    CubicLine = 35,

    // ---------------------------------------------------------------------------------
    // Quadratic polygon
    // ---------------------------------------------------------------------------------

    /// <summary>
    /// Represents a 2D n-sided (2*n nodes) parabolic polygon. The polygon cannot have any
    /// internal holes and cannot self-intersect. The cell includes a mid-edge node for each
    /// of the n edges of the cell.
    /// </summary>
    QuadraticPolygon = 36,

    /// <summary>
    /// Represents a second order 3D iso-parametric 19-node pyramid. The cell includes 5
    /// corner nodes, 8 mid-edge nodes, 5 mid-face nodes, and 1 volumetric centroid node.
    /// </summary>
    TriquadraticPyramid = 37,

    // ---------------------------------------------------------------------------------
    // Special class of cells formed by convex group of points
    // ---------------------------------------------------------------------------------

    /// <summary>
    /// Represents a cell formed by a convex group of points.
    /// </summary>
    ConvexPointSet = 41,

    // ---------------------------------------------------------------------------------
    // Polyhedron cell (consisting of polygonal faces)
    // ---------------------------------------------------------------------------------

    /// <summary>
    /// Represents a 3D cell defined by a set of polygonal faces.
    /// </summary>
    Polyhedron = 42,

    // ---------------------------------------------------------------------------------
    // Higher order cells in parametric form
    // ---------------------------------------------------------------------------------

    /// <summary>
    /// Higher order parametric curve.
    /// </summary>
    ParametricCurve = 51,

    /// <summary>
    /// Higher order parametric surface.
    /// </summary>
    ParametricSurface = 52,

    /// <summary>
    /// Higher order parametric triangular surface.
    /// </summary>
    ParametricTriSurface = 53,

    /// <summary>
    /// Higher order parametric quadrilateral surface.
    /// </summary>
    ParametricQuadSurface = 54,

    /// <summary>
    /// Higher order parametric tetrahedral region.
    /// </summary>
    ParametricTetraRegion = 55,

    /// <summary>
    /// Higher order parametric hexahedral region.
    /// </summary>
    ParametricHexRegion = 56,

    // ---------------------------------------------------------------------------------
    // Higher order cells
    // ---------------------------------------------------------------------------------

    /// <summary>
    /// Higher order edge.
    /// </summary>
    HigherOrderEdge = 60,

    /// <summary>
    /// Higher order triangle.
    /// </summary>
    HigherOrderTriangle = 61,

    /// <summary>
    /// Higher order quadrilateral.
    /// </summary>
    HigherOrderQuad = 62,

    /// <summary>
    /// Higher order polygon.
    /// </summary>
    HigherOrderPolygon = 63,

    /// <summary>
    /// Higher order tetrahedron.
    /// </summary>
    HigherOrderTetrahedron = 64,

    /// <summary>
    /// Higher order wedge.
    /// </summary>
    HigherOrderWedge = 65,

    /// <summary>
    /// Higher order pyramid.
    /// </summary>
    HigherOrderPyramid = 66,

    /// <summary>
    /// Higher order hexahedron.
    /// </summary>
    HigherOrderHexahedron = 67,

    // ---------------------------------------------------------------------------------
    // Arbitrary order Lagrange elements
    // ---------------------------------------------------------------------------------

    /// <summary>
    /// Arbitrary order Lagrange curve.
    /// </summary>
    LagrangeCurve = 68,

    /// <summary>
    /// Arbitrary order Lagrange triangle.
    /// </summary>
    LagrangeTriangle = 69,

    /// <summary>
    /// Arbitrary order Lagrange quadrilateral.
    /// </summary>
    LagrangeQuadrilateral = 70,

    /// <summary>
    /// Arbitrary order Lagrange tetrahedron.
    /// </summary>
    LagrangeTetrahedron = 71,

    /// <summary>
    /// Arbitrary order Lagrange hexahedron.
    /// </summary>
    LagrangeHexahedron = 72,

    /// <summary>
    /// Arbitrary order Lagrange wedge.
    /// </summary>
    LagrangeWedge = 73,

    /// <summary>
    /// Arbitrary order Lagrange pyramid.
    /// </summary>
    LagrangePyramid = 74,

    // ---------------------------------------------------------------------------------
    // Arbitrary order Bezier elements
    // ---------------------------------------------------------------------------------

    /// <summary>
    /// Arbitrary order Bezier curve.
    /// </summary>
    BezierCurve = 75,

    /// <summary>
    /// Arbitrary order Bezier triangle.
    /// </summary>
    BezierTriangle = 76,

    /// <summary>
    /// Arbitrary order Bezier quadrilateral.
    /// </summary>
    BezierQuadrilateral = 77,

    /// <summary>
    /// Arbitrary order Bezier tetrahedron.
    /// </summary>
    BezierTetrahedron = 78,

    /// <summary>
    /// Arbitrary order Bezier hexahedron.
    /// </summary>
    BezierHexahedron = 79,

    /// <summary>
    /// Arbitrary order Bezier wedge.
    /// </summary>
    BezierWedge = 80,

    /// <summary>
    /// Arbitrary order Bezier pyramid.
    /// </summary>
    BezierPyramid = 81,
}
