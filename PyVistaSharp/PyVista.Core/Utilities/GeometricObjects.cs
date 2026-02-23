using PyVista.Core;

namespace PyVista.Core.Utilities;

/// <summary>
/// Provides factory methods for creating common geometric objects.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.core.utilities.geometric_objects</c>
/// module. Each method returns a fully constructed <see cref="PolyData"/> instance representing
/// a geometric shape. Parameters mirror the Python API.
/// </para>
/// </summary>
public static class GeometricObjects
{
    /// <summary>
    /// Creates the surface of a capsule (cylinder with hemispherical caps).
    /// </summary>
    /// <param name="center">Center of the capsule in <c>[x, y, z]</c>.</param>
    /// <param name="direction">Direction the capsule points to in <c>[x, y, z]</c>.</param>
    /// <param name="radius">Radius of the capsule.</param>
    /// <param name="cylinderLength">Length of the cylindrical section.</param>
    /// <param name="resolution">Number of points on the circular face.</param>
    /// <returns>A <see cref="PolyData"/> representing the capsule surface.</returns>
    public static PolyData Capsule(
        double[] center = null!,
        double[] direction = null!,
        double radius = 0.5,
        double cylinderLength = 1.0,
        int resolution = 30)
    {
        center ??= [0.0, 0.0, 0.0];
        direction ??= [1.0, 0.0, 0.0];

        var source = new CylinderSource
        {
            Center = center,
            Direction = direction,
            Radius = radius,
            Height = cylinderLength,
            Resolution = resolution,
            Capping = true,
        };
        return source.GetOutput();
    }

    /// <summary>
    /// Creates the surface of a cylinder.
    /// </summary>
    /// <param name="center">Center of the cylinder in <c>[x, y, z]</c>.</param>
    /// <param name="direction">Direction the cylinder axis points to in <c>[x, y, z]</c>.</param>
    /// <param name="radius">Radius of the cylinder.</param>
    /// <param name="height">Height of the cylinder.</param>
    /// <param name="resolution">Number of points on the circular face.</param>
    /// <param name="capping">When <c>true</c>, cap the cylinder ends with polygons.</param>
    /// <returns>A <see cref="PolyData"/> representing the cylinder surface.</returns>
    public static PolyData Cylinder(
        double[] center = null!,
        double[] direction = null!,
        double radius = 0.5,
        double height = 1.0,
        int resolution = 100,
        bool capping = true)
    {
        center ??= [0.0, 0.0, 0.0];
        direction ??= [1.0, 0.0, 0.0];

        var source = new CylinderSource
        {
            Center = center,
            Direction = direction,
            Radius = radius,
            Height = height,
            Resolution = resolution,
            Capping = capping,
        };
        return source.GetOutput();
    }

    /// <summary>
    /// Creates an arrow mesh.
    /// </summary>
    /// <param name="start">Start location of the arrow in <c>[x, y, z]</c>.</param>
    /// <param name="direction">Direction the arrow points to in <c>[x, y, z]</c>.</param>
    /// <param name="tipLength">Length of the arrow tip as a fraction of total length.</param>
    /// <param name="tipRadius">Radius of the arrow tip.</param>
    /// <param name="tipResolution">Number of faces around the arrow tip.</param>
    /// <param name="shaftRadius">Radius of the arrow shaft.</param>
    /// <param name="shaftResolution">Number of faces around the arrow shaft.</param>
    /// <param name="scaleFactor">Scaling factor applied to the entire arrow.</param>
    /// <returns>A <see cref="PolyData"/> representing the arrow mesh.</returns>
    public static PolyData Arrow(
        double[] start = null!,
        double[] direction = null!,
        double tipLength = 0.25,
        double tipRadius = 0.1,
        int tipResolution = 20,
        double shaftRadius = 0.05,
        int shaftResolution = 20,
        double scaleFactor = 1.0)
    {
        start ??= [0.0, 0.0, 0.0];
        direction ??= [1.0, 0.0, 0.0];

        var source = new ArrowSource
        {
            TipLength = tipLength,
            TipRadius = tipRadius,
            TipResolution = tipResolution,
            ShaftRadius = shaftRadius,
            ShaftResolution = shaftResolution,
        };
        var mesh = source.GetOutput();

        if (Math.Abs(scaleFactor - 1.0) > 1e-12)
        {
            for (int i = 0; i < mesh.NPoints; i++)
            {
                mesh.Points[i * 3] *= scaleFactor;
                mesh.Points[i * 3 + 1] *= scaleFactor;
                mesh.Points[i * 3 + 2] *= scaleFactor;
            }
        }

        return mesh;
    }

    /// <summary>
    /// Creates a sphere mesh.
    /// </summary>
    /// <param name="radius">Radius of the sphere.</param>
    /// <param name="center">Center of the sphere in <c>[x, y, z]</c>.</param>
    /// <param name="thetaResolution">Number of points in the azimuthal direction.</param>
    /// <param name="phiResolution">Number of points in the polar direction.</param>
    /// <param name="startTheta">Starting azimuthal angle in degrees.</param>
    /// <param name="endTheta">Ending azimuthal angle in degrees.</param>
    /// <param name="startPhi">Starting polar angle in degrees.</param>
    /// <param name="endPhi">Ending polar angle in degrees.</param>
    /// <returns>A <see cref="PolyData"/> representing the sphere surface.</returns>
    public static PolyData Sphere(
        double radius = 0.5,
        double[] center = null!,
        int thetaResolution = 30,
        int phiResolution = 30,
        double startTheta = 0.0,
        double endTheta = 360.0,
        double startPhi = 0.0,
        double endPhi = 180.0)
    {
        center ??= [0.0, 0.0, 0.0];

        var source = new SphereSource
        {
            Radius = radius,
            Center = center,
            ThetaResolution = thetaResolution,
            PhiResolution = phiResolution,
            StartTheta = startTheta,
            EndTheta = endTheta,
            StartPhi = startPhi,
            EndPhi = endPhi,
        };
        return source.GetOutput();
    }

    /// <summary>
    /// Creates a plane mesh.
    /// </summary>
    /// <param name="center">Center of the plane in <c>[x, y, z]</c>.</param>
    /// <param name="direction">Normal direction of the plane in <c>[x, y, z]</c>.</param>
    /// <param name="xSize">Length of the plane in the x direction before rotation.</param>
    /// <param name="ySize">Length of the plane in the y direction before rotation.</param>
    /// <param name="xResolution">Number of points in the x direction.</param>
    /// <param name="yResolution">Number of points in the y direction.</param>
    /// <returns>A <see cref="PolyData"/> representing the plane.</returns>
    public static PolyData Plane(
        double[] center = null!,
        double[] direction = null!,
        double xSize = 1.0,
        double ySize = 1.0,
        int xResolution = 10,
        int yResolution = 10)
    {
        center ??= [0.0, 0.0, 0.0];
        direction ??= [0.0, 0.0, 1.0];

        var source = new PlaneSource
        {
            Center = center,
            Direction = direction,
            XSize = xSize,
            YSize = ySize,
            XResolution = xResolution,
            YResolution = yResolution,
        };
        return source.GetOutput();
    }

    /// <summary>
    /// Creates a line mesh between two points.
    /// </summary>
    /// <param name="pointa">Start point of the line in <c>[x, y, z]</c>.</param>
    /// <param name="pointb">End point of the line in <c>[x, y, z]</c>.</param>
    /// <param name="resolution">Number of pieces to divide the line into.</param>
    /// <returns>A <see cref="PolyData"/> representing the line.</returns>
    public static PolyData Line(
        double[] pointa = null!,
        double[] pointb = null!,
        int resolution = 1)
    {
        pointa ??= [0.0, 0.0, 0.0];
        pointb ??= [1.0, 0.0, 0.0];

        var source = new LineSource
        {
            PointA = pointa,
            PointB = pointb,
            Resolution = resolution,
        };
        return source.GetOutput();
    }

    /// <summary>
    /// Creates a box mesh defined by axis-aligned bounds.
    /// </summary>
    /// <param name="bounds">
    /// Bounding box in the form <c>[xMin, xMax, yMin, yMax, zMin, zMax]</c>.
    /// </param>
    /// <param name="level">Level of subdivision of the box faces.</param>
    /// <param name="quads">When <c>true</c>, faces are quads; otherwise triangles.</param>
    /// <returns>A <see cref="PolyData"/> representing the box surface.</returns>
    public static PolyData Box(
        double[] bounds = null!,
        int level = 0,
        bool quads = true)
    {
        bounds ??= [-1.0, 1.0, -1.0, 1.0, -1.0, 1.0];

        var source = new BoxSource
        {
            Bounds = bounds,
            Level = level,
            Quads = quads,
        };
        return source.GetOutput();
    }

    /// <summary>
    /// Creates a cone mesh.
    /// </summary>
    /// <param name="center">Center of the cone in <c>[x, y, z]</c>.</param>
    /// <param name="direction">Direction the cone tip points to in <c>[x, y, z]</c>.</param>
    /// <param name="height">Height of the cone.</param>
    /// <param name="radius">Base radius of the cone.</param>
    /// <param name="capping">When <c>true</c>, cap the base of the cone.</param>
    /// <param name="angle">Angle of the cone in degrees. Overrides <paramref name="radius"/> if set.</param>
    /// <param name="resolution">Number of points on the circular face.</param>
    /// <returns>A <see cref="PolyData"/> representing the cone surface.</returns>
    public static PolyData Cone(
        double[] center = null!,
        double[] direction = null!,
        double height = 1.0,
        double radius = 0.5,
        bool capping = true,
        double? angle = null,
        int resolution = 6)
    {
        center ??= [0.0, 0.0, 0.0];
        direction ??= [1.0, 0.0, 0.0];

        var source = new ConeSource
        {
            Center = center,
            Direction = direction,
            Height = height,
            Radius = radius,
            Capping = capping,
            Resolution = resolution,
        };
        if (angle.HasValue)
        {
            source.Angle = angle.Value;
        }
        return source.GetOutput();
    }

    /// <summary>
    /// Creates a regular polygon mesh.
    /// </summary>
    /// <param name="center">Center of the polygon in <c>[x, y, z]</c>.</param>
    /// <param name="radius">Radius of the polygon.</param>
    /// <param name="normal">Normal direction of the polygon face in <c>[x, y, z]</c>.</param>
    /// <param name="numberOfSides">Number of sides of the polygon.</param>
    /// <param name="fill">When <c>true</c>, create a filled polygon; otherwise just the edges.</param>
    /// <returns>A <see cref="PolyData"/> representing the polygon.</returns>
    public static PolyData Polygon(
        double[] center = null!,
        double radius = 1.0,
        double[] normal = null!,
        int numberOfSides = 6,
        bool fill = true)
    {
        center ??= [0.0, 0.0, 0.0];
        normal ??= [0.0, 0.0, 1.0];

        var source = new PolygonSource
        {
            Center = center,
            Radius = radius,
            Normal = normal,
            NumberOfSides = numberOfSides,
            Fill = fill,
        };
        return source.GetOutput();
    }

    /// <summary>
    /// Creates a disc (annulus) mesh.
    /// </summary>
    /// <param name="center">Center of the disc in <c>[x, y, z]</c>.</param>
    /// <param name="innerRadius">Inner radius of the disc.</param>
    /// <param name="outerRadius">Outer radius of the disc.</param>
    /// <param name="normal">Normal direction of the disc in <c>[x, y, z]</c>.</param>
    /// <param name="radialResolution">Number of points in the radial direction.</param>
    /// <param name="circumferentialResolution">Number of points in the circumferential direction.</param>
    /// <returns>A <see cref="PolyData"/> representing the disc.</returns>
    public static PolyData Disc(
        double[] center = null!,
        double innerRadius = 0.25,
        double outerRadius = 0.5,
        double[] normal = null!,
        int radialResolution = 1,
        int circumferentialResolution = 6)
    {
        center ??= [0.0, 0.0, 0.0];
        normal ??= [0.0, 0.0, 1.0];

        var source = new DiscSource
        {
            Center = center,
            InnerRadius = innerRadius,
            OuterRadius = outerRadius,
            Normal = normal,
            RadialResolution = radialResolution,
            CircumferentialResolution = circumferentialResolution,
        };
        return source.GetOutput();
    }

    /// <summary>
    /// Creates a cube mesh.
    /// </summary>
    /// <param name="center">Center of the cube in <c>[x, y, z]</c>.</param>
    /// <param name="xLength">Length of the cube in the x direction.</param>
    /// <param name="yLength">Length of the cube in the y direction.</param>
    /// <param name="zLength">Length of the cube in the z direction.</param>
    /// <param name="bounds">
    /// Optional explicit bounds <c>[xMin, xMax, yMin, yMax, zMin, zMax]</c>.
    /// When provided, <paramref name="center"/> and length parameters are ignored.
    /// </param>
    /// <returns>A <see cref="PolyData"/> representing the cube surface.</returns>
    public static PolyData Cube(
        double[] center = null!,
        double xLength = 1.0,
        double yLength = 1.0,
        double zLength = 1.0,
        double[]? bounds = null)
    {
        center ??= [0.0, 0.0, 0.0];

        var source = new CubeSource
        {
            Center = center,
            XLength = xLength,
            YLength = yLength,
            ZLength = zLength,
        };
        if (bounds != null)
        {
            source.Bounds = bounds;
        }
        return source.GetOutput();
    }

    /// <summary>
    /// Creates a circle (polygon) mesh.
    /// </summary>
    /// <param name="radius">Radius of the circle.</param>
    /// <param name="resolution">Number of points on the circle.</param>
    /// <returns>A <see cref="PolyData"/> representing the circle.</returns>
    public static PolyData Circle(double radius = 0.5, int resolution = 100)
    {
        return Polygon(
            center: [0.0, 0.0, 0.0],
            radius: radius,
            normal: [0.0, 0.0, 1.0],
            numberOfSides: resolution,
            fill: true);
    }

    /// <summary>
    /// Creates a superquadric mesh.
    /// </summary>
    /// <param name="center">Center of the superquadric in <c>[x, y, z]</c>.</param>
    /// <param name="scale">Scale factors in <c>[x, y, z]</c>.</param>
    /// <param name="size">Size of the superquadric.</param>
    /// <param name="thetaRoundness">Theta roundness parameter. Values less than 1 make the shape more square.</param>
    /// <param name="phiRoundness">Phi roundness parameter. Values less than 1 make the shape more square.</param>
    /// <param name="thetaResolution">Number of points in the theta direction.</param>
    /// <param name="phiResolution">Number of points in the phi direction.</param>
    /// <param name="toroidal">When <c>true</c>, creates a toroidal superquadric.</param>
    /// <param name="thickness">Thickness of the toroidal shape (only used when <paramref name="toroidal"/> is <c>true</c>).</param>
    /// <returns>A <see cref="PolyData"/> representing the superquadric surface.</returns>
    public static PolyData Superquadric(
        double[] center = null!,
        double[] scale = null!,
        double size = 0.5,
        double thetaRoundness = 1.0,
        double phiRoundness = 1.0,
        int thetaResolution = 16,
        int phiResolution = 16,
        bool toroidal = false,
        double thickness = 1.0 / 3.0)
    {
        center ??= [0.0, 0.0, 0.0];
        scale ??= [1.0, 1.0, 1.0];

        var source = new SuperquadricSource
        {
            Center = center,
            Scale = scale,
            Size = size,
            ThetaRoundness = thetaRoundness,
            PhiRoundness = phiRoundness,
            ThetaResolution = thetaResolution,
            PhiResolution = phiResolution,
            Toroidal = toroidal,
            Thickness = thickness,
        };
        return source.GetOutput();
    }

    /// <summary>
    /// Creates a platonic solid mesh.
    /// </summary>
    /// <param name="kind">
    /// Kind of platonic solid to create. One of <c>"tetrahedron"</c>, <c>"cube"</c>,
    /// <c>"octahedron"</c>, <c>"icosahedron"</c>, or <c>"dodecahedron"</c>.
    /// </param>
    /// <param name="radius">Radius of the circumscribing sphere.</param>
    /// <param name="center">Center of the solid in <c>[x, y, z]</c>.</param>
    /// <returns>A <see cref="PolyData"/> representing the platonic solid.</returns>
    /// <exception cref="ArgumentException">Thrown when <paramref name="kind"/> is not recognized.</exception>
    public static PolyData PlatonicSolid(
        string kind = "tetrahedron",
        double radius = 1.0,
        double[] center = null!)
    {
        center ??= [0.0, 0.0, 0.0];

        var source = new PlatonicSolidSource
        {
            Kind = kind,
            Radius = radius,
            Center = center,
        };
        return source.GetOutput();
    }

    /// <summary>
    /// Creates a tetrahedron mesh.
    /// </summary>
    /// <param name="radius">Radius of the circumscribing sphere.</param>
    /// <param name="center">Center of the tetrahedron in <c>[x, y, z]</c>.</param>
    /// <returns>A <see cref="PolyData"/> representing the tetrahedron.</returns>
    public static PolyData Tetrahedron(double radius = 1.0, double[] center = null!)
    {
        return PlatonicSolid("tetrahedron", radius, center ?? [0.0, 0.0, 0.0]);
    }

    /// <summary>
    /// Creates an octahedron mesh.
    /// </summary>
    /// <param name="radius">Radius of the circumscribing sphere.</param>
    /// <param name="center">Center of the octahedron in <c>[x, y, z]</c>.</param>
    /// <returns>A <see cref="PolyData"/> representing the octahedron.</returns>
    public static PolyData Octahedron(double radius = 1.0, double[] center = null!)
    {
        return PlatonicSolid("octahedron", radius, center ?? [0.0, 0.0, 0.0]);
    }

    /// <summary>
    /// Creates a dodecahedron mesh.
    /// </summary>
    /// <param name="radius">Radius of the circumscribing sphere.</param>
    /// <param name="center">Center of the dodecahedron in <c>[x, y, z]</c>.</param>
    /// <returns>A <see cref="PolyData"/> representing the dodecahedron.</returns>
    public static PolyData Dodecahedron(double radius = 1.0, double[] center = null!)
    {
        return PlatonicSolid("dodecahedron", radius, center ?? [0.0, 0.0, 0.0]);
    }

    /// <summary>
    /// Creates an icosahedron mesh.
    /// </summary>
    /// <param name="radius">Radius of the circumscribing sphere.</param>
    /// <param name="center">Center of the icosahedron in <c>[x, y, z]</c>.</param>
    /// <returns>A <see cref="PolyData"/> representing the icosahedron.</returns>
    public static PolyData Icosahedron(double radius = 1.0, double[] center = null!)
    {
        return PlatonicSolid("icosahedron", radius, center ?? [0.0, 0.0, 0.0]);
    }

    /// <summary>
    /// Creates a triangle mesh from three points.
    /// </summary>
    /// <param name="points">
    /// A 3×3 array of points defining the triangle vertices.
    /// When <c>null</c>, a default equilateral triangle is created.
    /// </param>
    /// <returns>A <see cref="PolyData"/> representing the triangle.</returns>
    /// <exception cref="ArgumentException">
    /// Thrown when <paramref name="points"/> does not contain exactly 9 values (3 points × 3 coordinates).
    /// </exception>
    public static PolyData Triangle(double[]? points = null)
    {
        points ??= [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        if (points.Length != 9)
        {
            throw new ArgumentException("Triangle requires exactly 3 points (9 values).", nameof(points));
        }

        var mesh = new PolyData();
        mesh.Points = points;
        return mesh;
    }

    /// <summary>
    /// Creates a rectangle mesh from four points.
    /// </summary>
    /// <param name="points">
    /// A 4×3 array (12 values) of corner points defining the rectangle.
    /// When <c>null</c>, a default unit rectangle in the xy-plane is created.
    /// </param>
    /// <returns>A <see cref="PolyData"/> representing the rectangle.</returns>
    /// <exception cref="ArgumentException">
    /// Thrown when <paramref name="points"/> does not contain exactly 12 values (4 points × 3 coordinates).
    /// </exception>
    public static PolyData Rectangle(double[]? points = null)
    {
        points ??= [1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0];
        if (points.Length != 12)
        {
            throw new ArgumentException("Rectangle requires exactly 4 points (12 values).", nameof(points));
        }

        var mesh = new PolyData();
        mesh.Points = points;
        return mesh;
    }

    /// <summary>
    /// Creates an ellipse mesh.
    /// </summary>
    /// <param name="semiMajorAxis">Length of the semi-major axis.</param>
    /// <param name="semiMinorAxis">Length of the semi-minor axis.</param>
    /// <param name="center">Center of the ellipse in <c>[x, y, z]</c>.</param>
    /// <param name="resolution">Number of points on the circumference.</param>
    /// <returns>A <see cref="PolyData"/> representing the ellipse.</returns>
    public static PolyData Ellipse(
        double semiMajorAxis = 0.5,
        double semiMinorAxis = 0.2,
        double[] center = null!,
        int resolution = 100)
    {
        center ??= [0.0, 0.0, 0.0];

        int numPoints = resolution + 1;
        var coords = new double[numPoints * 3];
        for (int i = 0; i < numPoints; i++)
        {
            double angle = 2.0 * Math.PI * i / resolution;
            coords[i * 3] = center[0] + semiMajorAxis * Math.Cos(angle);
            coords[i * 3 + 1] = center[1] + semiMinorAxis * Math.Sin(angle);
            coords[i * 3 + 2] = center[2];
        }

        var mesh = new PolyData();
        mesh.Points = coords;
        return mesh;
    }
}
