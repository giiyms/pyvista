using PyVista.Core;

namespace PyVista.Core.Utilities;

/// <summary>
/// Source class for generating binary images of an ellipsoid.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.ImageEllipsoidSource</c>.
/// Configure properties and call <see cref="GetOutput"/> to produce an <see cref="ImageData"/>.
/// </para>
/// </summary>
public class ImageEllipsoidSource
{
    /// <summary>Gets or sets the extent of the whole output image as <c>[x0, x1, y0, y1, z0, z1]</c>.</summary>
    public int[] WholeExtent { get; set; } = [0, 20, 0, 20, 0, 0];

    /// <summary>Gets or sets the center of the ellipsoid in <c>[x, y, z]</c>.</summary>
    public double[] Center { get; set; } = [10.0, 10.0, 0.0];

    /// <summary>Gets or sets the radii of the ellipsoid in <c>[rx, ry, rz]</c>.</summary>
    public double[] Radius { get; set; } = [3.0, 4.0, 5.0];

    /// <summary>
    /// Generates the ellipsoid image.
    /// </summary>
    /// <returns>An <see cref="ImageData"/> containing the binary ellipsoid image.</returns>
    public ImageData GetOutput()
    {
        int nx = WholeExtent[1] - WholeExtent[0] + 1;
        int ny = WholeExtent[3] - WholeExtent[2] + 1;
        int nz = Math.Max(WholeExtent[5] - WholeExtent[4] + 1, 1);

        var image = new ImageData();
        image.Dimensions = (nx, ny, nz);

        double rx2 = Radius[0] * Radius[0];
        double ry2 = Radius[1] * Radius[1];
        double rz2 = Radius[2] * Radius[2];

        var data = new double[nx * ny * nz];
        int idx = 0;
        for (int k = WholeExtent[4]; k <= WholeExtent[5] || (k == WholeExtent[4] && nz == 1); k++)
        {
            for (int j = WholeExtent[2]; j <= WholeExtent[3]; j++)
            {
                for (int i = WholeExtent[0]; i <= WholeExtent[1]; i++)
                {
                    double dx = i - Center[0];
                    double dy = j - Center[1];
                    double dz = nz > 1 ? k - Center[2] : 0;

                    double val = (dx * dx) / rx2 + (dy * dy) / ry2 + (rz2 > 0 ? (dz * dz) / rz2 : 0);
                    data[idx++] = val <= 1.0 ? 1.0 : 0.0;
                }
            }
            if (nz == 1) break;
        }

        image.PointData.SetArray(data, "ImageScalars");
        return image;
    }
}

/// <summary>
/// Source class for generating images of the Mandelbrot set.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.ImageMandelbrotSource</c>.
/// </para>
/// </summary>
public class ImageMandelbrotSource
{
    /// <summary>Gets or sets the extent of the whole output image as <c>[x0, x1, y0, y1, z0, z1]</c>.</summary>
    public int[] WholeExtent { get; set; } = [0, 200, 0, 200, 0, 0];

    /// <summary>Gets or sets the maximum number of iterations for computing the Mandelbrot set.</summary>
    public int MaxIter { get; set; } = 100;

    /// <summary>
    /// Generates the Mandelbrot set image.
    /// </summary>
    /// <returns>An <see cref="ImageData"/> containing the Mandelbrot set visualization.</returns>
    public ImageData GetOutput()
    {
        int nx = WholeExtent[1] - WholeExtent[0] + 1;
        int ny = WholeExtent[3] - WholeExtent[2] + 1;

        var image = new ImageData();
        image.Dimensions = (nx, ny, 1);

        var data = new double[nx * ny];

        double xMin = -2.0, xMax = 0.5;
        double yMin = -1.25, yMax = 1.25;

        int idx = 0;
        for (int j = 0; j < ny; j++)
        {
            double ci = yMin + (yMax - yMin) * j / (ny - 1);
            for (int i = 0; i < nx; i++)
            {
                double cr = xMin + (xMax - xMin) * i / (nx - 1);
                double zr = 0, zi = 0;
                int iter = 0;
                while (zr * zr + zi * zi <= 4.0 && iter < MaxIter)
                {
                    double tmp = zr * zr - zi * zi + cr;
                    zi = 2.0 * zr * zi + ci;
                    zr = tmp;
                    iter++;
                }
                data[idx++] = iter;
            }
        }

        image.PointData.SetArray(data, "Iterations");
        return image;
    }
}

/// <summary>
/// Source class for generating uniform random noise images.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.ImageNoiseSource</c>.
/// </para>
/// </summary>
public class ImageNoiseSource
{
    /// <summary>Gets or sets the extent of the whole output image as <c>[x0, x1, y0, y1, z0, z1]</c>.</summary>
    public int[] WholeExtent { get; set; } = [0, 200, 0, 200, 0, 0];

    /// <summary>Gets or sets the minimum noise value.</summary>
    public double Minimum { get; set; }

    /// <summary>Gets or sets the maximum noise value.</summary>
    public double Maximum { get; set; } = 1.0;

    /// <summary>Gets or sets the random seed for reproducibility.</summary>
    public int Seed { get; set; } = 1;

    /// <summary>
    /// Generates the noise image.
    /// </summary>
    /// <returns>An <see cref="ImageData"/> containing uniform random noise.</returns>
    public ImageData GetOutput()
    {
        int nx = WholeExtent[1] - WholeExtent[0] + 1;
        int ny = WholeExtent[3] - WholeExtent[2] + 1;
        int nz = Math.Max(WholeExtent[5] - WholeExtent[4] + 1, 1);

        var image = new ImageData();
        image.Dimensions = (nx, ny, nz);

        var rng = new Random(Seed);
        var data = new double[nx * ny * nz];
        double range = Maximum - Minimum;
        for (int i = 0; i < data.Length; i++)
        {
            data[i] = Minimum + rng.NextDouble() * range;
        }

        image.PointData.SetArray(data, "ImageScalars");
        return image;
    }
}

/// <summary>
/// Source class for generating sinusoidal images.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.ImageSinusoidSource</c>.
/// </para>
/// </summary>
public class ImageSinusoidSource
{
    /// <summary>Gets or sets the extent of the whole output image as <c>[x0, x1, y0, y1, z0, z1]</c>.</summary>
    public int[] WholeExtent { get; set; } = [0, 200, 0, 200, 0, 0];

    /// <summary>Gets or sets the direction of the sinusoid in <c>[x, y, z]</c>.</summary>
    public double[] Direction { get; set; } = [1.0, 0.0, 0.0];

    /// <summary>Gets or sets the period of the sinusoidal wave.</summary>
    public double Period { get; set; } = 20.0;

    /// <summary>Gets or sets the phase of the sinusoidal wave.</summary>
    public double Phase { get; set; }

    /// <summary>Gets or sets the amplitude of the sinusoidal wave.</summary>
    public double Amplitude { get; set; } = 255.0;

    /// <summary>
    /// Generates the sinusoidal image.
    /// </summary>
    /// <returns>An <see cref="ImageData"/> containing the sinusoidal pattern.</returns>
    public ImageData GetOutput()
    {
        int nx = WholeExtent[1] - WholeExtent[0] + 1;
        int ny = WholeExtent[3] - WholeExtent[2] + 1;
        int nz = Math.Max(WholeExtent[5] - WholeExtent[4] + 1, 1);

        var image = new ImageData();
        image.Dimensions = (nx, ny, nz);

        // Normalize direction
        double len = Math.Sqrt(Direction[0] * Direction[0] + Direction[1] * Direction[1] + Direction[2] * Direction[2]);
        double dx = len > 0 ? Direction[0] / len : 1;
        double dy = len > 0 ? Direction[1] / len : 0;
        double dz = len > 0 ? Direction[2] / len : 0;

        var data = new double[nx * ny * nz];
        int idx = 0;
        for (int k = 0; k < nz; k++)
        {
            for (int j = 0; j < ny; j++)
            {
                for (int i = 0; i < nx; i++)
                {
                    double dot = i * dx + j * dy + k * dz;
                    data[idx++] = Amplitude * Math.Sin(2.0 * Math.PI * dot / Period + Phase);
                }
            }
        }

        image.PointData.SetArray(data, "ImageScalars");
        return image;
    }
}

/// <summary>
/// Source class for generating Gaussian images.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.ImageGaussianSource</c>.
/// </para>
/// </summary>
public class ImageGaussianSource
{
    /// <summary>Gets or sets the center of the Gaussian in <c>[x, y, z]</c>.</summary>
    public double[] Center { get; set; } = [0.0, 0.0, 0.0];

    /// <summary>Gets or sets the extent of the whole output image as <c>[x0, x1, y0, y1, z0, z1]</c>.</summary>
    public int[] WholeExtent { get; set; } = [-30, 30, -30, 30, -30, 30];

    /// <summary>Gets or sets the maximum value of the Gaussian function.</summary>
    public double Maximum { get; set; } = 1.0;

    /// <summary>Gets or sets the standard deviation of the Gaussian function.</summary>
    public double StandardDeviation { get; set; } = 100.0;

    /// <summary>
    /// Generates the Gaussian image.
    /// </summary>
    /// <returns>An <see cref="ImageData"/> containing the Gaussian distribution image.</returns>
    public ImageData GetOutput()
    {
        int nx = WholeExtent[1] - WholeExtent[0] + 1;
        int ny = WholeExtent[3] - WholeExtent[2] + 1;
        int nz = Math.Max(WholeExtent[5] - WholeExtent[4] + 1, 1);

        var image = new ImageData();
        image.Dimensions = (nx, ny, nz);

        double sigma2 = StandardDeviation * StandardDeviation;
        var data = new double[nx * ny * nz];
        int idx = 0;

        for (int k = WholeExtent[4]; k <= WholeExtent[5]; k++)
        {
            for (int j = WholeExtent[2]; j <= WholeExtent[3]; j++)
            {
                for (int i = WholeExtent[0]; i <= WholeExtent[1]; i++)
                {
                    double dx = i - Center[0];
                    double dy = j - Center[1];
                    double dz = k - Center[2];
                    double r2 = dx * dx + dy * dy + dz * dz;
                    data[idx++] = Maximum * Math.Exp(-r2 / (2.0 * sigma2));
                }
            }
        }

        image.PointData.SetArray(data, "ImageScalars");
        return image;
    }
}

/// <summary>
/// Source class for generating grid pattern images.
/// <para>
/// This is the C# equivalent of the Python <c>pyvista.ImageGridSource</c>.
/// </para>
/// </summary>
public class ImageGridSource
{
    /// <summary>Gets or sets the origin of the grid in <c>[x, y, z]</c>.</summary>
    public double[] Origin { get; set; } = [0.0, 0.0, 0.0];

    /// <summary>Gets or sets the extent of the grid as <c>[x0, x1, y0, y1, z0, z1]</c>.</summary>
    public int[] Extent { get; set; } = [0, 10, 0, 10, 0, 0];

    /// <summary>Gets or sets the spacing of the grid lines in <c>[x, y, z]</c>.</summary>
    public double[] Spacing { get; set; } = [1.0, 1.0, 1.0];

    /// <summary>
    /// Generates the grid pattern image.
    /// </summary>
    /// <returns>An <see cref="ImageData"/> containing the grid pattern.</returns>
    public ImageData GetOutput()
    {
        int nx = Extent[1] - Extent[0] + 1;
        int ny = Extent[3] - Extent[2] + 1;
        int nz = Math.Max(Extent[5] - Extent[4] + 1, 1);

        var image = new ImageData();
        image.Dimensions = (nx, ny, nz);
        image.Origin = (Origin[0], Origin[1], Origin[2]);
        image.Spacing = (Spacing[0], Spacing[1], Spacing[2]);

        var data = new double[nx * ny * nz];
        int idx = 0;
        for (int k = 0; k < nz; k++)
        {
            for (int j = 0; j < ny; j++)
            {
                for (int i = 0; i < nx; i++)
                {
                    bool onGrid = (i % (int)Math.Max(Spacing[0], 1) == 0)
                               || (j % (int)Math.Max(Spacing[1], 1) == 0)
                               || (nz > 1 && k % (int)Math.Max(Spacing[2], 1) == 0);
                    data[idx++] = onGrid ? 1.0 : 0.0;
                }
            }
        }

        image.PointData.SetArray(data, "ImageScalars");
        return image;
    }
}
