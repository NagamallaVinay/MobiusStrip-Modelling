#Mobius Strip is a one-sided (one Face) surface with only one edge.

# Import Necessary Packages to Perform Operations
import numpy as np  # For numerical calculations
import matplotlib.pyplot as plt  # For plotting
from mpl_toolkits.mplot3d import Axes3D  # Handles 3D visualization


# Create the MobiusStrip Class
class MobiusStrip:
    def __init__(self, R, w, n):
        """
        Initializes the Mobius strip.
        Parameters:
        R : Radius of the strip
        w : Width of the strip
        n : Resolution (number of points)
        """
        self.R = R  # Radius
        self.w = w  # width
        self.n = n  # resolution
        self.points = None

    ## Generate Points for Mobius Strip
    def generate_points(self):
        """Creates 3D coordinates for the Mobius strip."""
        u = np.linspace(0, 2*np.pi, self.n)  # Angle around the strip
        v = np.linspace(-self.w/2, self.w/2, self.n)  # Width variation

        u_grid, v_grid = np.meshgrid(u, v)  # Create a 2D grid

        # Applying the Mobius strip parametric equations
        x = (self.R + v_grid * np.cos(u_grid/2)) * np.cos(u_grid)
        y = (self.R + v_grid * np.cos(u_grid/2)) * np.sin(u_grid)
        z = v_grid * np.sin(u_grid/2)

        self.points = np.stack((x, y, z), axis=-1)  # Stores points

    ## Calculate Surface Area of Mobius Strip
    def calculate_surface_area(self):
        if self.points is None:
            self.generate_points()

        total_area = 0.0

        # Loop through the grid (except edges)
        for i in range(self.n - 1):
            for j in range(self.n - 1):
                # Grab four neighboring points
                p00 = self.points[i, j]
                p01 = self.points[i, j+1]
                p10 = self.points[i+1, j]
                p11 = self.points[i+1, j+1]

                # Divide the patch into two triangles
                v1 = p10 - p00
                v2 = p01 - p00
                area = np.linalg.norm(np.cross(v1, v2)) / 2

                v1 = p11 - p10
                v2 = p01 - p10
                area += np.linalg.norm(np.cross(v1, v2)) / 2

                total_area += area  # sum up all tiny triangle areas for total area

        return total_area

    ## Calculate Edge Length
    def calculate_edge_length(self):
        """Finds the total length of the edge."""
        if self.points is None:
            self.generate_points()

        # Get edge points (v = -w/2)
        edge_points = self.points[0, :, :]  # First row of points

        # Sum distances between consecutive edge points
        lengths = np.linalg.norm(edge_points[1:] - edge_points[:-1], axis=1)
        total_length = np.sum(lengths)

        return total_length * 2  # Both edges have the same length

    ## Plot the Mobius Strip
    def plot(self):
        """Mobius strip Representation in 3D."""
        if self.points is None:
            self.generate_points()

        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')

        x = self.points[..., 0]
        y = self.points[..., 1]
        z = self.points[..., 2]

        ax.plot_surface(x, y, z, color='royalblue', alpha=0.6) #alpha indicates transperancy
        ax.set_title('Mobius Strip Representation')
        plt.show()


## Run the Code by Providing Sample inputs
if __name__ == "__main__":
    R = float(input("Enter Radius :"))
    w = float(input("Enter Width :"))
    n = int(input("Enter Resolution :"))
    mobius = MobiusStrip(R, w, n)  # Create a Mobius strip
    mobius.generate_points()  # Generate 3D points
    print(f"Surface area: {mobius.calculate_surface_area():.2f}")
    print(f"Edge length: {mobius.calculate_edge_length():.2f}")
    mobius.plot()  # Display the Mobius strip



###  Code Structuring
"""
Initially we Created "MobiusStrip" and included methods for generating points,
 computing surface area, calculating edge length, and visualizing the strip.
"""

#### Challenges
"""
      =>The surface isn’t flat, so breaking it into triangles generates errors.
      =>While calculating Edge,Identifying the correct edge points without
           redundancy was difficult.
      =>Choosing the right color while Plotting Mobius Strip was tricky
          as we need to ensure clean visibility of representation
    """


### Surface Area Approximation
""" => As we know that Mobius strip is a continuous curved surface,
          we approximate its area by dividing it into tiny triangular segments.
      => Re-calling the concept of vector cross-products in Vector calculus,
           we find area of each triangle and sum them up to get surface area.
    """
