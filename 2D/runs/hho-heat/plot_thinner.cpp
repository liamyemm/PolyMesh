#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <bits/stdc++.h>

int sign(double x0, double y0, double x1, double y1, double x2, double y2)
{
    double val = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);

    return (0 < val) - (val < 0);
}

int main(int argc, char *argv[])
{
    std::ifstream inFile;

    inFile.open(argv[1]);
    if (!inFile)
    {
        std::cerr << "Unable to open file\n";
        exit(1);
    }

    std::cout << "Reading data from file\n";

    std::array<double, 3> coord;
    std::vector<std::array<double, 3>> all_points;
    double xcoord , ycoord , zcoord;

    while (inFile >> xcoord >> ycoord >> zcoord)
    {
        coord[0] = xcoord;
        coord[1] = ycoord;
        coord[2] = zcoord;
        all_points.push_back(coord);
    }

    inFile.close();

    std::cout << "Sucessfully read " << all_points.size() << " points\n";

    std::vector<std::array<double, 3>> thinned_points;

    std::array<double, 3> point_on_hull = all_points[0];
    for (size_t i = 1; i < all_points.size(); ++i)
    {
        if (all_points[i][0] < point_on_hull[0])
        {
            point_on_hull = all_points[i];
        }
    }

    size_t count = 0;

    std::cout << "Forming convex hull\n";

    do
    {
        thinned_points.push_back(point_on_hull);
        std::array<double, 3> end_point = all_points[0];
        size_t index = 0;
        for (size_t j = 0; j < all_points.size(); ++j)
        {
            if( (end_point == all_points[j]) && (end_point != thinned_points[0]) )
            {
                continue;
            }
            if (sign(point_on_hull[0], point_on_hull[1], end_point[0], end_point[1], all_points[j][0], all_points[j][1]) > 0)
            {
                end_point = all_points[j];
                index = j;
            }
        }
        all_points.erase(all_points.begin() + index);
        point_on_hull = end_point;
        ++count;
    } while ( (point_on_hull != thinned_points[0]) && count < 5000);

    std::cout << "Outputting file\n";

    std::ofstream sol_out("solution_plot.tsv");
    for (size_t i = 0; i < thinned_points.size(); ++i)
    {
        sol_out << std::setprecision(12) << std::setw(20) << std::left << thinned_points[i][0] << std::setprecision(12) << std::setw(20) << std::left << thinned_points[i][1] << std::setprecision(12) << std::setw(20) << thinned_points[i][2] << std::endl;
    }
    sol_out.close();


    return 0;
}
