#include "misc.h"
#include "Box.h"
#include "Angular.h"
#include "Laser.h"
#include "BSpline.h"
#include "Atom.h"

#include <chrono>

int main()
{   

    nlohmann::json J = loadJson("input.json");
    std::array<int,3> initial_state{1,1,0};
    int rank = 0;

    Box box{J};
    Angular angular{J};
    Laser laser{J};
    Atom atom{J};
    BSpline bspline{J,box};


    std::string filename = "test.txt";
    std::ofstream outFile(filename);

    for (int spline{0}; spline < bspline.NBSpline(); ++spline)
    {
        for (int ridx{0}; ridx < box.Nr(); ++ridx)
        {   

            auto val1 = bspline.dB(spline,bspline.ecs_x(ridx * box.GridSpacing()));
            auto val2 = bspline.dBTest(spline,bspline.ecs_x(ridx * box.GridSpacing()));
            auto diff = std::abs(val1 - val2);
            outFile << diff << " ";
        }
    }
    outFile.close();



    angular.buildMaps(laser,initial_state);



    bspline.dumpTo(box,"misc",rank);
    angular.dumpTo("misc",rank);
    laser.dumpTo("misc",rank);

  

    return 0;
}